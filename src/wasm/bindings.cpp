// Thin embind layer (D6): exposes dmw_core to JS. Geometry logic lives in the core; this file
// only moves data across the boundary and owns the results so JS can view them.
#include <emscripten/bind.h>
#include <emscripten/val.h>

#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <sstream>
#include <string>
#include <vector>

#include "core/curvature.h"
#include "core/cusps.h"
#include "core/margin.h"
#include "core/bvh.h"
#include "core/undercut.h"
#include "core/generate.h"
#include "core/geodesic.h"
#include "core/halfedge.h"
#include "core/handles.h"
#include "core/io.h"
#include "core/topology.h"
#include "core/weld.h"

namespace {

using namespace dmw;
using emscripten::val;

// A JS typed array aliasing WASM memory: no copy. Valid only until the vector changes or
// WASM memory grows (which replaces the underlying ArrayBuffer), so JS copies it at once.
template <typename T>
val view(const std::vector<T>& v) {
    return val(emscripten::typed_memory_view(v.size(), v.data()));
}

std::string json_escape(const std::string& s) {
    std::string out;
    for (char ch : s) {
        if (ch == '"' || ch == '\\') out += '\\';
        if (static_cast<unsigned char>(ch) >= 0x20) out += ch;
    }
    return out;
}

class Session {
public:
    // `ptr` is the offset of `len` bytes JS already copied into WASM memory (malloc +
    // HEAPU8.set). Returns "" on success, else the loader's error message.
    std::string load(std::uintptr_t ptr, std::size_t len, const std::string& format, double weld_epsilon) {
        const auto* data = reinterpret_cast<const std::uint8_t*>(ptr);
        LoadResult r;
        if (format == "stl") {
            r = parse_stl(std::span<const std::uint8_t>(data, len));
            if (r.ok()) r.mesh = weld_vertices(r.mesh, weld_epsilon);  // STL stores a soup
        } else if (format == "obj") {
            r = parse_obj({reinterpret_cast<const char*>(data), len});
        } else {
            return "unknown format '" + format + "'";
        }
        if (!r.ok()) return r.error;
        set_mesh(std::move(r.mesh));
        return {};
    }

    std::string generate(const std::string& name) {
        TriMesh m;
        if (name == "icosphere") m = make_icosphere(4);
        else if (name == "torus") m = make_torus(72, 24);
        else if (name == "cylinder") m = make_cylinder(48, 15, 0.5, 1.0);
        else if (name == "grid") m = make_grid(16, 16);
        else if (name == "grid_holes") m = make_grid_with_holes(12, 12, {{3, 3}, {8, 6}});
        else if (name == "plate_handle") m = make_plate_with_handle(16, 24, 0.35);
        else if (name == "mobius") m = make_mobius(48);
        else if (name == "defects") m = make_defect_showcase();
        else if (name == "brick") m = make_brick_grid(30, 0.45);
        else return "unknown preset '" + name + "'";
        set_mesh(std::move(m));
        return {};
    }

    val positions() const { return view(positions_); }   // float32 xyz per vertex
    val indices() const { return view(indices_); }       // uint32, 3 per face (all faces)
    val mean() const { return view(mean_); }             // float32 per vertex, NaN = no data
    val gaussian() const { return view(gaussian_); }     // float32 per vertex, NaN = no data
    val kmin() const { return view(kmin_); }             // float32 min principal curvature (creases < 0)
    val faceComponent() const { return view(face_component_); }  // uint32, kInvalid = excluded
    val boundaryEdges() const { return view(boundary_edges_); }  // uint32 vertex pairs
    val nonmanifoldEdges() const { return view(nonmanifold_edges_); }
    val misorientedEdges() const { return view(misoriented_edges_); }
    val nonmanifoldVertices() const { return view(topo_.nonmanifold_vertices); }
    val handleEdges() const { return view(handle_edges_); }

    // Cusp tips with detector C at the operating point chosen on the training set (D68).
    std::string detectCusps() {
        if (!he_mesh_) return "cusp detection needs a manifold analysis view (" + he_error_ + ")";
        cusps_ = detect_cusps(*he_mesh_, cusp_operating_point()).vertices;
        return {};
    }
    val cusps() const { return view(cusps_); }  // uint32 vertex ids, strongest first

    // Tooth-gingiva margin at the training-chosen operating point (D73).
    std::string detectMargin() {
        if (!he_mesh_) return "margin detection needs a manifold analysis view (" + he_error_ + ")";
        const MarginResult r = detect_margin(*he_mesh_, margin_operating_point());
        predicted_tooth_ = r.tooth;
        margin_edges_.clear();
        for (const Edge& e : r.margin) margin_edges_.insert(margin_edges_.end(), {e.v0, e.v1});
        return {};
    }
    val marginEdges() const { return view(margin_edges_); }  // uint32 vertex pairs

    // Compare the predicted margin with ground-truth tooth flags (1 byte per vertex, 1 = tooth) that JS
    // copied into WASM memory at `ptr`. Same metric code as tools/margin_eval (D72). Returns JSON.
    std::string compareMargin(std::uintptr_t ptr, std::size_t len) {
        if (!he_mesh_ || predicted_tooth_.empty()) return "{\"error\":\"run detectMargin first\"}";
        if (len != mesh_.positions.size()) return "{\"error\":\"label count does not match vertex count\"}";
        const std::span<const std::uint8_t> truth(reinterpret_cast<const std::uint8_t*>(ptr), len);
        const auto gt_edges = label_boundary_edges(*he_mesh_, truth);
        truth_edges_.clear();
        for (const Edge& e : gt_edges) truth_edges_.insert(truth_edges_.end(), {e.v0, e.v1});
        std::vector<Edge> pred;
        for (std::size_t i = 0; i + 1 < margin_edges_.size(); i += 2) pred.push_back({margin_edges_[i], margin_edges_[i + 1]});
        const BoundaryMetrics b = compare_boundaries(edge_midpoints(*he_mesh_, pred), edge_midpoints(*he_mesh_, gt_edges));
        const double iou = region_iou(predicted_tooth_, truth, build_dec(*he_mesh_).star0);
        std::ostringstream o;
        o.precision(6);
        o << "{\"assd\":" << b.assd << ",\"hd95\":" << b.hd95 << ",\"f1_025\":" << b.f1_025 << ",\"f1_050\":" << b.f1_050
          << ",\"iou\":" << iou << "}";
        return o.str();
    }
    // Undercut map (M10b, D90) for the teeth (the detected tooth region if the margin was detected, else the whole
    // mesh), along the occlusal axis or, with `best`, along the best path of insertion within 25 degrees of it.
    // Per face: 0 outside the region, 1 reachable, 2 undercut. Returns JSON (axis, tilt, fraction).
    std::string computeUndercut(bool best) {
        if (!he_mesh_) return "{\"error\":\"undercut analysis needs a manifold analysis view\"}";
        if (!bvh_) bvh_ = std::make_unique<Bvh>(mesh_.positions, mesh_.triangles);
        std::vector<std::uint8_t> region;
        if (!predicted_tooth_.empty()) {
            region.assign(mesh_.triangles.size(), 0);
            for (std::size_t f = 0; f < mesh_.triangles.size(); ++f) {
                const auto& t = mesh_.triangles[f];
                region[f] = predicted_tooth_[t[0]] && predicted_tooth_[t[1]] && predicted_tooth_[t[2]];
            }
        }
        const Vec3 axis = occlusal_axis(*he_mesh_, largest_component_mask(*he_mesh_));
        Vec3 used = axis;
        UndercutResult r;
        std::size_t evaluations = 1;
        if (best) {
            InsertionAxis b = best_insertion_axis(mesh_, *bvh_, region, axis, 25.0);
            used = b.axis, r = std::move(b.result), evaluations = b.evaluations;
        } else {
            r = undercut_map(mesh_, *bvh_, axis, region);
        }
        undercut_faces_.assign(mesh_.triangles.size(), 0);
        for (std::size_t f = 0; f < undercut_faces_.size(); ++f)
            if (region.empty() || region[f]) undercut_faces_[f] = r.undercut[f] ? 2 : 1;
        const double tilt = std::acos(std::min(1.0, std::max(-1.0, used.x * axis.x + used.y * axis.y + used.z * axis.z))) * 180.0 / std::acos(-1.0);
        std::ostringstream o;
        o.precision(6);
        o << "{\"axis\":[" << used.x << "," << used.y << "," << used.z << "],\"tilt\":" << tilt << ",\"fraction\":" << r.fraction()
          << ",\"evaluations\":" << evaluations << ",\"teethOnly\":" << (region.empty() ? "false" : "true") << "}";
        return o.str();
    }
    val undercutFaces() const { return view(undercut_faces_); }  // uint8 per face: 0 outside, 1 reachable, 2 undercut

    val truthMarginEdges() const { return view(truth_edges_); }  // uint32 vertex pairs of all handle loops

    // Switch curvature and geodesics between the cotan and the intrinsic Delaunay Laplacian.
    void setIntrinsicDelaunay(bool on) {
        intrinsic_delaunay_ = on;
        update_curvature();
        geodesics_.reset();  // rebuilt with the new operators on the next query
        distance_.clear();
    }

    // Heat-method distance from one vertex (D51). The solver (and its two factorizations) is
    // built on the first query and reused for every later one on the same mesh.
    std::string geodesic(std::uint32_t source) {
        if (!he_mesh_) return "geodesic distance needs a consistently oriented manifold (" + he_error_ + ")";
        if (!geodesics_) {
            geodesics_ = std::make_unique<HeatGeodesics>(*he_mesh_, 1.0, intrinsic_delaunay_);
            if (!geodesics_->ok()) return geodesics_->error();
        }
        const std::uint32_t src[1] = {source};
        const GeodesicResult r = geodesics_->distance(src);
        if (!r.ok()) return r.error;
        distance_.assign(r.distance.begin(), r.distance.end());  // double -> float, NaN preserved
        return {};
    }
    val distance() const { return view(distance_); }  // float32 per vertex, NaN = unreachable

    // Summary as JSON (built by hand: no dependency for one small object; finite numbers only).
    std::string stats() const {
        std::ostringstream o;
        o.precision(17);
        auto opt_range = [&](const std::vector<float>& v) {
            float lo = std::numeric_limits<float>::infinity(), hi = -lo;
            for (float x : v) {
                if (std::isfinite(x)) lo = std::min(lo, x), hi = std::max(hi, x);
            }
            if (lo > hi) o << "null";
            else o << '[' << lo << ',' << hi << ']';
        };
        o << "{\"vertices\":" << mesh_.positions.size() << ",\"faces\":" << mesh_.triangles.size()
          << ",\"edges\":" << topo_.edges.size() << ",\"edgeKinds\":{\"boundary\":"
          << topo_.count(EdgeKind::Boundary) << ",\"manifold\":" << topo_.count(EdgeKind::Manifold)
          << ",\"misoriented\":" << topo_.count(EdgeKind::Misoriented)
          << ",\"nonmanifold\":" << topo_.count(EdgeKind::NonManifold) << "}"
          << ",\"invalidFaces\":" << topo_.invalid_faces.size()
          << ",\"duplicateFaces\":" << topo_.duplicate_faces.size()
          << ",\"isolatedVertices\":" << topo_.isolated_vertices.size()
          << ",\"nonmanifoldVertices\":" << topo_.nonmanifold_vertices.size() << ",\"components\":[";
        for (std::size_t i = 0; i < topo_.components.size(); ++i) {
            const auto& c = topo_.components[i];
            o << (i ? "," : "") << "{\"V\":" << c.num_vertices << ",\"E\":" << c.num_edges
              << ",\"F\":" << c.num_faces << ",\"chi\":" << c.euler_characteristic
              << ",\"b\":" << c.boundary_loops << ",\"manifold\":" << (c.manifold ? "true" : "false")
              << ",\"orientable\":" << (c.orientable ? "true" : "false")
              << ",\"consistent\":" << (c.consistently_oriented ? "true" : "false") << ",\"genus\":";
            if (c.genus) o << *c.genus; else o << "null";
            o << ",\"betti\":";
            if (c.betti) o << '[' << (*c.betti)[0] << ',' << (*c.betti)[1] << ',' << (*c.betti)[2] << ']';
            else o << "null";
            o << '}';
        }
        o << "],\"handleLoops\":{\"count\":" << handle_lengths_.size() << ",\"lengths\":[";
        for (std::size_t i = 0; i < handle_lengths_.size() && i < 64; ++i) o << (i ? "," : "") << handle_lengths_[i];
        o << "]}";
        o << ",\"excludedFaces\":" << excluded_faces_;
        o << ",\"intrinsicDelaunay\":" << (intrinsic_delaunay_ ? "true" : "false")
          << ",\"flips\":" << curvature_.intrinsic_flips << ",\"halfedge\":{\"ok\":" << (he_error_.empty() ? "true" : "false") << ",\"error\":\""
          << json_escape(he_error_) << "\"},\"curvature\":";
        if (!he_error_.empty()) {
            o << "null";
        } else {
            double total_defect = 0.0;
            for (double d : curvature_.angle_defect) total_defect += d;
            o << "{\"totalAngleDefect\":" << total_defect << ",\"meanRange\":";
            opt_range(mean_);
            o << ",\"gaussianRange\":";
            opt_range(gaussian_);
            o << ",\"degenerateFaces\":" << curvature_.degenerate_faces.size() << '}';
        }
        o << '}';
        return o.str();
    }

private:
    void set_mesh(TriMesh m) {
        mesh_ = std::move(m);
        topo_ = analyze_topology(mesh_);
        // Analysis view (D63, D69): faces at reported defects are excluded so the half-edge-based
        // algorithms run on real scans; the topology report above still describes the full input.
        AnalysisMesh analysis = manifold_analysis_mesh(mesh_);
        excluded_faces_ = analysis.excluded_faces;
        BuildResult he;
        if (analysis.manifold) he.mesh = std::move(analysis.halfedge);
        else he.error = "not a manifold even after excluding faces at defects";
        he_error_ = he.error;
        geodesics_.reset();  // rebuilt lazily for each new mesh
        he_mesh_ = he.ok() ? std::make_unique<HalfEdgeMesh>(std::move(he.mesh)) : nullptr;
        distance_.clear();
        cusps_.clear();
        predicted_tooth_.clear();
        margin_edges_.clear();
        bvh_.reset();
        undercut_faces_.clear();
        truth_edges_.clear();
        update_curvature();
        // Handle loops (M8b): 2g generator cycles per component, shortest first (lengths are upper
        // bounds on handle size, D64). Stored as line segments for the viewer.
        handle_edges_.clear();
        handle_lengths_.clear();
        if (he_mesh_) {
            for (const auto& loop : handle_loops(*he_mesh_)) {
                handle_lengths_.push_back(loop.length);
                for (std::size_t i = 0; i < loop.vertices.size(); ++i) {
                    handle_edges_.push_back(loop.vertices[i]);
                    handle_edges_.push_back(loop.vertices[(i + 1) % loop.vertices.size()]);
                }
            }
        }

        // GPU-facing buffers: float32 positions (D12: double in the core, float at the GPU).
        positions_.clear();
        for (const auto& p : mesh_.positions) {
            positions_.insert(positions_.end(), {static_cast<float>(p.x), static_cast<float>(p.y),
                                                 static_cast<float>(p.z)});
        }
        indices_.clear();
        for (const auto& t : mesh_.triangles) indices_.insert(indices_.end(), t.begin(), t.end());
        face_component_ = topo_.face_component;
        boundary_edges_.clear();
        nonmanifold_edges_.clear();
        misoriented_edges_.clear();
        for (std::size_t e = 0; e < topo_.edges.size(); ++e) {
            std::vector<std::uint32_t>* dst = nullptr;
            switch (topo_.edge_kind[e]) {
                case EdgeKind::Boundary: dst = &boundary_edges_; break;
                case EdgeKind::NonManifold: dst = &nonmanifold_edges_; break;
                case EdgeKind::Misoriented: dst = &misoriented_edges_; break;
                case EdgeKind::Manifold: break;
            }
            if (dst) dst->insert(dst->end(), {topo_.edges[e].v0, topo_.edges[e].v1});
        }
    }

    void update_curvature() {
        curvature_ = he_mesh_ ? compute_curvature(*he_mesh_, intrinsic_delaunay_) : CurvatureField{};
        auto to_f32 = [](const std::vector<double>& v, std::size_t n) {
            std::vector<float> out(n, std::numeric_limits<float>::quiet_NaN());
            for (std::size_t i = 0; i < v.size() && i < n; ++i) out[i] = static_cast<float>(v[i]);
            return out;
        };
        mean_ = to_f32(curvature_.mean, mesh_.positions.size());
        gaussian_ = to_f32(curvature_.gaussian, mesh_.positions.size());
        kmin_ = to_f32(curvature_.k2, mesh_.positions.size());
    }

    bool intrinsic_delaunay_ = false;
    std::size_t excluded_faces_ = 0;
    TriMesh mesh_;
    TopologyReport topo_;
    std::string he_error_;
    CurvatureField curvature_;
    std::unique_ptr<HalfEdgeMesh> he_mesh_;
    std::unique_ptr<HeatGeodesics> geodesics_;
    std::vector<float> distance_;
    std::vector<float> positions_, mean_, gaussian_, kmin_;
    std::vector<std::uint32_t> indices_, face_component_, boundary_edges_, nonmanifold_edges_,
        misoriented_edges_, handle_edges_;
    std::vector<double> handle_lengths_;
    std::vector<std::uint32_t> cusps_, margin_edges_, truth_edges_;
    std::vector<std::uint8_t> predicted_tooth_;
    std::unique_ptr<Bvh> bvh_;
    std::vector<std::uint8_t> undercut_faces_;
};

}  // namespace

EMSCRIPTEN_BINDINGS(dmw) {
    emscripten::class_<Session>("Session")
        .constructor<>()
        .function("load", &Session::load)
        .function("generate", &Session::generate)
        .function("positions", &Session::positions)
        .function("indices", &Session::indices)
        .function("mean", &Session::mean)
        .function("gaussian", &Session::gaussian)
        .function("kmin", &Session::kmin)
        .function("faceComponent", &Session::faceComponent)
        .function("boundaryEdges", &Session::boundaryEdges)
        .function("nonmanifoldEdges", &Session::nonmanifoldEdges)
        .function("misorientedEdges", &Session::misorientedEdges)
        .function("nonmanifoldVertices", &Session::nonmanifoldVertices)
        .function("handleEdges", &Session::handleEdges)
        .function("detectCusps", &Session::detectCusps)
        .function("cusps", &Session::cusps)
        .function("detectMargin", &Session::detectMargin)
        .function("marginEdges", &Session::marginEdges)
        .function("compareMargin", &Session::compareMargin)
        .function("truthMarginEdges", &Session::truthMarginEdges)
        .function("computeUndercut", &Session::computeUndercut)
        .function("undercutFaces", &Session::undercutFaces)
        .function("setIntrinsicDelaunay", &Session::setIntrinsicDelaunay)
        .function("geodesic", &Session::geodesic)
        .function("distance", &Session::distance)
        .function("stats", &Session::stats);
}
