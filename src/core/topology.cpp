#include "core/topology.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "detail/disjoint_sets.h"
#include "detail/hash.h"

namespace dmw {

std::size_t TopologyReport::count(EdgeKind kind) const {
    return static_cast<std::size_t>(std::count(edge_kind.begin(), edge_kind.end(), kind));
}

namespace {

using detail::DisjointSets;
using detail::pack_pair;

// Per undirected edge: how many valid faces use it, and for the first two, which face and
// whether that face traverses it "forward" (v0 -> v1, i.e. low index to high index).
struct EdgeIncidence {
    std::uint32_t count = 0;
    std::array<std::uint32_t, 2> face{kInvalid, kInvalid};
    std::array<bool, 2> forward{false, false};
};

// Compressed adjacency ("CSR"): the neighbors of item i are items[offset[i] .. offset[i+1]).
// One flat array instead of a vector per item: two allocations total, cache-friendly scans.
template <typename T>
struct Csr {
    std::vector<std::uint32_t> offset;
    std::vector<T> items;
};

}  // namespace

TopologyReport analyze_topology(const TriMesh& mesh) {
    TopologyReport r;
    const std::size_t nv = mesh.positions.size();
    const std::size_t nf = mesh.triangles.size();
    r.vertex_component.assign(nv, kInvalid);
    r.face_component.assign(nf, kInvalid);

    // --- 1. Valid faces (D31): indices in range, pairwise distinct, first copy of their
    // vertex set. Everything below looks only at valid faces.
    std::vector<std::uint32_t> valid;
    valid.reserve(nf);
    std::unordered_set<std::array<std::uint32_t, 3>, detail::Hash3<std::uint32_t>> seen;
    seen.reserve(nf);
    for (std::size_t f = 0; f < nf; ++f) {
        auto t = mesh.triangles[f];
        const auto fid = static_cast<std::uint32_t>(f);
        if (t[0] >= nv || t[1] >= nv || t[2] >= nv || t[0] == t[1] || t[1] == t[2] || t[0] == t[2]) {
            r.invalid_faces.push_back(fid);
            continue;
        }
        std::sort(t.begin(), t.end());  // vertex set, ignoring order and orientation
        if (!seen.insert(t).second) {
            r.duplicate_faces.push_back(fid);
            continue;
        }
        valid.push_back(fid);
    }

    // --- 2. Edge incidence and classification.
    std::unordered_map<std::uint64_t, std::uint32_t> edge_id;
    edge_id.reserve(3 * valid.size());
    std::vector<EdgeIncidence> inc;
    for (std::uint32_t f : valid) {
        const auto& t = mesh.triangles[f];
        for (std::size_t k = 0; k < 3; ++k) {
            const std::uint32_t u = t[k], v = t[(k + 1) % 3];
            const std::uint32_t lo = std::min(u, v), hi = std::max(u, v);
            const auto [it, inserted] =
                edge_id.try_emplace(pack_pair(lo, hi), static_cast<std::uint32_t>(r.edges.size()));
            if (inserted) {
                r.edges.push_back({lo, hi});
                inc.emplace_back();
            }
            EdgeIncidence& e = inc[it->second];
            if (e.count < 2) {
                e.face[e.count] = f;
                e.forward[e.count] = (u < v);
            }
            ++e.count;
        }
    }
    r.edge_kind.reserve(r.edges.size());
    for (const auto& e : inc) {
        if (e.count == 1) {
            r.edge_kind.push_back(EdgeKind::Boundary);
        } else if (e.count == 2) {
            // Consistent winding: the two faces cross the edge in opposite directions.
            r.edge_kind.push_back(e.forward[0] != e.forward[1] ? EdgeKind::Manifold
                                                               : EdgeKind::Misoriented);
        } else {
            r.edge_kind.push_back(EdgeKind::NonManifold);
        }
    }

    // --- 3. Connected components by vertex connectivity (D28): union the corners of each
    // valid face. Component ids are assigned in order of lowest vertex index.
    DisjointSets ds(nv);
    std::vector<bool> used(nv, false);
    for (std::uint32_t f : valid) {
        const auto& t = mesh.triangles[f];
        ds.unite(t[0], t[1]);
        ds.unite(t[1], t[2]);
        used[t[0]] = used[t[1]] = used[t[2]] = true;
    }
    std::vector<std::uint32_t> component_of_root(nv, kInvalid);
    for (std::uint32_t v = 0; v < nv; ++v) {
        if (!used[v]) {
            r.isolated_vertices.push_back(v);
            continue;
        }
        std::uint32_t& c = component_of_root[ds.find(v)];
        if (c == kInvalid) {
            c = static_cast<std::uint32_t>(r.components.size());
            r.components.emplace_back();
        }
        r.vertex_component[v] = c;
        ++r.components[c].num_vertices;
    }
    for (std::uint32_t f : valid) {
        const std::uint32_t c = r.vertex_component[mesh.triangles[f][0]];
        r.face_component[f] = c;
        ++r.components[c].num_faces;
    }
    for (std::size_t e = 0; e < r.edges.size(); ++e) {
        ComponentTopology& c = r.components[r.vertex_component[r.edges[e].v0]];
        ++c.num_edges;
        if (r.edge_kind[e] == EdgeKind::NonManifold) c.manifold = false;
        if (r.edge_kind[e] == EdgeKind::Misoriented) c.consistently_oriented = false;
    }

    // --- 4. Non-manifold (bowtie) vertices. Around v, two incident faces belong to the same
    // fan if they share an edge (v, w). Collect (w, face) for both edges at v of every incident
    // face, sort by w, and union faces with equal w. More than one resulting set = bowtie.
    Csr<std::uint32_t> vf;  // vertex -> incident valid faces
    vf.offset.assign(nv + 1, 0);
    for (std::uint32_t f : valid) {
        for (std::uint32_t v : mesh.triangles[f]) ++vf.offset[v + 1];
    }
    for (std::size_t v = 0; v < nv; ++v) vf.offset[v + 1] += vf.offset[v];  // prefix sums
    vf.items.resize(vf.offset[nv]);
    {
        std::vector<std::uint32_t> fill(vf.offset.begin(), vf.offset.end() - 1);
        for (std::uint32_t f : valid) {
            for (std::uint32_t v : mesh.triangles[f]) vf.items[fill[v]++] = f;
        }
    }
    DisjointSets local;
    std::vector<std::pair<std::uint32_t, std::uint32_t>> spokes;  // (w, local face index)
    for (std::uint32_t v = 0; v < nv; ++v) {
        const std::uint32_t begin = vf.offset[v], end = vf.offset[v + 1];
        const std::uint32_t deg = end - begin;
        if (deg < 2) continue;  // 0 or 1 face: trivially a single fan
        spokes.clear();
        for (std::uint32_t i = 0; i < deg; ++i) {
            for (std::uint32_t w : mesh.triangles[vf.items[begin + i]]) {
                if (w != v) spokes.emplace_back(w, i);
            }
        }
        std::sort(spokes.begin(), spokes.end());
        local.reset(deg);
        std::uint32_t fans = deg;
        for (std::size_t s = 1; s < spokes.size(); ++s) {
            if (spokes[s].first == spokes[s - 1].first && local.unite(spokes[s].second, spokes[s - 1].second)) {
                --fans;
            }
        }
        if (fans > 1) {
            r.nonmanifold_vertices.push_back(v);
            r.components[r.vertex_component[v]].manifold = false;
        }
    }

    // --- 5. Boundary loops: connected components of the graph of boundary edges. On a manifold
    // component every boundary vertex has exactly two boundary edges, so these are disjoint
    // cycles and counting sets counts loops.
    DisjointSets bds(nv);
    std::vector<bool> on_boundary(nv, false);
    for (std::size_t e = 0; e < r.edges.size(); ++e) {
        if (r.edge_kind[e] != EdgeKind::Boundary) continue;
        bds.unite(r.edges[e].v0, r.edges[e].v1);
        on_boundary[r.edges[e].v0] = on_boundary[r.edges[e].v1] = true;
    }
    for (std::uint32_t v = 0; v < nv; ++v) {
        if (on_boundary[v] && bds.find(v) == v) ++r.components[r.vertex_component[v]].boundary_loops;
    }

    // --- 6. Orientability (D30): try to choose flip[f] in {0,1} so that every two-face edge
    // becomes consistently wound. Across an edge whose faces already disagree (Misoriented),
    // exactly one of the pair must flip; across a Manifold edge, both or neither. Breadth-first
    // propagation from one face per component; a contradiction means no choice works.
    Csr<std::pair<std::uint32_t, bool>> adj;  // face -> (neighbor, must_differ)
    adj.offset.assign(nf + 1, 0);
    for (const auto& e : inc) {
        if (e.count != 2) continue;
        ++adj.offset[e.face[0] + 1];
        ++adj.offset[e.face[1] + 1];
    }
    for (std::size_t f = 0; f < nf; ++f) adj.offset[f + 1] += adj.offset[f];
    adj.items.resize(adj.offset[nf]);
    {
        std::vector<std::uint32_t> fill(adj.offset.begin(), adj.offset.end() - 1);
        for (const auto& e : inc) {
            if (e.count != 2) continue;
            const bool must_differ = (e.forward[0] == e.forward[1]);
            adj.items[fill[e.face[0]]++] = {e.face[1], must_differ};
            adj.items[fill[e.face[1]]++] = {e.face[0], must_differ};
        }
    }
    for (auto& c : r.components) c.orientable = c.manifold;
    std::vector<std::int8_t> flip(nf, -1);  // -1 = not yet visited
    std::vector<std::uint32_t> queue;
    for (std::uint32_t start : valid) {
        if (flip[start] != -1) continue;
        flip[start] = 0;
        queue.assign(1, start);
        for (std::size_t q = 0; q < queue.size(); ++q) {
            const std::uint32_t f = queue[q];
            for (std::uint32_t i = adj.offset[f]; i < adj.offset[f + 1]; ++i) {
                const auto [g, must_differ] = adj.items[i];
                const auto want = static_cast<std::int8_t>(flip[f] ^ (must_differ ? 1 : 0));
                if (flip[g] == -1) {
                    flip[g] = want;
                    queue.push_back(g);
                } else if (flip[g] != want) {
                    r.components[r.face_component[f]].orientable = false;
                }
            }
        }
    }

    // --- 7. Euler characteristic, genus, Betti numbers (D29). For a compact connected
    // orientable surface of genus g with b boundary loops, chi = 2 - 2g - b.
    for (auto& c : r.components) {
        c.euler_characteristic = std::int64_t{c.num_vertices} - std::int64_t{c.num_edges} +
                                 std::int64_t{c.num_faces};
        if (!c.manifold || !c.orientable) continue;
        const std::int64_t twice_genus = 2 - c.euler_characteristic - std::int64_t{c.boundary_loops};
        if (twice_genus < 0 || twice_genus % 2 != 0) continue;  // impossible for a valid surface
        const auto g = static_cast<std::uint32_t>(twice_genus / 2);
        const std::uint32_t b = c.boundary_loops;
        c.genus = g;
        c.betti = std::array<std::uint32_t, 3>{1, b > 0 ? 2 * g + b - 1 : 2 * g, b == 0 ? 1u : 0u};
    }
    return r;
}

}  // namespace dmw

namespace dmw {

TriMesh without_excluded_faces(const TriMesh& mesh, const TopologyReport& report) {
    std::vector<bool> drop(mesh.triangles.size(), false);
    for (std::uint32_t f : report.invalid_faces) drop[f] = true;
    for (std::uint32_t f : report.duplicate_faces) drop[f] = true;
    TriMesh out;
    out.positions = mesh.positions;
    out.triangles.reserve(mesh.triangles.size());
    for (std::size_t f = 0; f < mesh.triangles.size(); ++f) {
        if (!drop[f]) out.triangles.push_back(mesh.triangles[f]);
    }
    return out;
}

}  // namespace dmw

namespace dmw {

AnalysisMesh manifold_analysis_mesh(const TriMesh& input, int max_passes) {
    AnalysisMesh out;
    out.mesh = input;
    for (out.passes = 0; out.passes <= max_passes; ++out.passes) {
        BuildResult built = build_halfedge(out.mesh);
        if (built.ok()) {
            out.halfedge = std::move(built.mesh);
            out.manifold = true;
            break;
        }
        if (out.passes == max_passes) break;
        const TopologyReport r = analyze_topology(out.mesh);
        std::vector<bool> drop(out.mesh.triangles.size(), false);
        for (std::uint32_t f : r.invalid_faces) drop[f] = true;
        for (std::uint32_t f : r.duplicate_faces) drop[f] = true;
        std::vector<bool> bad_vertex(out.mesh.positions.size(), false);
        for (std::uint32_t v : r.nonmanifold_vertices) bad_vertex[v] = true;
        // Edges to cut around: non-manifold (3+ faces) and misoriented (winding flips there).
        std::unordered_set<std::uint64_t> bad_edge;
        for (std::size_t e = 0; e < r.edges.size(); ++e) {
            if (r.edge_kind[e] == EdgeKind::NonManifold || r.edge_kind[e] == EdgeKind::Misoriented) {
                bad_edge.insert(detail::pack_pair(r.edges[e].v0, r.edges[e].v1));
            }
        }
        for (std::size_t f = 0; f < out.mesh.triangles.size(); ++f) {
            if (drop[f]) continue;
            const auto& t = out.mesh.triangles[f];
            for (std::size_t k = 0; k < 3 && !drop[f]; ++k) {
                const std::uint32_t a = t[k], b = t[(k + 1) % 3];
                if (bad_vertex[a] || bad_edge.count(detail::pack_pair(std::min(a, b), std::max(a, b))) != 0) drop[f] = true;
            }
        }
        std::vector<std::array<std::uint32_t, 3>> kept;
        kept.reserve(out.mesh.triangles.size());
        for (std::size_t f = 0; f < out.mesh.triangles.size(); ++f) {
            if (!drop[f]) kept.push_back(out.mesh.triangles[f]);
        }
        out.mesh.triangles = std::move(kept);
    }
    out.excluded_faces = input.triangles.size() - out.mesh.triangles.size();
    return out;
}

}  // namespace dmw
