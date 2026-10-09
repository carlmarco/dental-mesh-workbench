// Wall thickness on real geometry (M10c, D92). Scans are open surfaces (no thickness), so each tooth crown
// (faces whose three vertices carry one FDI label) is turned into a solid the naive way: the crown surface, a copy
// moved inward by a nominal thickness along the vertex normals (orientation flipped), and a side wall along the
// boundary. The thickness check then shows where this naive offset fails (concave regions fold the inner surface),
// the motivation for offsets from a signed distance field (M10d-e).
//   thickness_eval <scan-dir> [stride] [nominal_mm] [sdf]
// sdf (D95): the shell is built from the crown's generalized signed distance (signed heat method, M10d) instead:
// outer wall = the phi = 0 iso-surface, inner wall = phi = -nominal, flipped; both closed, no stitching.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "core/bvh.h"
#include "core/isosurface.h"
#include "core/sdf.h"
#include "core/io.h"
#include "core/margin.h"
#include "core/thickness.h"
#include "dataset.h"

using namespace dmw;
using Clock = std::chrono::steady_clock;

namespace {
double median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}
const char* type_name(int unit) {
    switch (unit) {
        case 1: case 2: return "incisors";
        case 3: return "canines";
        case 4: case 5: return "premolars";
        default: return "molars";
    }
}
// Shell from the faces of `m` in `faces`: outer = the faces, inner = vertices moved by -t along the normal, flipped.
// Returns the shell and the number of outer vertices (they come first).
TriMesh naive_shell(const TriMesh& m, const std::vector<std::uint32_t>& faces, double t, std::size_t& outer_count,
                    std::vector<std::uint8_t>& near_rim) {
    std::map<std::uint32_t, std::uint32_t> local;
    TriMesh s;
    for (std::uint32_t f : faces)
        for (auto v : m.triangles[f])
            if (local.emplace(v, static_cast<std::uint32_t>(s.positions.size())).second) s.positions.push_back(m.positions[v]);
    const auto n = static_cast<std::uint32_t>(s.positions.size());
    outer_count = n;
    std::vector<Vec3> normal(n, Vec3{});
    std::map<std::pair<std::uint32_t, std::uint32_t>, int> edge_use;  // directed edges of the outer faces
    for (std::uint32_t f : faces) {
        const std::uint32_t a = local[m.triangles[f][0]], b = local[m.triangles[f][1]], c = local[m.triangles[f][2]];
        s.triangles.push_back({a, b, c});
        const Vec3 &pa = s.positions[a], &pb = s.positions[b], &pc = s.positions[c];
        const Vec3 u{pb.x - pa.x, pb.y - pa.y, pb.z - pa.z}, w{pc.x - pa.x, pc.y - pa.y, pc.z - pa.z};
        const Vec3 nn{u.y * w.z - u.z * w.y, u.z * w.x - u.x * w.z, u.x * w.y - u.y * w.x};
        for (std::uint32_t v : {a, b, c}) normal[v].x += nn.x, normal[v].y += nn.y, normal[v].z += nn.z;
        for (auto [p, q] : {std::pair{a, b}, std::pair{b, c}, std::pair{c, a}}) ++edge_use[{p, q}];
    }
    for (std::uint32_t v = 0; v < n; ++v) {
        const Vec3& nv = normal[v];
        const double l = std::sqrt(nv.x * nv.x + nv.y * nv.y + nv.z * nv.z);
        const Vec3& p = s.positions[v];
        s.positions.push_back(l > 0 ? Vec3{p.x - t * nv.x / l, p.y - t * nv.y / l, p.z - t * nv.z / l} : p);
    }
    const std::size_t outer_faces = s.triangles.size();
    for (std::size_t i = 0; i < outer_faces; ++i) {
        const auto tri = s.triangles[i];
        s.triangles.push_back({tri[0] + n, tri[2] + n, tri[1] + n});  // inner, flipped
    }
    // Outer vertices within 1 mm (straight line) of the crown boundary: their cone rays can reach the side wall.
    std::vector<Vec3> rim;
    for (const auto& [e, count] : edge_use)
        if (!edge_use.count({e.second, e.first})) rim.push_back(s.positions[e.first]);
    near_rim.assign(n, 0);
    for (std::uint32_t v = 0; v < n; ++v)
        for (const Vec3& q : rim) {
            const Vec3& p = s.positions[v];
            if ((p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y) + (p.z - q.z) * (p.z - q.z) < 1.0) {
                near_rim[v] = 1;
                break;
            }
        }
    for (const auto& [e, count] : edge_use) {
        if (edge_use.count({e.second, e.first})) continue;  // interior edge
        const std::uint32_t a = e.first, b = e.second;      // boundary edge a -> b of an outer face
        s.triangles.push_back({b, a, a + n});
        s.triangles.push_back({b, a + n, b + n});
    }
    return s;
}
}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: thickness_eval <scan-dir> [stride] [nominal_mm]\n");
        return 2;
    }
    const std::size_t stride = argc > 2 ? std::stoul(argv[2]) : 30;
    const double nominal = argc > 3 ? std::stod(argv[3]) : 0.8;
    const bool sdf = argc > 4 && std::string(argv[4]) == "sdf";
    const auto objs = dataset::index_files({argv[1]}, ".obj");
    const auto labels = dataset::index_files({argv[1]}, ".json");
    struct PerType {
        std::vector<double> thin, thin_normal, min_t, ms;
    };
    std::map<std::string, PerType> types;
    std::vector<double> rays_per_s;
    std::size_t index = 0, scans = 0, teeth = 0;
    for (const auto& [stem, obj] : objs) {
        if (index++ % stride != 0) continue;
        const auto lab = labels.find(stem);
        if (lab == labels.end()) continue;
        const JsonResult j = parse_json(dataset::read_text(lab->second));
        const LoadResult r = parse_obj(dataset::read_text(obj));
        const JsonValue* arr = j.ok() ? j.value.find("labels") : nullptr;
        if (!r.ok() || !arr || arr->array.size() != r.mesh.positions.size()) continue;
        std::map<int, std::vector<std::uint32_t>> crowns;
        for (std::uint32_t f = 0; f < r.mesh.triangles.size(); ++f) {
            const auto& t = r.mesh.triangles[f];
            const int l0 = int(arr->array[t[0]].number);
            if (l0 != 0 && l0 == int(arr->array[t[1]].number) && l0 == int(arr->array[t[2]].number)) crowns[l0].push_back(f);
        }
        for (const auto& [fdi, faces] : crowns) {
            std::size_t outer = 0;
            std::vector<std::uint8_t> near_rim;
            TriMesh shell = naive_shell(r.mesh, faces, nominal, outer, near_rim);
            if (sdf) {
                // Same crown, SDF shell. Measure on outer-wall vertices that lie on the original crown surface (within
                // 0.1 mm of its vertices) and more than 1 mm from its rim, matching the naive rows.
                TriMesh crown;
                crown.positions.assign(shell.positions.begin(), shell.positions.begin() + static_cast<std::ptrdiff_t>(outer));
                for (std::size_t f = 0; f < faces.size(); ++f) crown.triangles.push_back(shell.triangles[f]);
                std::vector<Vec3> rim_pts, crown_pts;
                for (std::size_t v = 0; v < outer; ++v) (near_rim[v] ? rim_pts : crown_pts).push_back(crown.positions[v]);
                SignedHeatParams prm;
                prm.h = 0.12, prm.padding = 2.5;
                const Grid3 phi = signed_heat_distance(crown, prm);
                TriMesh outer_wall = extract_isosurface(phi, 0.0), inner_wall = extract_isosurface(phi, -nominal);
                for (auto& t : inner_wall.triangles) std::swap(t[1], t[2]);  // normals out of the solid (into the hollow)
                shell = outer_wall;
                const auto offset = static_cast<std::uint32_t>(shell.positions.size());
                shell.positions.insert(shell.positions.end(), inner_wall.positions.begin(), inner_wall.positions.end());
                for (auto t : inner_wall.triangles) shell.triangles.push_back({t[0] + offset, t[1] + offset, t[2] + offset});
                outer = outer_wall.positions.size();
                const auto to_crown = nearest_point_distances(outer_wall.positions, crown_pts, 0.5);
                const auto to_rim = nearest_point_distances(outer_wall.positions, rim_pts, 0.5);
                near_rim.assign(outer, 0);
                for (std::size_t v = 0; v < outer; ++v) near_rim[v] = to_crown[v] > 0.1 || to_rim[v] < 0.1 ? 1 : 0;  // excluded
            }
            const auto t0 = Clock::now();
            const Bvh bvh(shell.positions, shell.triangles);
            ThicknessParams prm;
            const ThicknessField th = wall_thickness(shell, bvh, prm);
            const double ms = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
            rays_per_s.push_back(double(shell.positions.size()) * double(prm.rays + 1) / (ms / 1000.0));
            std::size_t measured = 0, thin = 0, thin_n = 0;
            double mn = 1e300;
            for (std::size_t v = 0; v < outer; ++v) {
                const double x = th.cone_min[v];
                if (near_rim[v] || !std::isfinite(x)) continue;  // away from the rim only
                ++measured, thin += x < 0.95 * nominal, mn = std::min(mn, x);
                thin_n += std::isfinite(th.along_normal[v]) && th.along_normal[v] < 0.95 * nominal;
            }
            if (measured == 0) continue;
            PerType& pt = types[type_name(fdi % 10)];
            pt.thin.push_back(100.0 * double(thin) / double(measured)), pt.min_t.push_back(mn), pt.ms.push_back(ms);
            pt.thin_normal.push_back(100.0 * double(thin_n) / double(measured));
            ++teeth;
        }
        ++scans;
        std::fprintf(stderr, "\r%zu scans, %zu teeth", scans, teeth);
    }
    std::fprintf(stderr, "\n");
    std::printf("%s %.2f mm shells from %zu crowns on %zu scans; thickness throughput median %.2f M rays/s (all threads)\n\n",
                sdf ? "signed-distance" : "naive", nominal,
                teeth, scans, median(rays_per_s) / 1e6);
    std::printf("outer vertices more than 1 mm from the crown rim only\n");
    std::printf("| tooth type | crowns | thinner than 95%% of nominal, cone minimum (median %%) | same, normal ray only (median %%) | minimum wall (median mm) | ms per crown (BVH + thickness) |\n");
    std::printf("|---|---:|---:|---:|---:|---:|\n");
    for (const auto& [name, pt] : types)
        std::printf("| %s | %zu | %.1f | %.1f | %.3f | %.0f |\n", name.c_str(), pt.thin.size(), median(pt.thin), median(pt.thin_normal),
                    median(pt.min_t), median(pt.ms));
    return 0;
}
