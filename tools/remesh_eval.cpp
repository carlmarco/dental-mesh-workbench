// Isotropic remeshing on real scans (M10f, D99): quality and fidelity before vs after, at the scan's own mean edge
// length. Fidelity: distance from the ORIGINAL vertices to the remeshed surface (BVH closest point), so lost detail
// shows up; the remeshed vertices themselves are projected onto the original and sit on it by construction.
//   remesh_eval <scan-dir> [stride]
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

#include "core/bvh.h"
#include "core/halfedge.h"
#include "core/io.h"
#include "core/remesh.h"
#include "core/topology.h"
#include "dataset.h"

using namespace dmw;

namespace {
double median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}
double angle_at(const Vec3& a, const Vec3& b, const Vec3& c) {
    const Vec3 u{b.x - a.x, b.y - a.y, b.z - a.z}, v{c.x - a.x, c.y - a.y, c.z - a.z};
    const double nu = std::sqrt(u.x * u.x + u.y * u.y + u.z * u.z), nv = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    return std::acos(std::clamp((u.x * v.x + u.y * v.y + u.z * v.z) / std::max(nu * nv, 1e-300), -1.0, 1.0)) * 180.0 / std::acos(-1.0);
}
struct Quality {
    double min_angle_q10, edge_cv;
};
Quality quality(const TriMesh& m) {
    std::vector<double> ang, len;
    for (const auto& t : m.triangles) {
        const Vec3 &a = m.positions[t[0]], &b = m.positions[t[1]], &c = m.positions[t[2]];
        ang.push_back(std::min({angle_at(a, b, c), angle_at(b, c, a), angle_at(c, a, b)}));
        for (int k = 0; k < 3; ++k) {
            const Vec3 &p = m.positions[t[static_cast<std::size_t>(k)]], &q = m.positions[t[static_cast<std::size_t>((k + 1) % 3)]];
            len.push_back(std::sqrt((p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y) + (p.z - q.z) * (p.z - q.z)));
        }
    }
    std::sort(ang.begin(), ang.end());
    double mean = 0.0, sq = 0.0;
    for (double l : len) mean += l;
    mean /= double(len.size());
    for (double l : len) sq += (l - mean) * (l - mean);
    return {ang[ang.size() / 10], std::sqrt(sq / double(len.size())) / mean};
}
}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: remesh_eval <scan-dir> [stride]\n");
        return 2;
    }
    const std::size_t stride = argc > 2 ? std::stoul(argv[2]) : 60;
    const auto objs = dataset::index_files({argv[1]}, ".obj");
    std::size_t index = 0, scans = 0, valid = 0, same_topology = 0;
    std::vector<double> q_before, q_after, cv_before, cv_after, mean_dev, max_dev, seconds, face_ratio;
    for (const auto& [stem, obj] : objs) {
        if (index++ % stride != 0) continue;
        const LoadResult r = parse_obj(dataset::read_text(obj));
        if (!r.ok() || !build_halfedge(r.mesh).ok()) continue;
        ++scans;
        const auto t0 = std::chrono::steady_clock::now();
        const TriMesh out = remesh_isotropic(r.mesh, RemeshParams{});
        seconds.push_back(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        valid += build_halfedge(out).ok();
        const auto a = analyze_topology(r.mesh), b = analyze_topology(out);
        same_topology += a.components.size() == b.components.size() && !a.components.empty() &&
                         a.components[0].genus == b.components[0].genus && a.components[0].boundary_loops == b.components[0].boundary_loops;
        const Quality qa = quality(r.mesh), qb = quality(out);
        q_before.push_back(qa.min_angle_q10), q_after.push_back(qb.min_angle_q10), cv_before.push_back(qa.edge_cv), cv_after.push_back(qb.edge_cv);
        face_ratio.push_back(double(out.triangles.size()) / double(r.mesh.triangles.size()));
        const Bvh bvh(out.positions, out.triangles);
        double s = 0.0, mx = 0.0;
        for (const Vec3& p : r.mesh.positions) {
            const double d = bvh.closest(p).distance;
            s += d, mx = std::max(mx, d);
        }
        mean_dev.push_back(s / double(r.mesh.positions.size())), max_dev.push_back(mx);
        std::fprintf(stderr, "\r%zu scans", scans);
    }
    std::fprintf(stderr, "\n");
    std::printf("scans: %zu (manifold); remeshed result a valid oriented manifold on %zu, same topology (components, genus,"
                " boundary loops) on %zu\n", scans, valid, same_topology);
    std::printf("10th-percentile minimum angle (median over scans): %.1f -> %.1f deg; edge-length CV: %.2f -> %.2f; faces x%.2f\n",
                median(q_before), median(q_after), median(cv_before), median(cv_after), median(face_ratio));
    std::printf("original vertices to the remeshed surface: mean %.4f mm, max %.3f mm (medians over scans); time median %.1f s\n",
                median(mean_dev), median(max_dev), median(seconds));
    return 0;
}
