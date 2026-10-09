#include <algorithm>
#include <cmath>
#include <map>
#include <random>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "core/bvh.h"
#include "core/generate.h"
#include "core/halfedge.h"
#include "core/remesh.h"
#include "core/topology.h"

using namespace dmw;

namespace {

double angle_at(const Vec3& a, const Vec3& b, const Vec3& c) {  // angle at a
    const Vec3 u{b.x - a.x, b.y - a.y, b.z - a.z}, v{c.x - a.x, c.y - a.y, c.z - a.z};
    const double nu = std::sqrt(u.x * u.x + u.y * u.y + u.z * u.z), nv = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    return std::acos(std::clamp((u.x * v.x + u.y * v.y + u.z * v.z) / (nu * nv), -1.0, 1.0)) * 180.0 / std::acos(-1.0);
}
// 10th percentile of the smallest angle per triangle (degrees).
double min_angle_q10(const TriMesh& m) {
    std::vector<double> v;
    for (const auto& t : m.triangles) {
        const Vec3 &a = m.positions[t[0]], &b = m.positions[t[1]], &c = m.positions[t[2]];
        v.push_back(std::min({angle_at(a, b, c), angle_at(b, c, a), angle_at(c, a, b)}));
    }
    std::sort(v.begin(), v.end());
    return v[v.size() / 10];
}
// Share of edges within [4/5 L, 4/3 L].
double in_band(const TriMesh& m, double L) {
    std::size_t in = 0, all = 0;
    for (const auto& t : m.triangles)
        for (int k = 0; k < 3; ++k) {
            const Vec3 &a = m.positions[t[static_cast<std::size_t>(k)]], &b = m.positions[t[static_cast<std::size_t>((k + 1) % 3)]];
            const double l = std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
            in += l >= 0.8 * L - 1e-12 && l <= 4.0 / 3.0 * L + 1e-12, ++all;
        }
    return double(in) / double(all);
}

}  // namespace

TEST_CASE("remesh: sphere - closed genus 0 kept, vertices on the surface, edge lengths near the target", "[remesh]") {
    const TriMesh in = make_icosphere(3);  // edges ~0.15
    RemeshParams prm;
    prm.target_edge = 0.08;
    const TriMesh out = remesh_isotropic(in, prm);
    REQUIRE(build_halfedge(out).ok());
    const auto topo = analyze_topology(out);
    REQUIRE(topo.components.size() == 1);
    CHECK(topo.components[0].boundary_loops == 0);
    CHECK(topo.components[0].genus.value_or(99u) == 0);
    CHECK(out.triangles.size() > 3 * in.triangles.size());  // finer target: more triangles
    const Bvh ref(in.positions, in.triangles);
    double worst = 0.0;
    for (const Vec3& p : out.positions) worst = std::max(worst, ref.closest(p).distance);
    CHECK(worst < 1e-9);  // relaxed vertices are projected back; split midpoints lie on flat facets
    CHECK(in_band(out, 0.08) > 0.85);
}

TEST_CASE("remesh: torus keeps genus 1; a plate with holes keeps its boundaries and boundary vertices", "[remesh]") {
    const TriMesh torus = remesh_isotropic(make_torus(40, 16), RemeshParams{});
    REQUIRE(build_halfedge(torus).ok());
    CHECK(analyze_topology(torus).components[0].genus.value_or(99u) == 1);
    const TriMesh plate_in = make_grid_with_holes(16, 16, {{4, 4}, {10, 9}});
    RemeshParams prm;
    prm.target_edge = 0.05;
    const TriMesh plate = remesh_isotropic(plate_in, prm);
    REQUIRE(build_halfedge(plate).ok());
    const auto topo = analyze_topology(plate);
    REQUIRE(topo.components.size() == 1);
    CHECK(topo.components[0].boundary_loops == 3);
    CHECK(topo.components[0].genus.value_or(99u) == 0);
    for (const Vec3& p : plate.positions) CHECK(std::abs(p.z) < 1e-12);  // stays planar
}

TEST_CASE("remesh: triangle quality improves on an irregular mesh", "[remesh]") {
    // A jittered sphere: irregular valences and slivers.
    TriMesh in = make_icosphere(3);
    std::mt19937 rng(9);
    std::normal_distribution<double> g(0.0, 0.25);
    // Jitter vertices tangentially (keeps them near the sphere) to create bad triangles.
    for (auto& p : in.positions) {
        const double dx = 0.03 * g(rng), dy = 0.03 * g(rng), dz = 0.03 * g(rng);
        p.x += dx, p.y += dy, p.z += dz;
        const double r = std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        p.x /= r, p.y /= r, p.z /= r;
    }
    const double before = min_angle_q10(in);
    const TriMesh out = remesh_isotropic(in, RemeshParams{});
    REQUIRE(build_halfedge(out).ok());
    const double after = min_angle_q10(out);
    INFO("10th-percentile minimum angle: " << before << " -> " << after << " degrees");
    CHECK(after > before + 5.0);
}

TEST_CASE("remesh: aggressive coarsening stays a valid manifold of the same genus (exercises the collapse and flip guards)", "[remesh]") {
    struct Case {
        TriMesh mesh;
        double target;
        std::uint32_t genus;
    };
    for (const auto& c : {Case{make_icosphere(3), 0.9, 0u}, Case{make_torus(40, 16, 1.0, 0.3), 0.55, 1u}, Case{make_torus(48, 20, 1.0, 0.25), 0.7, 1u}}) {
        RemeshParams prm;
        prm.target_edge = c.target;
        prm.iterations = 8;
        const TriMesh out = remesh_isotropic(c.mesh, prm);
        INFO("target " << c.target << ", genus " << c.genus << ", faces " << out.triangles.size());
        REQUIRE(build_halfedge(out).ok());  // no edge on three faces, consistent orientation
        const auto topo = analyze_topology(out);
        REQUIRE(topo.components.size() == 1);
        CHECK(topo.components[0].genus.value_or(99u) == c.genus);
        CHECK(topo.components[0].boundary_loops == 0);
        for (const auto& t : out.triangles) {  // no degenerate faces
            const Vec3 &a = out.positions[t[0]], &b = out.positions[t[1]], &d = out.positions[t[2]];
            const Vec3 u{b.x - a.x, b.y - a.y, b.z - a.z}, v{d.x - a.x, d.y - a.y, d.z - a.z};
            const Vec3 n{u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x};
            CHECK(n.x * n.x + n.y * n.y + n.z * n.z > 1e-20);
        }
    }
}

TEST_CASE("remesh: boundary vertices stay on the original boundary", "[remesh]") {
    const TriMesh in = make_grid_with_holes(16, 16, {{4, 4}, {10, 9}});
    RemeshParams prm;
    prm.target_edge = 0.05;
    const TriMesh out = remesh_isotropic(in, prm);
    // Input boundary segments.
    std::vector<std::pair<Vec3, Vec3>> segs;
    {
        std::map<std::pair<std::uint32_t, std::uint32_t>, int> use;
        for (const auto& t : in.triangles)
            for (int k = 0; k < 3; ++k) {
                const auto a = t[static_cast<std::size_t>(k)], b = t[static_cast<std::size_t>((k + 1) % 3)];
                ++use[{std::min(a, b), std::max(a, b)}];
            }
        for (const auto& [e, n] : use)
            if (n == 1) segs.push_back({in.positions[e.first], in.positions[e.second]});
    }
    auto dist_to_boundary = [&](const Vec3& p) {
        double best = 1e300;
        for (const auto& [a, b] : segs) {
            const Vec3 ab{b.x - a.x, b.y - a.y, b.z - a.z}, ap{p.x - a.x, p.y - a.y, p.z - a.z};
            const double t = std::clamp((ab.x * ap.x + ab.y * ap.y + ab.z * ap.z) / (ab.x * ab.x + ab.y * ab.y + ab.z * ab.z), 0.0, 1.0);
            const Vec3 d{ap.x - t * ab.x, ap.y - t * ab.y, ap.z - t * ab.z};
            best = std::min(best, std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z));
        }
        return best;
    };
    std::map<std::pair<std::uint32_t, std::uint32_t>, int> use;
    for (const auto& t : out.triangles)
        for (int k = 0; k < 3; ++k) {
            const auto a = t[static_cast<std::size_t>(k)], b = t[static_cast<std::size_t>((k + 1) % 3)];
            ++use[{std::min(a, b), std::max(a, b)}];
        }
    std::size_t checked = 0;
    for (const auto& [e, n] : use)
        if (n == 1) {
            CHECK(dist_to_boundary(out.positions[e.first]) < 1e-12);
            CHECK(dist_to_boundary(out.positions[e.second]) < 1e-12);
            ++checked;
        }
    CHECK(checked > 50);
}

TEST_CASE("remesh: a sheared grid full of slivers becomes well shaped away from the (fixed) boundary", "[remesh]") {
    TriMesh in = make_grid(30, 30);
    for (auto& p : in.positions) p.x *= 3.0, p.y = 0.3 * p.y + 0.9 * p.x;  // sheared: angles of a few degrees
    // Quality of triangles with no boundary vertex: boundary vertices are kept fixed by design (D99), which limits the
    // triangles that touch them.
    auto interior_q10 = [](const TriMesh& m) {
        std::map<std::pair<std::uint32_t, std::uint32_t>, int> use;
        for (const auto& t : m.triangles)
            for (int k = 0; k < 3; ++k) {
                const auto a = t[static_cast<std::size_t>(k)], b = t[static_cast<std::size_t>((k + 1) % 3)];
                ++use[{std::min(a, b), std::max(a, b)}];
            }
        std::vector<bool> on_boundary(m.positions.size(), false);
        for (const auto& [e, n] : use)
            if (n == 1) on_boundary[e.first] = on_boundary[e.second] = true;
        TriMesh inner;
        inner.positions = m.positions;
        for (const auto& t : m.triangles)
            if (!on_boundary[t[0]] && !on_boundary[t[1]] && !on_boundary[t[2]]) inner.triangles.push_back(t);
        return min_angle_q10(inner);
    };
    const double before = interior_q10(in);
    const TriMesh out = remesh_isotropic(in, RemeshParams{});
    REQUIRE(build_halfedge(out).ok());
    const double after = interior_q10(out);
    INFO("interior 10th-percentile minimum angle: " << before << " -> " << after << " degrees (all triangles: "
                                                    << min_angle_q10(in) << " -> " << min_angle_q10(out) << ")");
    CHECK(before < 15.0);
    // Measured 31.8 deg on ARM and 29.5 on x86: rounding changes which collapses and flips fire. The threshold leaves
    // margin on both (D99); the claim is a large improvement from ~3 deg, not a platform-exact value.
    CHECK(after > 25.0);
}
