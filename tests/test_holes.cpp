#include <algorithm>
#include <cmath>
#include <functional>
#include <random>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/generate.h"
#include "core/halfedge.h"
#include "core/holes.h"
#include "core/topology.h"

using namespace dmw;
using Catch::Matchers::WithinAbs;
using Tri = std::array<std::uint32_t, 3>;

namespace {

// Remove faces whose centroid satisfies `cut`, then drop unused vertices.
TriMesh punch(const TriMesh& m, const std::function<bool(const Vec3&)>& cut) {
    TriMesh out;
    std::vector<std::uint32_t> remap(m.positions.size(), 0xFFFFFFFFu);
    for (Tri t : m.triangles) {
        const Vec3 c{(m.positions[t[0]].x + m.positions[t[1]].x + m.positions[t[2]].x) / 3.0,
                     (m.positions[t[0]].y + m.positions[t[1]].y + m.positions[t[2]].y) / 3.0,
                     (m.positions[t[0]].z + m.positions[t[1]].z + m.positions[t[2]].z) / 3.0};
        if (cut(c)) continue;
        for (auto& v : t) {
            if (remap[v] == 0xFFFFFFFFu) remap[v] = static_cast<std::uint32_t>(out.positions.size()), out.positions.push_back(m.positions[v]);
            v = remap[v];
        }
        out.triangles.push_back(t);
    }
    return out;
}

// All triangulations of the polygon i..j (chord (i, j)), as lists of triangles (i, j, m).
void all_triangulations(std::size_t i, std::size_t j, std::vector<std::vector<Tri>>& out) {
    if (j < i + 2) {
        out.push_back({});
        return;
    }
    for (std::size_t m = i + 1; m < j; ++m) {
        std::vector<std::vector<Tri>> left, right;
        all_triangulations(i, m, left), all_triangulations(m, j, right);
        for (const auto& l : left)
            for (const auto& r : right) {
                std::vector<Tri> t = l;
                t.insert(t.end(), r.begin(), r.end());
                t.push_back({static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(j), static_cast<std::uint32_t>(m)});
                out.push_back(std::move(t));
            }
    }
}

}  // namespace

TEST_CASE("holes: boundary loops of meshes with known boundaries", "[holes]") {
    CHECK(boundary_loops(make_icosphere(2)).empty());
    CHECK(boundary_loops(make_cylinder(12, 4)).size() == 2);
    CHECK(boundary_loops(make_grid_with_holes(8, 8, {{2, 2}, {5, 4}})).size() == 3);  // outer rim + 2 holes
}

TEST_CASE("holes: DP is exactly optimal for area; for (dihedral, area) it is Liepa's approximation (gap measured)", "[holes]") {
    std::mt19937 rng(17);
    std::uniform_real_distribution<double> jitter(-0.25, 0.25);
    int exact_dihedral = 0, trials = 0;
    double worst_gap = 0.0;
    for (int trial = 0; trial < 120; ++trial) {
        const std::size_t n = 4 + static_cast<std::size_t>(trial % 5);  // 4..8 vertices: up to 132 triangulations
        std::vector<Vec3> loop, outside;
        for (std::size_t i = 0; i < n; ++i) {
            const double a = -2.0 * std::acos(-1.0) * double(i) / double(n);  // clockwise: the hole lies to the right
            loop.push_back({std::cos(a) + jitter(rng), std::sin(a) + jitter(rng), jitter(rng)});
        }
        for (std::size_t i = 0; i < n; ++i) {  // rim faces outside the loop, sloping down
            const Vec3 &a = loop[i], &b = loop[(i + 1) % n];
            outside.push_back({1.6 * (a.x + b.x) / 2, 1.6 * (a.y + b.y) / 2, -0.3 + jitter(rng)});
        }
        std::vector<std::vector<Tri>> all;
        all_triangulations(0, n - 1, all);
        // Area only: decomposable, so the DP must hit the exhaustive minimum exactly.
        const auto area_dp = triangulate_loop(loop, outside, false);
        REQUIRE(area_dp.size() == n - 2);
        double min_area = 1e300;
        for (const auto& t : all) min_area = std::min(min_area, triangulation_weight(loop, outside, t).second);
        INFO("trial " << trial << ", n " << n);
        CHECK_THAT(triangulation_weight(loop, outside, area_dp).second, WithinAbs(min_area, 1e-12));
        // Dihedral first: record how often the DP matches the exhaustive optimum and the largest angle gap.
        const auto wdp = triangulation_weight(loop, outside, triangulate_loop(loop, outside, true));
        double best_angle = 1e300;
        for (const auto& t : all) best_angle = std::min(best_angle, triangulation_weight(loop, outside, t).first);
        CHECK(wdp.first >= best_angle - 1e-12);  // sanity: never better than the true optimum
        exact_dihedral += wdp.first <= best_angle + 1e-9;
        worst_gap = std::max(worst_gap, wdp.first - best_angle);
        ++trials;
    }
    WARN("dihedral DP matched the exhaustive optimum on " << exact_dihedral << " of " << trials << " polygons; largest gap "
                                                           << worst_gap * 180.0 / std::acos(-1.0) << " degrees");
    CHECK(exact_dihedral * 10 >= trials * 8);  // documented property: an approximation that is usually exact
}

TEST_CASE("holes: a planar hole fills flat; the result is closed off and consistently oriented", "[holes]") {
    // Grid with a square hole; the outer rim is longer than max_loop and stays open.
    const auto grid = punch(make_grid(30, 30), [](const Vec3& c) { return std::abs(c.x - 0.5) < 0.2 && std::abs(c.y - 0.5) < 0.2; });
    HoleFillParams prm;
    prm.max_loop = 100;
    const auto r = fill_holes(grid, prm);
    CHECK(r.loops == 2);
    CHECK(r.filled == 1);
    CHECK_FALSE(r.patch_vertices.empty());  // refinement added vertices
    for (auto v : r.patch_vertices) CHECK_THAT(r.mesh.positions[v].z, WithinAbs(0.0, 1e-9));  // fairing keeps it planar
    CHECK(build_halfedge(r.mesh).ok());
    CHECK(analyze_topology(r.mesh).components[0].boundary_loops == 1);
}

TEST_CASE("holes: a sphere with a cap removed closes to genus 0 and the patch follows the sphere", "[holes]") {
    // Off-centre on purpose: a sphere at the origin is symmetric under p -> -p, which would hide a sign error.
    const Vec3 centre{2.0, 3.0, 1.0};
    TriMesh sphere = make_icosphere(4);
    for (auto& p : sphere.positions) p.x += centre.x, p.y += centre.y, p.z += centre.z;
    const auto holed = punch(sphere, [&](const Vec3& c) { return c.z - centre.z > 0.8; });
    const auto r = fill_holes(holed);
    REQUIRE(r.filled == 1);
    const auto topo = analyze_topology(r.mesh);
    REQUIRE(topo.components.size() == 1);
    CHECK(topo.components[0].boundary_loops == 0);
    CHECK(topo.components[0].genus.value_or(99u) == 0);
    CHECK(build_halfedge(r.mesh).ok());
    // Thin-plate fairing continues the curvature: patch vertices sit near the unit sphere (a flat cap would be 0.2 off).
    double worst = 0.0;
    for (auto v : r.patch_vertices) {
        const Vec3 p{r.mesh.positions[v].x - centre.x, r.mesh.positions[v].y - centre.y, r.mesh.positions[v].z - centre.z};
        worst = std::max(worst, std::abs(std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z) - 1.0));
        CHECK(p.z > 0.7);  // inside the removed cap
    }
    INFO("largest radial deviation of patch vertices: " << worst);
    CHECK(worst < 0.05);
}

TEST_CASE("holes: a chord that already exists as a mesh edge is never used (the result stays manifold)", "[holes]") {
    // A square hole whose two opposite corners are already joined by an edge outside the hole (a thin "bridge"
    // triangle fan): the obvious triangulation would use that diagonal and put three faces on one edge.
    const std::vector<Vec3> loop{{0, 0, 0}, {0, 1, 0}, {1, 1, 0}, {1, 0, 0}};  // clockwise seen from +z
    const std::vector<Vec3> outside{{-1, 0.5, 0}, {0.5, 2, 0}, {2, 0.5, 0}, {0.5, -1, 0}};
    // Forbid the chord (0, 2): only the other diagonal (1, 3) remains.
    const auto tris = triangulate_loop(loop, outside, true, [](std::uint32_t a, std::uint32_t b) { return !(a == 0 && b == 2); });
    REQUIRE(tris.size() == 2);
    for (const auto& t : tris) {
        const bool has0 = std::find(t.begin(), t.end(), 0u) != t.end(), has2 = std::find(t.begin(), t.end(), 2u) != t.end();
        CHECK_FALSE((has0 && has2));
    }
    // Forbid both diagonals: no triangulation exists, and none is returned.
    CHECK(triangulate_loop(loop, outside, true, [](std::uint32_t, std::uint32_t) { return false; }).empty());
}
