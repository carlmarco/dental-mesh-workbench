#include <algorithm>
#include <cmath>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/bvh.h"
#include "core/generate.h"
#include "core/thickness.h"

using namespace dmw;
using Catch::Matchers::WithinAbs;

namespace {

// Closed slab [-2, 2] x [-2, 2] x [0, h] built from two grids and four side strips, outward normals.
TriMesh slab(double h, std::uint32_t n = 40) {
    TriMesh top = make_grid(n, n), bottom = make_grid(n, n);
    for (auto& p : top.positions) p.x = 4.0 * p.x - 2.0, p.y = 4.0 * p.y - 2.0, p.z = h;
    for (auto& p : bottom.positions) p.x = 4.0 * p.x - 2.0, p.y = 4.0 * p.y - 2.0, p.z = 0.0;
    for (auto& t : bottom.triangles) std::swap(t[1], t[2]);  // face down
    TriMesh m = top;
    append(m, bottom);
    return m;  // the side walls are not needed for interior vertices (rays never reach them)
}

// Spherical shell: outer unit sphere (outward normals) and inner sphere of radius r with normals flipped (pointing
// out of the solid, i.e. towards the centre).
TriMesh shell(double r, int level = 4) {
    TriMesh m = make_icosphere(level), inner = make_icosphere(level);
    for (auto& p : inner.positions) p.x *= r, p.y *= r, p.z *= r;
    for (auto& t : inner.triangles) std::swap(t[1], t[2]);
    append(m, inner);
    return m;
}

}  // namespace

TEST_CASE("thickness: slab - the normal ray and the cone minimum give the slab thickness", "[thickness]") {
    const double h = 0.8;
    const auto m = slab(h);
    const Bvh bvh(m.positions, m.triangles);
    const auto t = wall_thickness(m, bvh);
    std::size_t checked = 0;
    for (std::size_t v = 0; v < m.positions.size(); ++v) {
        const Vec3& p = m.positions[v];
        if (std::abs(p.x) > 1.0 || std::abs(p.y) > 1.0) continue;  // interior: cones stay inside the slab footprint
        CHECK_THAT(t.along_normal[v], WithinAbs(h, 1e-9));
        CHECK_THAT(t.cone_min[v], WithinAbs(h, 1e-9));
        CHECK(t.cone_median[v] >= h - 1e-9);                       // oblique rays are longer: h / cos(angle)
        CHECK(t.cone_median[v] <= h / std::cos(30.0 * std::acos(-1.0) / 180.0) + 1e-9);
        ++checked;
    }
    CHECK(checked > 100);
}

TEST_CASE("thickness: spherical shell - thickness R_out - R_in from both walls", "[thickness]") {
    const double r = 0.7;
    const auto m = shell(r);
    const Bvh bvh(m.positions, m.triangles);
    const auto t = wall_thickness(m, bvh);
    double worst = 0.0;
    for (std::size_t v = 0; v < m.positions.size(); ++v) {
        REQUIRE(std::isfinite(t.along_normal[v]));
        worst = std::max(worst, std::abs(t.along_normal[v] - (1.0 - r)));
        CHECK(t.cone_min[v] <= t.along_normal[v] + 1e-12);
    }
    CHECK(worst < 0.01);  // flat facets of a level-4 icosphere: the radius varies by < 0.5%
}

TEST_CASE("thickness: an open surface has no thickness (NaN), a single sheet facing the wrong way neither", "[thickness]") {
    const auto m = make_grid(10, 10);  // one open sheet
    const Bvh bvh(m.positions, m.triangles);
    const auto t = wall_thickness(m, bvh);
    for (double x : t.along_normal) CHECK(std::isnan(x));
    for (double x : t.cone_median) CHECK(std::isnan(x));
}

TEST_CASE("thickness: a wall hit from the wrong side is not a thickness (two sheets facing the same way)", "[thickness]") {
    // Two parallel sheets both facing +z do not bound a solid: the downward ray from the upper sheet meets the lower
    // sheet's front side (its normal opposes the ray), which must be rejected.
    TriMesh m = make_grid(12, 12);
    TriMesh upper = make_grid(12, 12);
    for (auto& p : upper.positions) p.x += 0.013, p.y += 0.017, p.z = 1.0;  // offset: rays avoid vertices
    const std::size_t lower_vertices = m.positions.size();
    append(m, upper);
    const Bvh bvh(m.positions, m.triangles);
    const auto t = wall_thickness(m, bvh);
    for (std::size_t v = lower_vertices; v < m.positions.size(); ++v) CHECK(std::isnan(t.along_normal[v]));
}
