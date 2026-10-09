#include <cmath>
#include <random>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/bvh.h"
#include "core/generate.h"

using namespace dmw;
using Catch::Matchers::WithinAbs;

namespace {

Vec3 random_unit(std::mt19937& rng) {
    std::normal_distribution<double> g(0.0, 1.0);
    const Vec3 v{g(rng), g(rng), g(rng)};
    const double n = std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    return {v.x / n, v.y / n, v.z / n};
}

void compare(const TriMesh& m, std::mt19937& rng, int rays, double spread) {
    const Bvh bvh(m.positions, m.triangles);
    std::uniform_real_distribution<double> u(-spread, spread);
    int hits = 0;
    for (int i = 0; i < rays; ++i) {
        const Vec3 o{u(rng), u(rng), u(rng)}, d = random_unit(rng);
        const RayHit a = bvh.intersect(o, d), b = intersect_brute_force(m.positions, m.triangles, o, d);
        INFO("ray " << i);
        REQUIRE(a.hit == b.hit);
        CHECK(bvh.occluded(o, d) == b.hit);
        if (!b.hit) continue;
        ++hits;
        CHECK_THAT(a.t, WithinAbs(b.t, 1e-12));
        if (a.face != b.face) CHECK_THAT(a.t, WithinAbs(b.t, 1e-12));  // a tie on a shared edge: same distance
    }
    CHECK(hits > rays / 10);  // the test exercises hits, not just misses
}

}  // namespace

TEST_CASE("bvh: one triangle - analytic hit, barycentrics, misses, parallel ray, range and ignore", "[bvh]") {
    const std::vector<Vec3> p{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    const std::vector<std::array<std::uint32_t, 3>> t{{0, 1, 2}};
    const Bvh bvh(p, t);
    const RayHit h = bvh.intersect({0.25, 0.5, 2.0}, {0, 0, -1});
    REQUIRE(h.hit);
    CHECK_THAT(h.t, WithinAbs(2.0, 1e-15));
    CHECK_THAT(h.u, WithinAbs(0.25, 1e-15));
    CHECK_THAT(h.v, WithinAbs(0.5, 1e-15));
    CHECK_FALSE(bvh.intersect({0.8, 0.8, 2.0}, {0, 0, -1}).hit);   // outside the triangle
    CHECK_FALSE(bvh.intersect({0.2, 0.2, 1.0}, {1, 0, 0}).hit);    // parallel to its plane
    CHECK_FALSE(bvh.intersect({0.2, 0.2, 2.0}, {0, 0, 1}).hit);    // pointing away
    CHECK_FALSE(bvh.intersect({0.2, 0.2, 2.0}, {0, 0, -1}, 0.0, 1.5).hit);  // beyond t_max
    CHECK_FALSE(bvh.intersect({0.2, 0.2, 2.0}, {0, 0, -1}, 0.0, 10.0, 0).hit);  // the only face ignored
}

TEST_CASE("bvh: matches brute force on random rays (sphere, torus, random triangle soup)", "[bvh]") {
    std::mt19937 rng(10);
    compare(make_icosphere(3), rng, 2000, 1.5);
    compare(make_torus(40, 16, 1.0, 0.35), rng, 2000, 1.5);
    TriMesh soup;  // overlapping, randomly oriented triangles of mixed size: stresses the SAH split and traversal
    std::uniform_real_distribution<double> u(-1.0, 1.0), s(0.01, 0.4);
    for (int i = 0; i < 3000; ++i) {
        const Vec3 c{u(rng), u(rng), u(rng)};
        const double r = s(rng);
        for (int k = 0; k < 3; ++k) {
            const Vec3 d = random_unit(rng);
            soup.positions.push_back({c.x + r * d.x, c.y + r * d.y, c.z + r * d.z});
        }
        const auto b = static_cast<std::uint32_t>(3 * i);
        soup.triangles.push_back({b, b + 1, b + 2});
    }
    compare(soup, rng, 3000, 1.2);
}

TEST_CASE("bvh: rays from inside a closed sphere always hit, at the radius", "[bvh]") {
    const auto m = make_icosphere(4);  // unit sphere, flat facets: hit distance within (inradius, 1]
    const Bvh bvh(m.positions, m.triangles);
    std::mt19937 rng(3);
    for (int i = 0; i < 1000; ++i) {
        const RayHit h = bvh.intersect({0, 0, 0}, random_unit(rng));
        REQUIRE(h.hit);
        CHECK(h.t <= 1.0 + 1e-12);
        CHECK(h.t > 0.99);
    }
    CHECK(bvh.node_count() < 2 * m.triangles.size());
}

TEST_CASE("bvh: axis-aligned rays from grid-aligned origins (zero direction components, origins on box planes)", "[bvh]") {
    // Two parallel grids: rays straight down from the upper grid's vertices start exactly on bounding-box planes
    // of the lower grid's nodes and pass exactly through its vertices.
    TriMesh m = make_grid(16, 16);
    TriMesh top = make_grid(16, 16);
    for (auto& p : top.positions) p.z = 1.0;
    append(m, top);
    const Bvh bvh(m.positions, m.triangles);
    for (std::size_t v = m.positions.size() / 2; v < m.positions.size(); ++v) {
        const Vec3 o = m.positions[v];
        const RayHit a = bvh.intersect(o, {0, 0, -1}, 1e-9), b = intersect_brute_force(m.positions, m.triangles, o, {0, 0, -1}, 1e-9);
        INFO("vertex " << v);
        REQUIRE(b.hit);
        REQUIRE(a.hit);
        CHECK_THAT(a.t, WithinAbs(1.0, 1e-12));
    }
}

TEST_CASE("bvh: closest point matches brute force (distance) and lies on the mesh", "[bvh]") {
    std::mt19937 rng(23);
    std::uniform_real_distribution<double> u(-2.0, 2.0);
    for (const TriMesh& m : {make_icosphere(3), make_torus(30, 12, 1.0, 0.35)}) {
        const Bvh bvh(m.positions, m.triangles);
        for (int i = 0; i < 1500; ++i) {
            const Vec3 p{u(rng), u(rng), u(rng)};
            const ClosestPoint a = bvh.closest(p), b = closest_brute_force(m.positions, m.triangles, p);
            INFO("query " << i);
            CHECK_THAT(a.distance, WithinAbs(b.distance, 1e-12));
            const Vec3 d{a.point.x - p.x, a.point.y - p.y, a.point.z - p.z};
            CHECK_THAT(std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z), WithinAbs(a.distance, 1e-12));
        }
    }
}

TEST_CASE("bvh: closest point on a triangle in each Voronoi region", "[bvh]") {
    const Vec3 a{0, 0, 0}, b{1, 0, 0}, c{0, 1, 0};
    auto near = [](const Vec3& x, const Vec3& y) { return std::abs(x.x - y.x) + std::abs(x.y - y.y) + std::abs(x.z - y.z) < 1e-12; };
    CHECK(near(closest_point_on_triangle({-1, -1, 2}, a, b, c), a));           // vertex a
    CHECK(near(closest_point_on_triangle({2, -0.5, 1}, a, b, c), b));          // vertex b
    CHECK(near(closest_point_on_triangle({-0.5, 2, -1}, a, b, c), c));         // vertex c
    CHECK(near(closest_point_on_triangle({0.5, -1, 0}, a, b, c), {0.5, 0, 0}));  // edge ab
    CHECK(near(closest_point_on_triangle({-1, 0.3, 0}, a, b, c), {0, 0.3, 0}));  // edge ac
    CHECK(near(closest_point_on_triangle({1, 1, 0}, a, b, c), {0.5, 0.5, 0}));   // edge bc
    CHECK(near(closest_point_on_triangle({0.2, 0.3, 5}, a, b, c), {0.2, 0.3, 0}));  // face
}
