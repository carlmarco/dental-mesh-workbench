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
