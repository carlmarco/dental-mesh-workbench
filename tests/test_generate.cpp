#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <set>
#include <utility>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/generate.h"

using namespace dmw;
using Catch::Matchers::WithinAbs;

namespace {

// Independent of the topology module on purpose: generators are checked against
// closed-form counts with the simplest possible edge counter.
std::size_t count_edges(const TriMesh& m) {
    std::set<std::pair<std::uint32_t, std::uint32_t>> edges;
    for (const auto& t : m.triangles) {
        for (int k = 0; k < 3; ++k) {
            const std::uint32_t a = t[static_cast<std::size_t>(k)];
            const std::uint32_t b = t[static_cast<std::size_t>((k + 1) % 3)];
            edges.insert({std::min(a, b), std::max(a, b)});
        }
    }
    return edges.size();
}

// (1/6) sum v0 . (v1 x v2): the enclosed volume for an outward-wound closed mesh
// (divergence theorem), negative if wound inward.
double signed_volume(const TriMesh& m) {
    double vol = 0.0;
    for (const auto& t : m.triangles) {
        const Vec3 &a = m.positions[t[0]], &b = m.positions[t[1]], &c = m.positions[t[2]];
        vol += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x) +
               a.z * (b.x * c.y - b.y * c.x);
    }
    return vol / 6.0;
}

bool indices_in_range(const TriMesh& m) {
    for (const auto& t : m.triangles) {
        for (std::uint32_t v : t) {
            if (v >= m.positions.size()) return false;
        }
    }
    return true;
}

}  // namespace

TEST_CASE("generate: grid counts and +z winding", "[generate]") {
    const auto m = make_grid(4, 3);
    REQUIRE(indices_in_range(m));
    CHECK(m.positions.size() == 5 * 4);
    CHECK(m.triangles.size() == 2 * 4 * 3);
    CHECK(count_edges(m) == 4 * 4 + 3 * 5 + 4 * 3);  // horizontal + vertical + diagonals
    // Every triangle's normal (b-a)x(c-a) has positive z.
    for (const auto& t : m.triangles) {
        const Vec3 &a = m.positions[t[0]], &b = m.positions[t[1]], &c = m.positions[t[2]];
        CHECK((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x) > 0.0);
    }
}

TEST_CASE("generate: grid with holes removes 2 faces and 1 edge per hole", "[generate]") {
    const auto m = make_grid_with_holes(6, 6, {{1, 1}, {4, 3}});
    REQUIRE(indices_in_range(m));
    CHECK(m.positions.size() == 7 * 7);  // corners of interior holes stay referenced
    CHECK(m.triangles.size() == 2 * 36 - 4);
    CHECK(count_edges(m) == (6 * 7 + 6 * 7 + 36) - 2);  // only the two diagonals vanish
}

TEST_CASE("generate: torus counts, outward winding, volume", "[generate]") {
    const auto m = make_torus(16, 8, 1.0, 0.3);
    REQUIRE(indices_in_range(m));
    CHECK(m.positions.size() == 16 * 8);
    CHECK(m.triangles.size() == 2 * 16 * 8);
    CHECK(count_edges(m) == 3 * 16 * 8);
    // Positive = outward. Sanity bound only: within 15% of the smooth 2 pi^2 R r^2. (A torus
    // is not convex, so "inscribed polyhedron has less volume" does not hold in general.)
    const double vol = signed_volume(m);
    const double smooth = 2.0 * std::numbers::pi * std::numbers::pi * 1.0 * 0.3 * 0.3;
    CHECK(vol > 0.85 * smooth);
    CHECK(vol < 1.15 * smooth);
}

TEST_CASE("generate: cylinder counts and outward winding", "[generate]") {
    const auto m = make_cylinder(12, 4, 0.5, 2.0);
    REQUIRE(indices_in_range(m));
    CHECK(m.positions.size() == 12 * 5);
    CHECK(m.triangles.size() == 2 * 12 * 4);
    // rings nu(nv+1) + verticals nu*nv + diagonals nu*nv = nu(3nv+1); chi = V - E + F = 0.
    CHECK(count_edges(m) == 12 * (3 * 4 + 1));
    // Outward: each face normal points away from the z axis.
    for (const auto& t : m.triangles) {
        const Vec3 &a = m.positions[t[0]], &b = m.positions[t[1]], &c = m.positions[t[2]];
        const double nx = (b.y - a.y) * (c.z - a.z) - (b.z - a.z) * (c.y - a.y);
        const double ny = (b.z - a.z) * (c.x - a.x) - (b.x - a.x) * (c.z - a.z);
        CHECK(nx * (a.x + b.x + c.x) + ny * (a.y + b.y + c.y) > 0.0);
    }
}

TEST_CASE("generate: icosphere counts, on the sphere, outward", "[generate]") {
    for (std::uint32_t s = 0; s <= 3; ++s) {
        const auto m = make_icosphere(s, 2.0);
        const std::size_t p = std::size_t{1} << (2 * s);  // 4^s
        REQUIRE(indices_in_range(m));
        CHECK(m.positions.size() == 10 * p + 2);
        CHECK(m.triangles.size() == 20 * p);
        CHECK(count_edges(m) == 30 * p);
        for (const auto& v : m.positions) {
            CHECK_THAT(std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z), WithinAbs(2.0, 1e-12));
        }
        CHECK(signed_volume(m) > 0.0);
    }
}

TEST_CASE("generate: Moebius strip counts", "[generate]") {
    const auto m = make_mobius(8);
    REQUIRE(indices_in_range(m));
    CHECK(m.positions.size() == 16);
    CHECK(m.triangles.size() == 16);
    CHECK(count_edges(m) == 32);  // n rungs + n diagonals + 2n boundary
}

TEST_CASE("generate: append is a disjoint union", "[generate]") {
    auto m = make_grid(1, 1);
    append(m, make_grid(1, 1));
    REQUIRE(indices_in_range(m));
    CHECK(m.positions.size() == 8);
    CHECK(m.triangles.size() == 4);
    CHECK(m.triangles[2][0] == 4);  // second copy's indices offset by 4
}

TEST_CASE("generate: invalid parameters give an empty mesh", "[generate]") {
    CHECK(make_grid(0, 3).triangles.empty());
    CHECK(make_torus(2, 8).triangles.empty());
    CHECK(make_icosphere(11).triangles.empty());
    CHECK(make_mobius(2).triangles.empty());
    CHECK(make_cylinder(2, 3).triangles.empty());
    CHECK(make_grid_with_holes(3, 3, {{3, 0}}).triangles.empty());  // cell out of range
}
