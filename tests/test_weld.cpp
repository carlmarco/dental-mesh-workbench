#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "core/weld.h"

using dmw::TriMesh;
using dmw::Vec3;
using dmw::weld_vertices;
using Tris = std::vector<std::array<std::uint32_t, 3>>;

namespace {

// Unwelded "soup": every triangle owns its 3 vertices, as STL stores them.
TriMesh soup(const std::vector<std::array<Vec3, 3>>& tris) {
    TriMesh m;
    for (const auto& t : tris) {
        const auto base = static_cast<std::uint32_t>(m.positions.size());
        m.positions.insert(m.positions.end(), t.begin(), t.end());
        m.triangles.push_back({base, base + 1, base + 2});
    }
    return m;
}

constexpr Vec3 O{0, 0, 0}, X{1, 0, 0}, Y{0, 1, 0}, Z{0, 0, 1};

}  // namespace

TEST_CASE("weld exact: tetrahedron soup (12 vertices) -> 4, first-appearance order", "[weld]") {
    // Same outward-wound faces as the OBJ tetrahedron fixture.
    const auto m = weld_vertices(soup({{O, Y, X}, {O, X, Z}, {O, Z, Y}, {X, Y, Z}}));
    REQUIRE(m.positions.size() == 4);
    CHECK(m.positions[1].y == 1.0);  // Y appeared second
    CHECK(m.triangles == Tris{{0, 1, 2}, {0, 2, 3}, {0, 3, 1}, {2, 1, 3}});
}

TEST_CASE("weld exact: nearly-equal points stay separate", "[weld]") {
    const Vec3 X_next{std::nextafter(1.0, 2.0), 0, 0};  // 1 ulp away from X
    const auto m = weld_vertices(soup({{O, X, Y}, {O, X_next, Y}}));
    CHECK(m.positions.size() == 4);
}

TEST_CASE("weld exact: -0.0 and +0.0 are the same point", "[weld]") {
    // Hashing raw bits would treat these as different: the sign bit differs.
    const auto m = weld_vertices(soup({{O, X, Y}, {Vec3{-0.0, 0, 0}, X, Y}}));
    CHECK(m.positions.size() == 3);
}

TEST_CASE("weld epsilon: points straddling a grid-cell boundary still merge", "[weld]") {
    // eps = 0.1 -> cell boundary at x = 0.1. A and B are 0.002 apart but in
    // different cells; a grid that only checks the point's own cell misses this.
    const Vec3 A{0.099, 0, 0}, B{0.101, 0, 0};
    const auto m = weld_vertices(soup({{A, X, Y}, {B, X, Y}}), 0.1);
    REQUIRE(m.positions.size() == 3);
    CHECK(m.positions[0].x == 0.099);  // representative keeps the first point's position
    CHECK(m.triangles == Tris{{0, 1, 2}, {0, 1, 2}});
}

TEST_CASE("weld epsilon: deterministic first-representative rule, not transitive", "[weld]") {
    // eps = 1: B (0.8 from A) merges into A. C is 0.8 from B but 1.6 from A,
    // so it is NOT merged: B is no longer a representative.
    const Vec3 A{0, 0, 0}, B{0.8, 0, 0}, C{1.6, 0, 0}, P{10, 0, 0}, Q{0, 10, 0};
    const auto m = weld_vertices(soup({{A, P, Q}, {B, P, Q}, {C, P, Q}}), 1.0);
    REQUIRE(m.positions.size() == 4);
    CHECK(m.triangles == Tris{{0, 1, 2}, {0, 1, 2}, {3, 1, 2}});
}

TEST_CASE("weld never drops triangles, even ones it collapses", "[weld]") {
    const auto m = weld_vertices(soup({{O, Vec3{0.1, 0, 0}, Y}}), 0.5);
    REQUIRE(m.triangles.size() == 1);
    CHECK(m.triangles[0][0] == m.triangles[0][1]);  // degenerate, left for diagnostics
}

TEST_CASE("weld: empty mesh stays empty", "[weld]") {
    const auto m = weld_vertices(TriMesh{});
    CHECK(m.positions.empty());
    CHECK(m.triangles.empty());
}
