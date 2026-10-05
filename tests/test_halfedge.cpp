#include <algorithm>
#include <array>
#include <cstdint>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "core/halfedge.h"

using dmw::build_halfedge;
using dmw::dest;
using dmw::HalfEdgeMesh;
using dmw::kInvalid;
using dmw::one_ring;
using dmw::TriMesh;
using Ring = std::vector<std::uint32_t>;

namespace {

TriMesh make(std::vector<dmw::Vec3> p, std::vector<std::array<std::uint32_t, 3>> t) {
    return TriMesh{std::move(p), std::move(t)};
}

// Unit tetrahedron, outward-wound (same as the OBJ/weld fixtures).
TriMesh tetra() {
    return make({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}},
                {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}});
}

// Unit cube, each quad fan-split as in the OBJ cube test. V=8, E=18, F=12.
TriMesh cube() {
    return make({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
                {{0, 3, 2}, {0, 2, 1}, {4, 5, 6}, {4, 6, 7}, {0, 1, 5}, {0, 5, 4},
                 {2, 3, 7}, {2, 7, 6}, {0, 4, 7}, {0, 7, 3}, {1, 2, 6}, {1, 6, 5}});
}

// Hand-built open mesh: 4 triangles fanned around vertex 0 in the z=0 plane,
// wound CCW seen from +z. Interior vertex 0; boundary square 1-2-3-4.
//        2
//      / | \
//     3--0--1
//      \ | /
//        4
TriMesh open_fan() {
    return make({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {-1, 0, 0}, {0, -1, 0}},
                {{0, 1, 2}, {0, 2, 3}, {0, 3, 4}, {0, 4, 1}});
}

// Checks every invariant documented in halfedge.h.
void check_invariants(const HalfEdgeMesh& m, const TriMesh& src) {
    const auto nh = static_cast<std::uint32_t>(m.origin.size());
    REQUIRE(nh == 3 * src.triangles.size());
    REQUIRE(m.twin.size() == nh);
    REQUIRE(m.vertex_halfedge.size() == src.positions.size());
    for (std::uint32_t h = 0; h < nh; ++h) {
        CHECK(dmw::next(dmw::next(dmw::next(h))) == h);
        CHECK(dmw::prev(dmw::next(h)) == h);
        CHECK(dmw::face(h) == h / 3);
        CHECK(m.origin[h] == src.triangles[h / 3][h % 3]);  // the numbering convention
        const std::uint32_t t = m.twin[h];
        if (t == kInvalid) continue;
        CHECK(t != h);
        CHECK(m.twin[t] == h);
        CHECK(m.origin[t] == dest(m, h));
        CHECK(dest(m, t) == m.origin[h]);
    }
    for (std::uint32_t v = 0; v < m.vertex_halfedge.size(); ++v) {
        const std::uint32_t h = m.vertex_halfedge[v];
        if (h != kInvalid) CHECK(m.origin[h] == v);
    }
}

std::size_t count_boundary(const HalfEdgeMesh& m) {
    return static_cast<std::size_t>(std::count(m.twin.begin(), m.twin.end(), kInvalid));
}

// Rotates a cyclic sequence to start at `first`, so cyclic orders compare with ==.
Ring starting_at(Ring r, std::uint32_t first) {
    const auto it = std::find(r.begin(), r.end(), first);
    if (it != r.end()) std::rotate(r.begin(), it, r.end());
    return r;
}

}  // namespace

TEST_CASE("half-edge: index arithmetic within a face", "[halfedge]") {
    CHECK(dmw::next(0) == 1);
    CHECK(dmw::next(2) == 0);  // wraps within face 0, never into face 1
    CHECK(dmw::prev(3) == 5);
    CHECK(dmw::face(5) == 1);
}

TEST_CASE("half-edge: tetrahedron is closed and satisfies all invariants", "[halfedge]") {
    const auto src = tetra();
    const auto r = build_halfedge(src);
    REQUIRE(r.ok());
    check_invariants(r.mesh, src);
    CHECK(count_boundary(r.mesh) == 0);  // closed surface: every half-edge has a twin
}

TEST_CASE("half-edge: tetrahedron one-ring is counter-clockwise from outside", "[halfedge]") {
    // At vertex 0 the outward direction is n = -(1,1,1). (Y - O) x (X - O) = -z and
    // -z . n > 0, so Y -> X is a CCW step about n: cyclic order (2, 1, 3).
    const auto r = build_halfedge(tetra());
    REQUIRE(r.ok());
    const Ring ring = one_ring(r.mesh, 0);
    REQUIRE(ring.size() == 3);
    CHECK(starting_at(ring, 2) == Ring{2, 1, 3});
}

TEST_CASE("half-edge: cube - closed, Euler characteristic 2, valences sum to 2E", "[halfedge]") {
    const auto src = cube();
    const auto r = build_halfedge(src);
    REQUIRE(r.ok());
    check_invariants(r.mesh, src);
    CHECK(count_boundary(r.mesh) == 0);
    // Closed: each edge has exactly two half-edges, so E = H / 2.
    const std::size_t E = r.mesh.origin.size() / 2;
    CHECK(8 - static_cast<long>(E) + 12 == 2);
    std::size_t valence_sum = 0;
    for (std::uint32_t v = 0; v < 8; ++v) valence_sum += one_ring(r.mesh, v).size();
    CHECK(valence_sum == 2 * E);  // handshake lemma
}

TEST_CASE("half-edge: open fan - boundary half-edges and interior one-ring", "[halfedge]") {
    const auto src = open_fan();
    const auto r = build_halfedge(src);
    REQUIRE(r.ok());
    check_invariants(r.mesh, src);
    CHECK(count_boundary(r.mesh) == 4);  // the outer square: 1->2, 2->3, 3->4, 4->1
    // Interior vertex: CCW about +z starting anywhere.
    CHECK(starting_at(one_ring(r.mesh, 0), 1) == Ring{1, 2, 3, 4});
}

TEST_CASE("half-edge: boundary vertex one-ring covers the whole open fan", "[halfedge]") {
    // Vertex 1 = (1,0,0) has 2 faces, so 3 neighbors. CCW about +z, by angle as seen
    // from vertex 1: 2 at 135 deg, 0 at 180 deg, 4 at 225 deg. A boundary ring is not
    // cyclic: it must start at one boundary neighbor (2, across 1->2, whose twin is
    // invalid) and end at the other (4). Starting anywhere else misses neighbors.
    const auto r = build_halfedge(open_fan());
    REQUIRE(r.ok());
    const std::uint32_t h = r.mesh.vertex_halfedge[1];
    REQUIRE(h != kInvalid);
    CHECK(r.mesh.twin[h] == kInvalid);  // boundary vertex stores its boundary half-edge
    CHECK(one_ring(r.mesh, 1) == Ring{2, 0, 4});
}

TEST_CASE("half-edge: isolated vertex is allowed and has an empty ring", "[halfedge]") {
    const auto r = build_halfedge(make({{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {5, 5, 5}}, {{0, 1, 2}}));
    REQUIRE(r.ok());
    CHECK(r.mesh.vertex_halfedge[3] == kInvalid);
    CHECK(one_ring(r.mesh, 3).empty());
}

TEST_CASE("half-edge: empty mesh builds", "[halfedge]") {
    const auto r = build_halfedge(TriMesh{});
    REQUIRE(r.ok());
    CHECK(r.mesh.origin.empty());
}

TEST_CASE("half-edge: invalid input is rejected", "[halfedge]") {
    const std::vector<dmw::Vec3> p{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}};

    SECTION("index out of range") {
        CHECK_FALSE(build_halfedge(make(p, {{0, 1, 9}})).ok());
    }
    SECTION("degenerate triangle (repeated index)") {
        CHECK_FALSE(build_halfedge(make(p, {{0, 0, 1}})).ok());
    }
    SECTION("inconsistent orientation: two faces both traverse 0->1") {
        CHECK_FALSE(build_halfedge(make(p, {{0, 1, 2}, {0, 1, 3}})).ok());
    }
    SECTION("non-manifold edge: three faces share edge 0-1") {
        CHECK_FALSE(build_halfedge(make(p, {{0, 1, 2}, {1, 0, 3}, {1, 0, 4}})).ok());
    }
    SECTION("non-manifold vertex: bowtie, two fans touching only at vertex 0") {
        // Every edge is manifold (no twins at all), but vertex 0's faces form two
        // separate fans, so a one-ring walk from vertex_halfedge[0] sees only one.
        CHECK_FALSE(build_halfedge(make(p, {{0, 1, 2}, {0, 3, 4}})).ok());
    }
}
