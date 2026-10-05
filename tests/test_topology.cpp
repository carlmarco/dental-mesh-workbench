#include <array>
#include <cstdint>
#include <optional>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "core/generate.h"
#include "core/halfedge.h"
#include "core/topology.h"

using namespace dmw;
using Ids = std::vector<std::uint32_t>;
using Betti = std::array<std::uint32_t, 3>;

namespace {

TriMesh tetra() {
    return {{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}}};
}

TriMesh cube() {
    return {{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}},
            {{0, 3, 2}, {0, 2, 1}, {4, 5, 6}, {4, 6, 7}, {0, 1, 5}, {0, 5, 4},
             {2, 3, 7}, {2, 7, 6}, {0, 4, 7}, {0, 7, 3}, {1, 2, 6}, {1, 6, 5}}};
}

// Every component with a genus must satisfy b0 - b1 + b2 = chi (Euler-Poincare).
void check_euler_poincare(const TopologyReport& r) {
    for (const auto& c : r.components) {
        if (!c.betti) continue;
        const auto& b = *c.betti;
        CHECK(std::int64_t{b[0]} - std::int64_t{b[1]} + std::int64_t{b[2]} == c.euler_characteristic);
    }
}

}  // namespace

TEST_CASE("topology: tetrahedron is a closed genus-0 surface", "[topology]") {
    const auto r = analyze_topology(tetra());
    REQUIRE(r.components.size() == 1);
    const auto& c = r.components[0];
    CHECK(c.num_vertices == 4);
    CHECK(c.num_edges == 6);
    CHECK(c.num_faces == 4);
    CHECK(c.euler_characteristic == 2);
    CHECK(c.boundary_loops == 0);
    CHECK(c.manifold);
    CHECK(c.orientable);
    CHECK(c.consistently_oriented);
    CHECK(c.genus == std::optional<std::uint32_t>{0});
    CHECK(c.betti == std::optional<Betti>{{1, 0, 1}});
    CHECK(r.count(EdgeKind::Manifold) == 6);
    check_euler_poincare(r);
}

TEST_CASE("topology: icosphere is genus 0 at every resolution", "[topology]") {
    for (std::uint32_t s = 0; s <= 3; ++s) {
        const auto r = analyze_topology(make_icosphere(s));
        REQUIRE(r.components.size() == 1);
        CHECK(r.components[0].euler_characteristic == 2);
        CHECK(r.components[0].genus == std::optional<std::uint32_t>{0});
    }
}

TEST_CASE("topology: torus has genus 1 and b1 = 2", "[topology]") {
    const auto r = analyze_topology(make_torus(12, 6));
    REQUIRE(r.components.size() == 1);
    const auto& c = r.components[0];
    CHECK(c.euler_characteristic == 0);
    CHECK(c.boundary_loops == 0);
    CHECK(c.genus == std::optional<std::uint32_t>{1});
    CHECK(c.betti == std::optional<Betti>{{1, 2, 1}});
    check_euler_poincare(r);
}

TEST_CASE("topology: grid is a disk with one boundary loop", "[topology]") {
    const auto r = analyze_topology(make_grid(4, 3));
    REQUIRE(r.components.size() == 1);
    const auto& c = r.components[0];
    CHECK(c.euler_characteristic == 1);
    CHECK(c.boundary_loops == 1);
    CHECK(c.genus == std::optional<std::uint32_t>{0});
    CHECK(c.betti == std::optional<Betti>{{1, 0, 0}});
    CHECK(r.count(EdgeKind::Boundary) == 2 * (4 + 3));  // the perimeter
    check_euler_poincare(r);
}

TEST_CASE("topology: grid with 2 holes has 3 boundary loops, b1 = 2", "[topology]") {
    const auto r = analyze_topology(make_grid_with_holes(6, 6, {{1, 1}, {4, 3}}));
    REQUIRE(r.components.size() == 1);
    const auto& c = r.components[0];
    CHECK(c.euler_characteristic == -1);
    CHECK(c.boundary_loops == 3);
    CHECK(c.genus == std::optional<std::uint32_t>{0});
    CHECK(c.betti == std::optional<Betti>{{1, 2, 0}});
    check_euler_poincare(r);
}

TEST_CASE("topology: plate with a handle - same chi as two holes, but genus 1", "[topology]") {
    // chi = -1 like the two-hole grid; b = 1 instead of 3, so g = (2 + 1 - 1) / 2 = 1.
    const auto r = analyze_topology(make_plate_with_handle());
    REQUIRE(r.components.size() == 1);
    const auto& c = r.components[0];
    CHECK(c.euler_characteristic == -1);
    CHECK(c.boundary_loops == 1);
    CHECK(c.manifold);
    CHECK(c.consistently_oriented);  // the mirrored ring order at hole B is right
    CHECK(c.genus == std::optional<std::uint32_t>{1});
    CHECK(c.betti == std::optional<Betti>{{1, 2, 0}});
    check_euler_poincare(r);
}

TEST_CASE("topology: disjoint union reports each component separately", "[topology]") {
    auto m = tetra();  // vertices 0-3
    append(m, make_torus(8, 5));
    const auto r = analyze_topology(m);
    REQUIRE(r.components.size() == 2);
    CHECK(r.components[0].genus == std::optional<std::uint32_t>{0});
    CHECK(r.components[1].genus == std::optional<std::uint32_t>{1});
    CHECK(r.vertex_component[3] == 0);
    CHECK(r.vertex_component[4] == 1);
    CHECK(r.face_component[0] == 0);
    CHECK(r.face_component[4] == 1);
}

TEST_CASE("topology: Moebius strip is manifold but non-orientable; no genus", "[topology]") {
    const auto r = analyze_topology(make_mobius(8));
    REQUIRE(r.components.size() == 1);
    const auto& c = r.components[0];
    CHECK(c.manifold);
    CHECK_FALSE(c.orientable);
    CHECK_FALSE(c.consistently_oriented);  // some edge must be misoriented
    CHECK(c.euler_characteristic == 0);
    CHECK(c.boundary_loops == 1);  // the famous single edge
    CHECK_FALSE(c.genus.has_value());
    CHECK_FALSE(c.betti.has_value());
}

TEST_CASE("topology: one flipped face is fixable winding, not a topology change", "[topology]") {
    auto m = cube();
    m.triangles[0] = {0, 2, 3};  // was {0, 3, 2}
    const auto r = analyze_topology(m);
    REQUIRE(r.components.size() == 1);
    const auto& c = r.components[0];
    CHECK(c.orientable);                   // flipping that face back fixes everything
    CHECK_FALSE(c.consistently_oriented);
    CHECK(r.count(EdgeKind::Misoriented) == 3);  // the flipped face's three edges
    CHECK(c.genus == std::optional<std::uint32_t>{0});  // winding doesn't change topology
}

TEST_CASE("topology: fin (3 faces on one edge) is non-manifold", "[topology]") {
    const TriMesh m{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}},
                    {{0, 1, 2}, {1, 0, 3}, {1, 0, 4}}};
    const auto r = analyze_topology(m);
    REQUIRE(r.components.size() == 1);
    CHECK(r.count(EdgeKind::NonManifold) == 1);
    CHECK_FALSE(r.components[0].manifold);
    CHECK_FALSE(r.components[0].genus.has_value());
}

TEST_CASE("topology: bowtie is ONE component with a non-manifold vertex", "[topology]") {
    // Two triangles touching only at vertex 0: path-connected, so b0 = 1 (D28).
    const TriMesh m{{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {-1, 0, 0}, {-1, -1, 0}},
                    {{0, 1, 2}, {0, 3, 4}}};
    const auto r = analyze_topology(m);
    CHECK(r.components.size() == 1);
    CHECK(r.nonmanifold_vertices == Ids{0});
    CHECK_FALSE(r.components[0].manifold);
    CHECK_FALSE(r.components[0].genus.has_value());
}

TEST_CASE("topology: holes sharing a corner create a bowtie vertex", "[topology]") {
    // Cells (1,1) and (2,2) touch at grid vertex (2,2) = index 2*5+2 = 12.
    const auto r = analyze_topology(make_grid_with_holes(4, 4, {{1, 1}, {2, 2}}));
    CHECK(r.nonmanifold_vertices == Ids{12});
}

TEST_CASE("topology: invalid, duplicate and isolated elements are reported and excluded", "[topology]") {
    auto m = tetra();
    m.positions.push_back({9, 9, 9});      // vertex 4: isolated
    m.triangles.push_back({0, 0, 1});      // face 4: repeated index
    m.triangles.push_back({0, 1, 7});      // face 5: out of range
    m.triangles.push_back({2, 1, 0});      // face 6: same vertex set as face 0 (reversed)
    const auto r = analyze_topology(m);
    CHECK(r.invalid_faces == Ids{4, 5});
    CHECK(r.duplicate_faces == Ids{6});
    CHECK(r.isolated_vertices == Ids{4});
    CHECK(r.vertex_component[4] == kInvalid);
    CHECK(r.face_component[6] == kInvalid);
    REQUIRE(r.components.size() == 1);
    CHECK(r.components[0].euler_characteristic == 2);  // still the clean tetrahedron
    CHECK(r.count(EdgeKind::Manifold) == 6);
}

TEST_CASE("topology: defect showcase - every defect detected and located", "[topology]") {
    const auto m = make_defect_showcase();
    const auto r = analyze_topology(m);
    const std::uint32_t apex = 25, isolated = 28, first_strip_vertex = 29;
    CHECK(r.invalid_faces == Ids{35});
    CHECK(r.duplicate_faces == Ids{34});
    CHECK(r.isolated_vertices == Ids{isolated});
    CHECK(r.nonmanifold_vertices == Ids{24});
    CHECK(r.count(EdgeKind::NonManifold) == 1);
    for (std::size_t e = 0; e < r.edges.size(); ++e) {
        if (r.edge_kind[e] == EdgeKind::NonManifold) {
            CHECK(r.edges[e].v0 == 6);
            CHECK(r.edges[e].v1 == 7);
        }
    }
    REQUIRE(r.components.size() == 2);
    CHECK(r.vertex_component[apex] == 0);
    CHECK_FALSE(r.components[0].manifold);
    CHECK_FALSE(r.components[0].consistently_oriented);
    CHECK(r.vertex_component[first_strip_vertex] == 1);
    CHECK(r.components[1].manifold);
    CHECK_FALSE(r.components[1].orientable);
}

TEST_CASE("topology: empty mesh has no components", "[topology]") {
    const auto r = analyze_topology(TriMesh{});
    CHECK(r.components.empty());
    CHECK(r.edges.empty());
}

TEST_CASE("topology: excluding collapsed and duplicate faces lets a scan-like mesh build", "[topology]") {
    // Real scans contain triangles like (v, v, v). The rest is a clean tetrahedron.
    auto m = tetra();
    m.triangles.push_back({2, 2, 2});  // collapsed
    m.triangles.push_back({3, 1, 0});  // duplicate of face 1 = {0, 1, 3}
    REQUIRE_FALSE(build_halfedge(m).ok());
    const auto cleaned = without_excluded_faces(m, analyze_topology(m));
    CHECK(cleaned.positions.size() == m.positions.size());  // vertex indexing preserved
    CHECK(cleaned.triangles.size() == 4);
    CHECK(build_halfedge(cleaned).ok());
}

TEST_CASE("topology: excluding faces does not hide non-manifold defects", "[topology]") {
    const auto m = make_defect_showcase();
    CHECK_FALSE(build_halfedge(without_excluded_faces(m, analyze_topology(m))).ok());  // fin, bowtie remain
}

TEST_CASE("topology: manifold analysis view of the defect showcase builds, excluding few faces", "[topology]") {
    const auto m = make_defect_showcase();
    const auto a = manifold_analysis_mesh(m);
    REQUIRE(a.manifold);
    CHECK(build_halfedge(a.mesh).ok());
    CHECK(a.mesh.positions.size() == m.positions.size());  // vertex indexing unchanged
    CHECK(a.excluded_faces > 0);
    CHECK(a.excluded_faces < 20);  // fin edge, bowtie, flipped face and neighbours, seam, invalid, duplicate
}

TEST_CASE("topology: a clean mesh passes through untouched", "[topology]") {
    const auto a = manifold_analysis_mesh(make_torus(8, 5));
    CHECK(a.manifold);
    CHECK(a.excluded_faces == 0);
    CHECK(a.passes == 0);
}
