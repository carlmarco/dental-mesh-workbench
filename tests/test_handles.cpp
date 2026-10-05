#include <algorithm>
#include <cmath>
#include <cstdint>
#include <set>
#include <utility>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/generate.h"
#include "core/handles.h"
#include "core/topology.h"

using namespace dmw;
using Catch::Matchers::WithinAbs;

namespace {

HalfEdgeMesh he(const TriMesh& m) {
    auto r = build_halfedge(m);
    REQUIRE(r.ok());
    return std::move(r.mesh);
}

using EdgeSet = std::set<std::pair<std::uint32_t, std::uint32_t>>;

EdgeSet loop_edges(const HandleLoop& l) {
    EdgeSet s;
    for (std::size_t i = 0; i < l.vertices.size(); ++i) {
        const std::uint32_t a = l.vertices[i], b = l.vertices[(i + 1) % l.vertices.size()];
        s.insert({std::min(a, b), std::max(a, b)});
    }
    return s;
}

// A loop is a closed walk along mesh edges: consecutive vertices are adjacent, no repeats.
void check_cycle(const HalfEdgeMesh& m, const HandleLoop& l) {
    REQUIRE(l.vertices.size() >= 3);
    const std::set<std::uint32_t> distinct(l.vertices.begin(), l.vertices.end());
    CHECK(distinct.size() == l.vertices.size());
    double length = 0.0;
    for (std::size_t i = 0; i < l.vertices.size(); ++i) {
        const std::uint32_t a = l.vertices[i], b = l.vertices[(i + 1) % l.vertices.size()];
        const auto ring = one_ring(m, a);
        CHECK(std::find(ring.begin(), ring.end(), b) != ring.end());
        const Vec3 &p = m.positions[a], &q = m.positions[b];
        length += std::sqrt((p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y) + (p.z - q.z) * (p.z - q.z));
    }
    CHECK_THAT(l.length, WithinAbs(length, 1e-9));
}

// Cutting the surface along the given edges leaves its faces connected iff the cut is
// non-separating. A loop that bounds a disk (homologically trivial) would separate.
bool faces_connected_after_cut(const HalfEdgeMesh& m, const EdgeSet& cut) {
    const std::size_t nf = m.origin.size() / 3;
    std::vector<bool> seen(nf, false);
    std::vector<std::uint32_t> queue{0};
    seen[0] = true;
    for (std::size_t q = 0; q < queue.size(); ++q) {
        for (std::uint32_t k = 0; k < 3; ++k) {
            const std::uint32_t h = 3 * queue[q] + k, t = m.twin[h];
            if (t == kInvalid) continue;
            const std::uint32_t a = m.origin[h], b = dest(m, h);
            if (cut.count({std::min(a, b), std::max(a, b)}) != 0) continue;
            if (!seen[face(t)]) seen[face(t)] = true, queue.push_back(face(t));
        }
    }
    return queue.size() == nf;
}

}  // namespace

TEST_CASE("handles: genus-0 surfaces have none, even with holes (boundaries are capped)", "[handles]") {
    CHECK(handle_loops(he(make_icosphere(2))).empty());
    CHECK(handle_loops(he(make_grid_with_holes(8, 8, {{2, 2}, {5, 4}}))).empty());  // b = 3
    CHECK(handle_loops(he(make_cylinder(12, 4))).empty());                         // b = 2
}

TEST_CASE("handles: torus has 2 independent non-separating loops, shortest = tube circle", "[handles]") {
    const auto mesh = he(make_torus(16, 8, 1.0, 0.3));
    const auto loops = handle_loops(mesh);
    REQUIRE(loops.size() == 2);  // 2g, g = 1
    for (const auto& l : loops) {
        check_cycle(mesh, l);
        CHECK(faces_connected_after_cut(mesh, loop_edges(l)));  // non-trivial in homology
    }
    EdgeSet both = loop_edges(loops[0]);
    const EdgeSet second = loop_edges(loops[1]);
    both.insert(second.begin(), second.end());
    CHECK(faces_connected_after_cut(mesh, both));  // independent: cutting both still connected
    // Shortest non-trivial cycle: the 8-gon around the tube, 8 * 2 r sin(pi / 8).
    CHECK_THAT(loops[0].length, WithinAbs(8 * 2 * 0.3 * std::sin(std::acos(-1.0) / 8), 1e-9));
    CHECK(loops[0].length <= loops[1].length);
}

TEST_CASE("handles: plate with a handle reports its 2 handle loops, not its boundary", "[handles]") {
    const std::uint32_t n = 12;
    const auto mesh = he(make_plate_with_handle(n));
    const auto loops = handle_loops(mesh);
    REQUIRE(loops.size() == 2);  // genus 1; the outer rim (b = 1) is not a generator
    for (const auto& l : loops) {
        check_cycle(mesh, l);
        CHECK(faces_connected_after_cut(mesh, loop_edges(l)));
    }
    // Lower bound only: no non-trivial cycle is shorter than the tube's square cross-section,
    // perimeter 4h with h = 1/n. The loops are NOT shortest representatives of their classes:
    // measured, the loop homologous to the tube's cross-section circles hole A on the plate at
    // length 1.158, vs 0.333 for the tight cycle (D64). Trimming the stem gives valid generators,
    // not tight ones.
    const double tube = 4.0 / n;
    CHECK(loops[0].length >= tube - 1e-9);
}

TEST_CASE("handles: loop count is 2g per component, summed over components", "[handles]") {
    auto m = make_torus(10, 6);
    append(m, make_icosphere(1));
    append(m, make_plate_with_handle(10));
    const auto loops = handle_loops(he(m));
    std::uint32_t total_genus = 0;
    for (const auto& c : analyze_topology(m).components) total_genus += c.genus.value_or(0);
    CHECK(total_genus == 2);
    CHECK(loops.size() == 2 * total_genus);
}
