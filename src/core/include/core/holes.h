#pragma once

#include <array>
#include <functional>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

#include "core/mesh.h"

namespace dmw {

// Hole filling after Liepa, "Filling Holes in Meshes" (SGP 2003), D98:
//   1. boundary loops of the (manifold, oriented) mesh;
//   2. minimum-weight triangulation of each loop by dynamic programming, O(n^3): weight = (largest dihedral angle
//      to neighbouring triangles, including the existing faces along the rim; then total area), compared
//      lexicographically, so the patch avoids creases first and is small second. As in Liepa's paper, the dihedral
//      angle across a chord is taken against the apex the sub-polygon chose for itself; that apex was optimized
//      without knowing the triangle on the other side, so optimal substructure fails and the result is not exactly
//      optimal for the dihedral term (an exact DP needs the outer apex in its state, O(n^4)). With area alone the
//      DP is exact (minimum-area triangulation, Barequet & Sharir 1995). Measured gap: tests/test_holes.cpp;
//   3. refinement: split patch triangles at their centroid while the centroid is farther than sqrt(2) times the
//      local edge-length scale from each corner, relaxing patch edges by Delaunay flips after each round, so the
//      patch matches the density of its surroundings;
//   4. fairing: the patch interior minimizes the discrete thin-plate energy |L x|^2 with the rim fixed, solved as
//      (K D^-1 K)[interior] x = -(K D^-1 K)[interior, fixed] x_fixed (K = graph Laplacian, D = degrees), SPD.
// Loops longer than `max_loop` (e.g. the open cut of an intraoral scan) are left open. Patches never create an edge
// the mesh already has (neither by a triangulation chord nor by a refinement flip), so the result stays manifold.

// Directed boundary loops: each loop lists vertices so that consecutive pairs are boundary edges a -> b of the
// existing faces (the hole lies to the right, the faces to the left).
std::vector<std::vector<std::uint32_t>> boundary_loops(const TriMesh& mesh);

// Minimum-weight triangulation of one loop (positions of the loop vertices in order). `outside[i]` is the third
// vertex of the existing face on loop edge (i, i+1), used for the rim dihedral angles. Returns triangles as index
// triples into the loop, oriented consistently with the surrounding faces.
// `chord_ok(i, j)` (loop indices, i < j) may forbid chords, e.g. pairs of loop vertices that are already joined by an
// edge of the mesh outside the hole (using such a chord would put three faces on one edge). Returns an empty list if
// no triangulation avoids the forbidden chords.
std::vector<std::array<std::uint32_t, 3>> triangulate_loop(std::span<const Vec3> loop, std::span<const Vec3> outside,
                                                           bool use_dihedral = true,
                                                           const std::function<bool(std::uint32_t, std::uint32_t)>& chord_ok = {});

// The weight used by triangulate_loop, exposed for tests: (largest dihedral angle in radians, total area) of a set of
// loop triangles, including the rim faces.
std::pair<double, double> triangulation_weight(std::span<const Vec3> loop, std::span<const Vec3> outside,
                                               std::span<const std::array<std::uint32_t, 3>> triangles);

struct HoleFillParams {
    std::size_t max_loop = 500;  // longer loops are skipped
    bool dihedral = true;        // weight (dihedral, area) as in Liepa; false = area only
    bool refine = true;
    bool fair = true;
};
struct HoleFillResult {
    TriMesh mesh;                       // input faces first, then the patches; new vertices appended
    std::size_t loops = 0, filled = 0;  // boundary loops found, loops filled
    std::vector<std::uint32_t> patch_vertices;  // vertices created by refinement (moved by fairing)
};
HoleFillResult fill_holes(const TriMesh& mesh, const HoleFillParams& params = {});

}  // namespace dmw
