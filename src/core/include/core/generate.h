#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "core/mesh.h"

namespace dmw {

// Synthetic test meshes (D23, D27). All are outward / +z wound (counter-clockwise seen
// from outside). Invalid parameters return an empty mesh (errors are values, D15).

// Unit square [0,1]^2 in the z = 0 plane, nx by ny cells, two triangles per cell.
// A topological disk: chi = 1, one boundary loop. Requires nx, ny >= 1.
TriMesh make_grid(std::uint32_t nx, std::uint32_t ny);

// make_grid with the listed cells {i, j} removed (both triangles). Each removed interior
// cell that shares no vertex with another removed cell adds one hole: chi = 1 - h,
// h + 1 boundary loops. Cells sharing only a corner produce a non-manifold (bowtie) vertex.
TriMesh make_grid_with_holes(std::uint32_t nx, std::uint32_t ny,
                             const std::vector<std::array<std::uint32_t, 2>>& hole_cells);

// Torus of revolution about z: nu cells around the major circle, nv around the tube.
// Closed, genus 1, chi = 0. Requires nu, nv >= 3 (fewer would duplicate edges).
TriMesh make_torus(std::uint32_t nu, std::uint32_t nv, double major_radius = 1.0,
                   double minor_radius = 0.3);

// Icosahedron subdivided `subdivisions` times (each triangle -> 4), vertices projected to
// the sphere. V = 10*4^s + 2, E = 30*4^s, F = 20*4^s. Closed, genus 0. Requires s <= 10.
TriMesh make_icosphere(std::uint32_t subdivisions, double radius = 1.0);

// Moebius strip: one row of `segments` quads around a circle with a half twist.
// Non-orientable, chi = 0, one boundary loop. V = 2n, E = 4n, F = 2n. Requires n >= 3.
TriMesh make_mobius(std::uint32_t segments, double radius = 1.0, double half_width = 0.3);

// Disjoint union: appends src to dst, offsetting src's indices.
void append(TriMesh& dst, const TriMesh& src);

}  // namespace dmw
