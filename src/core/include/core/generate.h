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

// A disk with one handle: an n x n unit grid with two square holes (cells (n/4, n/2) and
// (3n/4 - 1, n/2)) whose rims are joined by an arched square tube of `tube_segments` segments.
// chi = -1, one boundary loop (the outer rim), genus 1. Requires n >= 8, tube_segments >= 2.
TriMesh make_plate_with_handle(std::uint32_t n = 12, std::uint32_t tube_segments = 16,
                               double arch_height = 0.35);

// Torus of revolution about z: nu cells around the major circle, nv around the tube.
// Closed, genus 1, chi = 0. Requires nu, nv >= 3 (fewer would duplicate edges).
TriMesh make_torus(std::uint32_t nu, std::uint32_t nv, double major_radius = 1.0,
                   double minor_radius = 0.3);

// Open cylinder of radius r about the z axis, z in [0, height]: nu segments around, nv rings
// of cells. chi = 0, two boundary loops. Interior vertices are exactly flat (K = 0) while
// H -> 1/(2r): the test case that separates mean from Gaussian curvature. nu >= 3, nv >= 1.
TriMesh make_cylinder(std::uint32_t nu, std::uint32_t nv, double radius = 1.0, double height = 1.0);

// Icosahedron subdivided `subdivisions` times (each triangle -> 4), vertices projected to
// the sphere. V = 10*4^s + 2, E = 30*4^s, F = 20*4^s. Closed, genus 0. Requires s <= 10.
TriMesh make_icosphere(std::uint32_t subdivisions, double radius = 1.0);

// Moebius strip: one row of `segments` quads around a circle with a half twist.
// Non-orientable, chi = 0, one boundary loop. V = 2n, E = 4n, F = 2n. Requires n >= 3.
TriMesh make_mobius(std::uint32_t segments, double radius = 1.0, double half_width = 0.3);

// One mesh containing every defect the topology report detects, for demos and tests:
// a 4x4 grid (vertex (i,j) = 5j + i) with
//   - triangle 12,13,18 flipped          -> 3 misoriented edges
//   - a fin on interior edge 6-7         -> 1 non-manifold edge
//   - a triangle touching corner 24 only -> bowtie (non-manifold) vertex 24
//   - a copy of face 24                  -> duplicate face
//   - face (0, 0, 1)                     -> invalid face
//   - one unreferenced position          -> isolated vertex
// plus a Moebius strip (second component, non-orientable) beside it.
TriMesh make_defect_showcase();

// Flat n x n unit grid with alternate interior rows shifted by +-shift cells: "brick"
// triangles with obtuse angles. At shift 0.45 the cotan weights go down to -0.9 and heat
// diffusion breaks the maximum principle; the intrinsic Delaunay Laplacian fixes it (M6b, D55).
TriMesh make_brick_grid(std::uint32_t n, double shift = 0.45);

// Disjoint union: appends src to dst, offsetting src's indices.
void append(TriMesh& dst, const TriMesh& src);

}  // namespace dmw
