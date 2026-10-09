#pragma once

#include "core/mesh.h"
#include "core/sdf.h"

namespace dmw {

// Offset shells (M10e, D95): a closed solid of wall thickness `wall` behind a surface (on the side opposite its
// normals). Two constructions, kept side by side because the comparison is the point:
//   sdf    outer wall = {phi = 0} and inner wall = {phi = -wall} of the generalized signed distance (signed heat
//          method), extracted by marching tetrahedra; both iso-surfaces are closed, so no stitching. An open input is
//          completed implicitly. Grid spacing: `grid_h`, or 0 = the bounding-box diagonal / `max_nodes` per axis.
//   naive  the surface plus a copy moved by -wall along the area-weighted vertex normals (orientation flipped) and a
//          side wall stitched along the boundary. Folds and thins in concave regions (D92).
TriMesh offset_shell_sdf(const TriMesh& surface, double wall, double grid_h = 0.0, std::size_t max_nodes = 96);
TriMesh offset_shell_naive(const TriMesh& surface, double wall);

}  // namespace dmw
