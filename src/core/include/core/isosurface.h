#pragma once

#include "core/mesh.h"
#include "core/sdf.h"

namespace dmw {

// Iso-surface {phi = iso} of a grid field by marching tetrahedra (M10e, D95). Each grid cube is split into the 6
// tetrahedra of the Kuhn (Freudenthal) subdivision along its main diagonal; the subdivision is consistent between
// neighbouring cubes, so the surface is watertight wherever it does not reach the grid boundary, with no ambiguous
// cases (unlike marching cubes). Vertices are linear interpolations on grid edges (shared by edge key); faces are
// oriented so their normals point towards increasing phi (outward for phi = signed distance).
TriMesh extract_isosurface(const Grid3& field, double iso);

}  // namespace dmw
