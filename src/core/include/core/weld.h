#pragma once

#include "core/mesh.h"

namespace dmw {

// Merges duplicate vertices (D20). Output vertices are renumbered in order of first
// appearance; triangles are never dropped, even if welding collapses them (diagnostics
// flag those, M3).
//
// epsilon == 0: exact match on coordinate values (+0.0 and -0.0 count as equal).
// epsilon  > 0: vertices are visited in order; each merges into the first earlier
//               representative within Euclidean distance <= epsilon, else becomes a new
//               representative that keeps its own position. Deterministic, but not
//               transitive: A~B and B~C does not imply A~C.
// Precondition: epsilon >= 0.
TriMesh weld_vertices(const TriMesh& soup, double epsilon = 0.0);

}  // namespace dmw
