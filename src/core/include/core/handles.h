#pragma once

#include <cstdint>
#include <vector>

#include "core/halfedge.h"

namespace dmw {

// One homology generator: a closed cycle of mesh edges (vertices in order; the last connects back
// to the first), with its length in mesh units (mm for scans).
struct HandleLoop {
    std::vector<std::uint32_t> vertices;
    double length = 0.0;
};

// Handle loops (M8b, D64): 2g generator cycles per component of genus g, shortest first.
//
// Tree-cotree decomposition (Eppstein, SODA 2003) with the greedy choice of Erickson & Whittlesey
// (SODA 2005): T = Dijkstra shortest-path tree from a basepoint x per component; C = maximum
// spanning tree of the dual graph over edges not in T, weighted by the loop length
// sigma(e) = d(x,u) + l(e) + d(x,v); the edges in neither (exactly 2g) close the shortest system of
// loops through x. Each loop is returned trimmed to its cycle u -> lca(u,v) -> v -> u (dropping
// the shared stem from x, which does not change its homology class).
//
// The count (2g), validity and independence are exact. The cycles are NOT the shortest in their
// homology classes: on the plate with a handle, the loop homologous to the tube's cross-section
// measures 1.158 vs 0.333 for the tight cycle (D64). Lengths are upper bounds on handle size.
//
// Boundary loops are capped by one virtual dual node each, so holes are not reported: only the
// 2g handle generators of the capped surface are. Requires a valid half-edge mesh (manifold,
// consistently oriented).
std::vector<HandleLoop> handle_loops(const HalfEdgeMesh& mesh);

}  // namespace dmw
