#pragma once

#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "core/cholesky.h"
#include "core/dec.h"
#include "core/halfedge.h"
#include "core/intrinsic.h"

namespace dmw {

struct GeodesicResult {
    std::vector<double> distance;  // per vertex; NaN where no source is reachable
    double time_step = 0.0;
    std::string error;  // empty on success
    bool ok() const { return error.empty(); }
};

// Heat method of
//   K. Crane, C. Weischedel, M. Wardetzky, "Geodesics in Heat: A New Approach to Computing
//   Distance Based on Heat Flow", ACM Transactions on Graphics 32(5), 2013.
//   I.   (*0 - t L) u = delta_sources            (one backward-Euler heat step, t = m h^2)
//   II.  X = -grad u / |grad u|                  (per face; unit, pointing away from sources)
//   III. L phi = div X, shifted so phi = 0 at the sources.
// Neumann (natural) boundary conditions (D48). Both systems are factored once in the
// constructor (envelope LDL^T, D49); each query is then two pairs of triangular solves.
// Iterative CG is NOT used for step I: it cannot resolve the exponentially small heat values
// far from the source (measured, D47).
//
// With intrinsic_delaunay = true the operators come from the intrinsic Delaunay triangulation
// of the mesh (Bobenko & Springborn 2007; flips per Fisher et al. 2007): same vertices and
// surface, all cotan weights >= 0 (M6b). Distances are still per original vertex.
class HeatGeodesics {
public:
    explicit HeatGeodesics(const HalfEdgeMesh& mesh, double time_factor = 1.0, bool intrinsic_delaunay = false);
    bool ok() const { return error_.empty(); }
    const std::string& error() const { return error_; }
    GeodesicResult distance(std::span<const std::uint32_t> sources) const;
    const DecOperators& operators() const { return ops_; }
    const IntrinsicTriangulation& triangulation() const { return tri_; }
    std::size_t factor_entries() const { return heat_.stored_entries() + poisson_.stored_entries(); }

private:
    IntrinsicTriangulation tri_;  // owned: no pointer back to the caller's mesh
    DecOperators ops_;
    double time_step_;
    EnvelopeLdlt heat_;                  // *0 - t L, SPD
    std::vector<std::uint32_t> pinned_;  // one vertex per component (declared before poisson_:
                                         // members initialize in declaration order)
    EnvelopeLdlt poisson_;               // -L with the pinned rows/columns replaced by identity, SPD
    std::string error_;
};

// Baseline: shortest paths along mesh edges (Dijkstra, Euclidean edge lengths). Never below
// the true geodesic distance; NaN where unreachable.
std::vector<double> dijkstra_distance(const HalfEdgeMesh& mesh, std::span<const std::uint32_t> sources);

}  // namespace dmw
