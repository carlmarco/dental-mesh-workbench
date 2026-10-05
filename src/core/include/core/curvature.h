#pragma once

#include <cstdint>
#include <vector>

#include "core/halfedge.h"
#include "core/mesh.h"

namespace dmw {

// Discrete curvature operators of
//   M. Meyer, M. Desbrun, P. Schroeder, A. H. Barr, "Discrete Differential-Geometry
//   Operators for Triangulated 2-Manifolds", Visualization and Mathematics III, 2003.
//
// Per vertex i with one-ring neighbors j:
//   A_i      mixed area (Voronoi area, falling back to area/2 or area/4 on obtuse triangles)
//   K(x_i) = 1/(2 A_i) * sum_j (cot a_ij + cot b_ij) (x_i - x_j)        [= 2 H n]
//   K_G    = (2 pi - sum of incident angles) / A_i
// Pointwise values are averages over A_i ("integrated quantity / area", D35).
struct CurvatureField {
    std::vector<std::uint8_t> boundary;  // 1 if the vertex lies on a boundary loop

    std::vector<double> mixed_area;    // A_mixed; sums to the total surface area
    std::vector<double> angle_defect;  // integrated K: 2pi - sum(theta) interior,
                                       // pi - sum(theta) boundary (geodesic-curvature term)
                                       // so that sum = 2 pi chi (discrete Gauss-Bonnet)

    // Interior vertices only; NaN (or the zero vector) on boundary and isolated vertices (D36).
    std::vector<double> gaussian;      // K = angle_defect / A_mixed
    std::vector<Vec3> mean_normal;     // Meyer's K(x_i) = 2 H n (mean-curvature normal)
    std::vector<double> mean;          // H = (1/2) K(x_i) . n_hat, signed; n_hat from
                                       // area-weighted face normals (D37). Sphere: +1/r.
    std::vector<double> k1, k2;        // principal: H +- sqrt(max(H^2 - K, 0)), k1 >= k2

    std::vector<std::uint32_t> degenerate_faces;  // zero area; excluded from cotan sums (D38)
    std::uint32_t intrinsic_flips = 0;            // edge flips used (intrinsic_delaunay only)
};

// Requires a valid half-edge mesh: the operators assume a 2-manifold (with boundary).
//
// intrinsic_delaunay = true (M6b, D54): A_i and K(x_i) come from the intrinsic Delaunay
// triangulation instead: K(x_i) = -(L_idt x)_i / A_i with A_i its mixed area. Angle defects
// (hence K_G) are unchanged by intrinsic flips: cone angles are intrinsic. Reduces the
// blow-ups of cotan mean curvature on obtuse triangles; does not restore pointwise convergence
// (measured, D55).
CurvatureField compute_curvature(const HalfEdgeMesh& mesh, bool intrinsic_delaunay = false);

}  // namespace dmw
