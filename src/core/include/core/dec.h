#pragma once

#include <span>
#include <vector>

#include "core/halfedge.h"
#include "core/sparse.h"
#include "core/topology.h"

namespace dmw {

// Discrete exterior calculus on a triangle mesh (D45). Discrete k-forms are values on
// k-simplices: Omega^0 (vertices) --d0--> Omega^1 (edges) --d1--> Omega^2 (faces).
//
//   d0 (E x V): (d0 f)_e = f(v1) - f(v0) for edge e = (v0 -> v1)   [discrete gradient]
//   d1 (F x E): (d1 w)_f = sum of w over f's boundary edges, +1 where the edge's orientation
//               agrees with the face's winding, -1 where it opposes   [discrete curl]
//   d1 * d0 = 0 exactly (a boundary has no boundary); ker d0 = constants per component (H^0).
//
//   *0 (V): lumped vertex area = one third of the incident triangle areas (barycentric)
//   *1 (E): (cot a + cot b) / 2 with a, b the angles opposite the edge (one term on a boundary)
//   L = -d0^T *1 d0: the cotan Laplacian, symmetric negative semidefinite, "integrated":
//       (L u)_i = sum_j (cot a_ij + cot b_ij) / 2 * (u_j - u_i).
//   Caveat: an angle > 90 degrees has a negative cotangent, so *1 can be negative on obtuse
//   triangles, and L then loses the maximum principle (motivation for M6b).
struct DecOperators {
    std::vector<Edge> edges;  // oriented v0 -> v1 with v0 < v1
    SparseMatrix d0;
    SparseMatrix d1;
    std::vector<double> star0;
    std::vector<double> star1;
    SparseMatrix laplacian;
    double mean_edge_length = 0.0;
};

DecOperators build_dec(const HalfEdgeMesh& mesh);

// Per-face gradient of a piecewise-linear vertex function:
//   grad u = 1/(2A) * sum_i u_i (N x e_i), e_i = the edge opposite corner i (counter-clockwise).
// Zero-area faces get the zero vector.
std::vector<Vec3> face_gradient(const HalfEdgeMesh& mesh, std::span<const double> u);

// Integrated divergence at vertices of a per-face vector field (Crane et al. 2013, Sec. 3.2):
//   (div X)_i = 1/2 * sum over faces at i of [cot t1 (e1 . X_f) + cot t2 (e2 . X_f)],
// e1, e2 the edges leaving i, t1, t2 the opposite angles. Defined so div(grad u) = L u.
std::vector<double> vertex_divergence(const HalfEdgeMesh& mesh, std::span<const Vec3> field);

}  // namespace dmw
