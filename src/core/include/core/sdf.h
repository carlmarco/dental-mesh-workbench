#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "core/mesh.h"

namespace dmw {

// Scalar field on a regular grid with cell-centred nodes: node (i, j, k) sits at origin + h * (i, j, k) and is the
// centre of a cell of size h (so the box spans origin - h/2 ... origin + h*(n - 1/2) along each axis).
struct Grid3 {
    Vec3 origin;
    double h = 1.0;
    std::array<std::size_t, 3> n{0, 0, 0};
    std::vector<double> values;  // x fastest: index = (k * n[1] + j) * n[0] + i

    std::size_t index(std::size_t i, std::size_t j, std::size_t k) const { return (k * n[1] + j) * n[0] + i; }
    Vec3 position(std::size_t i, std::size_t j, std::size_t k) const;
    double sample(const Vec3& p) const;  // trilinear, clamped to the grid
};

// Neumann-box solvers by the orthonormal DCT-II along each axis (D94): the discrete Laplacian with zero-flux
// boundaries, L = -D^T D (D = forward differences on the n - 1 edges of each line), is diagonal in that basis with
// eigenvalues -sum_axes (2 - 2 cos(pi k / n)) / h^2. Both solves are exact up to rounding (no iteration).
//   screened:  (I - t L) u = f
//   poisson:   L u = f, with the constant mode removed (f's mean is ignored; the solution has mean 0)
void dct_solve_screened(Grid3& field, double t);
void dct_solve_poisson(Grid3& field);

// Generalized signed distance by the signed heat method (Feng & Crane, ACM TOG 43(4), 2024), D94:
//   1. diffuse the area-weighted surface normals (splatted onto the grid) by one backward-Euler heat step
//      (I - t L) Y = N, t = t_coef * h^2, componentwise;
//   2. normalize X = Y / |Y|;
//   3. integrate: least squares min sum_edges (D phi - X_edge)^2, i.e. D^T D phi = D^T X (the Neumann Poisson problem);
//   4. shift phi so its area-weighted mean over the input surface is 0.
// Positive outside (normals point out). Robust to holes: the diffused normals complete the shape implicitly.
struct SignedHeatParams {
    double h = 0.0;        // grid spacing; 0 = bounding-box diagonal / 64
    double padding = 0.0;  // margin around the bounding box; 0 = 8 h
    double t_coef = 1.0;   // diffusion time t = t_coef * h^2
};
Grid3 signed_heat_distance(const TriMesh& surface, const SignedHeatParams& params = {});

}  // namespace dmw
