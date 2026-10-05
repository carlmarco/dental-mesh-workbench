#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <vector>

#include "core/dec.h"
#include "core/halfedge.h"

namespace dmw {

// An intrinsic triangulation: the original vertex set, a (possibly different) connectivity,
// and one length per edge. Geometry is a piecewise-flat metric: each face is the Euclidean
// triangle with those side lengths. Everything intrinsic (angles, areas, the Laplacian,
// gradients, divergence, angle defects) is computed from lengths alone (D52).
struct IntrinsicTriangulation {
    HalfEdgeMesh connectivity;   // positions are the original (extrinsic) vertex positions
    std::vector<double> length;  // per half-edge; twins carry equal lengths
    std::uint32_t flips = 0;     // edge flips performed by flip_to_delaunay
    std::uint32_t skipped = 0;   // non-Delaunay edges left unflipped (D53: would need a loop or
                                 // a duplicate edge, which this connectivity cannot represent)
};

// The input mesh itself, with Euclidean edge lengths.
IntrinsicTriangulation intrinsic_from_mesh(const HalfEdgeMesh& mesh);

// Intrinsic Delaunay condition for the edge of half-edge h: cot(alpha) + cot(beta) >= 0, where
// alpha, beta are the angles opposite the edge (equivalently alpha + beta <= pi). Boundary
// edges are always Delaunay.
bool is_delaunay(const IntrinsicTriangulation& t, std::uint32_t h);

// Edge flipping (Fisher, Springborn, Bobenko, Schroeder 2007): flip non-Delaunay interior edges
// until none remain. Each flip replaces the diagonal of the two adjacent triangles by the other
// diagonal of the unfolded quadrilateral, whose length is computed by laying it out in the
// plane. The metric (the surface) is unchanged; only which intrinsic lines are edges changes.
void flip_to_delaunay(IntrinsicTriangulation& t);

// Intrinsic Delaunay triangulation of a mesh (intrinsic_from_mesh + flip_to_delaunay).
IntrinsicTriangulation intrinsic_delaunay(const HalfEdgeMesh& mesh);

// DEC operators from edge lengths only (same definitions as build_dec, D45). For the input
// mesh with Euclidean lengths this reproduces build_dec(mesh).
DecOperators build_dec(const IntrinsicTriangulation& t);

// Per-vertex mixed area (Meyer et al.) from lengths: Voronoi split on non-obtuse triangles,
// area/2 and area/4 on obtuse ones.
std::vector<double> mixed_area(const IntrinsicTriangulation& t);

// Interior angle sum at each vertex (cone angle); invariant under intrinsic flips.
std::vector<double> angle_sums(const IntrinsicTriangulation& t);

// Gradient and divergence from lengths. Each face gets its own 2D frame: corner 0 at the
// origin, corner 1 on the +x axis, corner 2 above it. Face vectors are expressed in that frame;
// divergence reads them back in the same frame, so div(grad u) = L u as in 3D.
using Vec2 = std::array<double, 2>;
std::vector<Vec2> intrinsic_gradient(const IntrinsicTriangulation& t, std::span<const double> u);
std::vector<double> intrinsic_divergence(const IntrinsicTriangulation& t, std::span<const Vec2> field);

}  // namespace dmw
