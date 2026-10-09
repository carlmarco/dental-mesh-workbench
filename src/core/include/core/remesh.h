#pragma once

#include "core/mesh.h"

namespace dmw {

// Isotropic remeshing after Botsch & Kobbelt, "A Remeshing Approach to Multiresolution Modeling" (SGP 2004), D99.
// With target edge length L, each iteration: split edges longer than 4/3 L; collapse edges shorter than 4/5 L when
// safe (link condition, no new edge longer than 4/3 L, no face flips, no valence below 3, boundary vertices never
// pulled inward); flip edges that bring the four vertices closer to valence 6 (4 on the boundary); move interior
// vertices tangentially towards their neighbours' centroid and project them back onto the original surface (closest
// point, BVH). Boundary vertices stay where they are. Requires a manifold, consistently oriented input.
struct RemeshParams {
    double target_edge = 0.0;  // L; 0 = the mean edge length of the input
    int iterations = 5;
    double relaxation = 0.5;   // fraction of the tangential step towards the centroid
};
TriMesh remesh_isotropic(const TriMesh& mesh, const RemeshParams& params = {});

}  // namespace dmw
