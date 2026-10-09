#pragma once

#include <cstdint>
#include <limits>
#include <vector>

#include "core/bvh.h"
#include "core/mesh.h"

namespace dmw {

// Wall thickness of a closed solid (M10c, D92): for each vertex, rays are cast INTO the solid (around the inward
// vertex normal) and only rays that leave through the opposite wall count (the hit face's outward normal points
// along the ray). Three statistics per vertex:
//   along_normal  the single ray along the inward normal (exact between parallel walls)
//   cone_min      the shortest accepted ray in the cone: the thinnest local wall (what a printing/milling check wants)
//   cone_median   the median of accepted rays: robust, after the shape diameter function (Shapira, Shamir,
//                 Cohen-Or 2008; they use an outlier-trimmed weighted mean, this uses the median)
// NaN where no ray was accepted (open surfaces, outward-facing normals). Outward orientation is assumed.
struct ThicknessParams {
    int rays = 30;               // rays in the cone (plus the normal ray)
    double cone_degrees = 30.0;  // half-angle of the cone around the inward normal
    double max_distance = std::numeric_limits<double>::infinity();
    int threads = 0;             // 0 = all hardware threads (native); WebAssembly builds are single-threaded
};
struct ThicknessField {
    std::vector<double> along_normal, cone_min, cone_median;
};
ThicknessField wall_thickness(const TriMesh& mesh, const Bvh& bvh, const ThicknessParams& params = {});

}  // namespace dmw
