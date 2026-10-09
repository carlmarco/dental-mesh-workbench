#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "core/bvh.h"
#include "core/mesh.h"

namespace dmw {

// Undercut analysis along a path of insertion (M10b, D90): a "virtual dental surveyor".
// For a unit direction d (the direction a restoration, appliance or tool is withdrawn), a face is UNDERCUT if it
// faces away from d (n . d < -tolerance, n the unit face normal) or if it is hidden along d: a ray from its
// centroid (lifted a hair off the surface) towards +d hits the mesh. Faces parallel to d count as reachable. Everything else is reachable along d. The border between the two on a crown is the
// survey line (height of contour). Only faces in `region` are classified (empty = all faces); every face of the
// mesh can occlude.
struct UndercutResult {
    std::vector<std::uint8_t> undercut;  // per face: 1 undercut, 0 reachable (0 outside the region)
    double undercut_area = 0.0, region_area = 0.0;
    double fraction() const { return region_area > 0.0 ? undercut_area / region_area : 0.0; }
};

// `back_facing_tolerance`: a face counts as back-facing only if n . d < -tolerance (n unit). 0 is the strict
// definition; a small positive value ignores faces that are numerically parallel to d.
// `threads`: 0 = all hardware threads (native builds; WebAssembly builds are single-threaded), 1 = serial.
UndercutResult undercut_map(const TriMesh& mesh, const Bvh& bvh, const Vec3& d, std::span<const std::uint8_t> region = {},
                            double back_facing_tolerance = 0.0, int threads = 0);

// Path of insertion minimizing the undercut area of `region`: directions within `max_tilt_degrees` of `hint`
// sampled on a Fibonacci spiral (`samples`), then a local pattern search on the sphere around the best sample
// (step halving down to `min_step_degrees`). Returns the axis and its undercut result.
// Coarse-to-fine (D91): with more than `coarse_faces` region faces (0 = off), the spiral and the main pattern search
// evaluate a regular subsample of the faces; a short pattern search at full resolution (steps of 1 degree down to
// `min_step_degrees`) and the final evaluation use every face.
struct InsertionAxis {
    Vec3 axis;
    UndercutResult result;
    std::size_t evaluations = 0;
};
InsertionAxis best_insertion_axis(const TriMesh& mesh, const Bvh& bvh, std::span<const std::uint8_t> region, const Vec3& hint,
                                  double max_tilt_degrees = 30.0, int samples = 120, double min_step_degrees = 0.25,
                                  double back_facing_tolerance = 0.0, std::size_t coarse_faces = 2000, int threads = 0);

}  // namespace dmw
