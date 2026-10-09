#pragma once

#include <array>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

#include "core/mesh.h"

namespace dmw {

// Ray casting against a triangle mesh (M10a, D89): a bounding-volume hierarchy built top-down with the binned
// surface-area heuristic (SAH; MacDonald & Booth 1990, binning as in Wald 2007), traversed near-child-first
// with an explicit stack; ray-triangle tests by Moller & Trumbore (1997). Infrastructure for the undercut
// (10b) and wall-thickness (10c) tools.
struct RayHit {
    bool hit = false;
    double t = std::numeric_limits<double>::infinity();  // hit point = origin + t * dir
    std::uint32_t face = 0xFFFFFFFFu;
    double u = 0.0, v = 0.0;  // barycentric coordinates of the hit (weights of vertices 1 and 2)
};

struct ClosestPoint {
    Vec3 point;
    double distance = std::numeric_limits<double>::infinity();
    std::uint32_t face = 0xFFFFFFFFu;
};

class Bvh {
public:
    // Copies what it needs; the mesh may change afterwards. Leaves hold at most `leaf_size` triangles.
    Bvh(std::span<const Vec3> positions, std::span<const std::array<std::uint32_t, 3>> triangles, std::uint32_t leaf_size = 4);

    // Nearest hit with t in (t_min, t_max). `ignore` skips one face (a ray leaving that face's surface).
    RayHit intersect(const Vec3& origin, const Vec3& dir, double t_min = 0.0,
                     double t_max = std::numeric_limits<double>::infinity(), std::uint32_t ignore = 0xFFFFFFFFu) const;
    // Any hit with t in (t_min, t_max): cheaper (stops at the first hit).
    bool occluded(const Vec3& origin, const Vec3& dir, double t_min = 0.0, double t_max = std::numeric_limits<double>::infinity(),
                  std::uint32_t ignore = 0xFFFFFFFFu) const;

    // Closest point on the mesh to p (D99): nodes are visited nearest-box-first and pruned by their distance lower
    // bound; per triangle, Ericson's region-based closest point (Real-Time Collision Detection, 5.1.5).
    ClosestPoint closest(const Vec3& p) const;

    std::size_t node_count() const { return nodes_.size(); }

private:
    struct Node {
        Vec3 lo, hi;          // bounds
        std::uint32_t first;  // leaf: first index into order_; inner: index of the left child (right = left + 1)
        std::uint32_t count;  // leaf: number of triangles; inner: 0
    };
    template <bool kAny>
    RayHit trace(const Vec3& origin, const Vec3& dir, double t_min, double t_max, std::uint32_t ignore) const;

    std::vector<Vec3> v0_, e1_, e2_;      // per triangle: vertex 0 and two edge vectors (Moller-Trumbore)
    std::vector<std::uint32_t> order_;    // triangle ids, grouped by leaf
    std::vector<Node> nodes_;
};

// Reference: the same query by testing every triangle (for tests and small meshes).
// Closest point on a triangle (a, b, c) to p (Ericson): exact, handles every vertex / edge / face region.
Vec3 closest_point_on_triangle(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c);
ClosestPoint closest_brute_force(std::span<const Vec3> positions, std::span<const std::array<std::uint32_t, 3>> triangles, const Vec3& p);

RayHit intersect_brute_force(std::span<const Vec3> positions, std::span<const std::array<std::uint32_t, 3>> triangles,
                             const Vec3& origin, const Vec3& dir, double t_min = 0.0,
                             double t_max = std::numeric_limits<double>::infinity(), std::uint32_t ignore = 0xFFFFFFFFu);

}  // namespace dmw
