#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "core/dec.h"
#include "core/halfedge.h"

namespace dmw {

// Cusp-tip detection on intraoral scans (M9, D67). Three detectors, each a scored set of vertices
// followed by non-maximum suppression (NMS):
//   RawCurvature       local maxima of Meyer mean curvature H (baseline: noise-sensitive)
//   SmoothedCurvature  local maxima of H diffused to scale sigma: (*0 - tL) h = *0 H, t = sigma^2 / 2
//   OcclusalProminence local maxima of p = h - diffuse(h, R), h = height along the occlusal axis,
//                      gated by convexity (smoothed H > 0) and by global height (h above a quantile)
// Only the largest connected component (the arch) is searched; debris is ignored.
struct CuspParams {
    enum class Method { RawCurvature, SmoothedCurvature, OcclusalProminence };
    Method method = Method::OcclusalProminence;
    double curvature_scale = 0.5;   // sigma (mm): smoothing of H (B, and C's convexity gate)
    double prominence_scale = 2.0;  // R (mm): scale of the local height baseline (C)
    double height_quantile = 0.5;   // C: candidates must lie above this quantile of arch height
    double nms_radius = 1.5;        // mm: at most one detection per ball of this radius
    double threshold = 0.0;         // minimum score (1/mm for A and B, mm for C)
};

struct CuspDetection {
    std::vector<std::uint32_t> vertices;  // detected cusp tips, strongest first
    std::vector<double> scores;
    Vec3 occlusal_axis{};                 // unit, pointing from the gingival cut toward the cusps
};

CuspDetection detect_cusps(const HalfEdgeMesh& mesh, const CuspParams& params);

// --- Building blocks (exposed for testing) -------------------------------------------------

// Unit occlusal axis: the smallest-variance principal axis of the vertices in `mask` (an arch is
// wide and long but shallow), signed to point away from the boundary (the gingival cut).
// (uint8_t, not bool: std::vector<bool> is bit-packed and cannot back a span.)
Vec3 occlusal_axis(const HalfEdgeMesh& mesh, std::span<const std::uint8_t> mask);

// Implicit heat smoothing of a vertex field to spatial scale sigma (one backward-Euler step,
// t = sigma^2 / 2). The solution is O(1) everywhere and used directly, so CG is adequate here
// (contrast with the heat-method step, D47).
std::vector<double> diffuse(const DecOperators& ops, std::span<const double> field, double sigma,
                            double tolerance = 1e-8);

// Vertices of the largest connected component (1) vs the rest (0).
std::vector<std::uint8_t> largest_component_mask(const HalfEdgeMesh& mesh);

// Vertices whose score exceeds `threshold` and is strictly greater than every one-ring neighbor's.
std::vector<std::uint32_t> local_maxima(const HalfEdgeMesh& mesh, std::span<const double> score, double threshold);

// Local prominence p = f - diffuse(f, sigma), computed by solving for p directly (D70):
//   (*0 - tL) base = *0 f  with base = f - p   =>   (*0 - tL) p = -t L f.
// Mathematically identical to subtracting the diffused field, but the right-hand side is small and
// smooth (for f = height, -tLf is a large-scale mean curvature term), so CG converges in far fewer
// iterations than when reconstructing the whole field.
std::vector<double> prominence(const DecOperators& ops, std::span<const double> field, double sigma,
                               double tolerance = 1e-8);

// Keep candidates in descending score; drop any within `radius` of an already kept one.
std::vector<std::uint32_t> non_max_suppression(std::span<const Vec3> positions, std::span<const double> score,
                                               std::span<const std::uint32_t> candidates, double radius);

}  // namespace dmw
