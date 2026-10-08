#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "core/dec.h"
#include "core/halfedge.h"
#include "core/learn.h"
#include "core/topology.h"

namespace dmw {

// Tooth-gingiva margin detection (M9, D71). A binary tooth/gingiva labelling of the arch; the margin
// line is its boundary (label_boundary_edges), the same definition used for the ground truth.
//   HeightPlane      tooth = height along the occlusal axis above a quantile (naive baseline)
//   GeodesicVoronoi  multi-source Dijkstra from tooth seeds (detected cusp tips + highest vertices) and
//                    gingiva seeds (the scan's cut boundary + lowest vertices); each vertex takes the label
//                    of the front that reaches it first. Edge cost l * (1 + valley_weight * s), where
//                    s = max(0, -kappa_min) of the smoothed minimum principal curvature (valley strength):
//                    fronts stall in concave creases, so they meet in the cervical crease.
//   GraphCut         same seeds and edge costs; minimum s-t cut (D78) of the energy
//                      E = sum_v A_v D_v(label) + mu * sum_{cut edges} l*_e / (1 + beta * s_e)
//                    with D from the two arrival distances (p_v = d_G / (d_T + d_G) is the tooth
//                    likelihood, D = -log p or -log(1 - p)), l*_e the cotan dual edge length, and seeds
//                    as hard constraints. A fake tooth region on flat gingiva pays mu per mm of its
//                    boundary, so it survives only if its unary evidence outweighs its perimeter.
//   PerToothCut      (D85) the GraphCut result, then cusp seeds grouped into teeth (group_cusps on its tooth
//                    region) and a multi-label cut (gingiva + one label per group) by alpha-expansion: unary
//                    p_l proportional to 1 / d_l (reduces to the GraphCut unary for two labels), the same
//                    boundary costs (Potts). Two touching crowns now need a boundary between them; where it runs
//                    down a valley (kappa_min < strip_kmin) its deeper vertex becomes gingiva, the one-vertex
//                    interdental strip of the labels (D84).
// Only the largest component (the arch) is labelled; everything else is gingiva.
struct MarginParams {
    enum class Method { HeightPlane, GeodesicVoronoi, GraphCut, PerToothCut };
    Method method = Method::GeodesicVoronoi;
    double plane_quantile = 0.5;     // HeightPlane: tooth above this height quantile
    double valley_weight = 10.0;     // alpha (mm); 0 = plain geodesic Voronoi
    double curvature_scale = 0.3;    // sigma (mm) for smoothing kappa_min
    double gingiva_quantile = 0.15;  // gingiva seeds: arch vertices below this height quantile
    double tooth_quantile = 0.9;     // tooth seeds: arch vertices above this height quantile
    double cut_smoothness = 1.0;     // GraphCut: mu (mm), boundary cost per mm of cut on flat surface
    double cut_crease = 100.0;       // GraphCut: beta (mm), how much cheaper a cut is along a crease
    double group_crease = -1.5;      // PerToothCut: tau for group_cusps (1/mm); >= 1e8 = no grouping (one label per tip)
    double strip_kmin = -1e300;      // PerToothCut: carve tooth|tooth boundary vertices with kappa_min below this (off:
                                     // every carving variant was worse on the sweep scans, D85)
    int expansion_sweeps = 3;        // PerToothCut: alpha-expansion sweeps over all labels (stops early if stable)
    double expansion_radius = 3.0;   // PerToothCut: a move to label l only considers vertices within ~this (mm, grid
                                     // cells) of vertices labelled l (0 = all vertices; D85 speed)
    double min_tooth_region = 0.0;   // relabel tooth components smaller than this area (mm^2) as gingiva (D80)
    bool cusp_seeds = true;          // also seed teeth at detected cusp tips (detect_cusps, D68 point)
    double cusp_seed_quantile = 0.0; // keep only cusp seeds above this arch-height quantile (D75; 0 = all)
    // Seed classifier (D76): drop cusp seeds with P(tooth) < seed_threshold. Non-owning; null = off.
    const LogisticModel* seed_model = nullptr;
    double seed_threshold = 0.5;
    // DIAGNOSTICS ONLY (oracle experiments, D77): vertices flagged 1 may not seed teeth. Null = off.
    const std::vector<std::uint8_t>* tooth_seed_veto = nullptr;
};

struct MarginResult {
    std::vector<std::uint8_t> tooth;  // per vertex: 1 tooth, 0 gingiva
    std::vector<Edge> margin;         // boundary edges between tooth and gingiva
};

MarginResult detect_margin(const HalfEdgeMesh& mesh, const MarginParams& params);

// The operating point chosen on training scans (D73). Single source of truth for viewer and tools.
MarginParams margin_operating_point();

// --- Building blocks (detect_margin composes these; the evaluation tool caches them per scan) ---
struct MarginInputs {
    std::vector<std::uint8_t> arch;          // largest component
    std::vector<double> height;              // along the occlusal axis
    std::vector<double> sorted_arch_height;  // ascending, arch vertices only
    std::vector<std::uint32_t> cut_vertices; // boundary of the arch (the scan's cut through the gingiva)
    std::vector<std::uint32_t> cusp_tips;    // detect_cusps at the D68 operating point (if requested)
    std::vector<double> seed_features;       // kSeedFeatureCount per cusp tip, row-major (D76)
    std::vector<double> kmin;                // minimum principal curvature (0 where undefined)
};
MarginInputs margin_inputs(const HalfEdgeMesh& mesh, bool with_cusps);

// Label-free features of each cusp seed for the seed classifier (D76), in this order:
//   prominence (mm), height quantile in the arch (0..1), vertex normal . occlusal axis,
//   distance to the scan's cut boundary (mm), smoothed mean curvature (1/mm), smoothed Gaussian (1/mm^2).
inline constexpr std::size_t kSeedFeatureCount = 6;
const std::vector<std::string>& seed_feature_names();

// s = max(0, -kappa_min) after smoothing kappa_min to scale sigma.
std::vector<double> valley_strength(const DecOperators& ops, std::span<const double> kmin, double sigma);

// The labelling step for either method, given precomputed inputs (valley ignored for HeightPlane).
std::vector<std::uint8_t> margin_labels(const HalfEdgeMesh& mesh, const MarginInputs& in, std::span<const double> valley,
                                        const MarginParams& params);

// Group cusp tips into teeth (D85): two tips share a group iff they are connected inside `region` through
// vertices with kappa_min > tau (no crease deeper than tau on the way). Adjacent crowns meet in a deep valley
// (D84), while cusps of one crown are joined by ridges. Returns a group id per tip (0..G-1), or -1 for a tip
// outside the region or itself at kappa_min <= tau.
std::vector<int> group_cusps(const HalfEdgeMesh& mesh, std::span<const std::uint32_t> tips, std::span<const std::uint8_t> region,
                             std::span<const double> kmin, double tau);

// Edges whose endpoints carry different binary labels (the margin line for tooth/gingiva labels).
std::vector<Edge> label_boundary_edges(const HalfEdgeMesh& mesh, std::span<const std::uint8_t> label);

// Boundary-to-boundary comparison in mm (D72): margin edges are sampled at their midpoints, and each
// sample is matched to its nearest sample on the other line (spatial hash).
struct BoundaryMetrics {
    double assd = 0.0;           // average symmetric surface distance
    double hd95 = 0.0;           // max of the two directed 95th-percentile distances
    double hausdorff = 0.0;      // max of the two directed max distances
    double f1_025 = 0.0, f1_050 = 0.0;  // boundary F1 at 0.25 mm / 0.5 mm
    std::size_t predicted = 0, truth = 0;  // sample counts
};
BoundaryMetrics compare_boundaries(std::span<const Vec3> predicted, std::span<const Vec3> truth);

std::vector<Vec3> edge_midpoints(const HalfEdgeMesh& mesh, std::span<const Edge> edges);

// Distance from each query point to its nearest target point (uniform grid of cell size `cell`,
// ring search). Infinity when `target` is empty.
std::vector<double> nearest_point_distances(std::span<const Vec3> query, std::span<const Vec3> target, double cell = 0.5);

// Area-weighted IoU of the label-1 regions (weights: lumped vertex areas).
double region_iou(std::span<const std::uint8_t> predicted, std::span<const std::uint8_t> truth, std::span<const double> area);

}  // namespace dmw
