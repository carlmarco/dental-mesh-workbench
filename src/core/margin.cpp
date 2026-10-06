#include "core/margin.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <queue>
#include <unordered_map>
#include <vector>

#include "core/curvature.h"
#include "core/dec.h"
#include "core/cusps.h"
#include "detail/hash.h"
#include "detail/vec.h"

namespace dmw {
namespace {

using namespace detail;

// Distance from each query point to its nearest point in `target`, via a uniform grid. Rings of cells
// are searched outward until the ring's inner distance exceeds the best distance found.
std::vector<double> nearest_distances(std::span<const Vec3> query, std::span<const Vec3> target, double cell) {
    std::vector<double> out(query.size(), std::numeric_limits<double>::infinity());
    if (target.empty()) return out;
    using Key = std::array<std::int64_t, 3>;
    auto key = [cell](const Vec3& p) {
        return Key{static_cast<std::int64_t>(std::floor(p.x / cell)), static_cast<std::int64_t>(std::floor(p.y / cell)),
                   static_cast<std::int64_t>(std::floor(p.z / cell))};
    };
    std::unordered_map<Key, std::vector<std::uint32_t>, Hash3<std::int64_t>> grid;
    for (std::uint32_t i = 0; i < target.size(); ++i) grid[key(target[i])].push_back(i);
    for (std::size_t q = 0; q < query.size(); ++q) {
        const Key c = key(query[q]);
        double best2 = std::numeric_limits<double>::infinity();
        for (std::int64_t ring = 0; ring < 100000; ++ring) {
            // Any point outside rings 0..ring is at least ring * cell away.
            const double inner = static_cast<double>(ring) * cell;
            if (ring > 0 && inner * inner > best2) break;
            for (std::int64_t dx = -ring; dx <= ring; ++dx)
                for (std::int64_t dy = -ring; dy <= ring; ++dy)
                    for (std::int64_t dz = -ring; dz <= ring; ++dz) {
                        if (std::max({std::abs(dx), std::abs(dy), std::abs(dz)}) != ring) continue;  // shell only
                        const auto it = grid.find({c[0] + dx, c[1] + dy, c[2] + dz});
                        if (it == grid.end()) continue;
                        for (std::uint32_t i : it->second) {
                            const Vec3 d = target[i] - query[q];
                            best2 = std::min(best2, dot(d, d));
                        }
                    }
        }
        out[q] = std::sqrt(best2);
    }
    return out;
}

}  // namespace

const std::vector<std::string>& seed_feature_names() {
    static const std::vector<std::string> names{"prominence_mm", "height_quantile", "normal_dot_axis",
                                                "distance_to_cut_mm", "smoothed_mean_curvature", "smoothed_gaussian_curvature"};
    return names;
}

std::vector<Edge> label_boundary_edges(const HalfEdgeMesh& m, std::span<const std::uint8_t> label) {
    std::vector<Edge> out;
    for (std::uint32_t h = 0; h < m.origin.size(); ++h) {
        const std::uint32_t t = m.twin[h];
        if (t != kInvalid && t < h) continue;  // each interior edge once
        const std::uint32_t a = m.origin[h], b = dest(m, h);
        if (label[a] != label[b]) out.push_back({std::min(a, b), std::max(a, b)});
    }
    return out;
}

std::vector<Vec3> edge_midpoints(const HalfEdgeMesh& m, std::span<const Edge> edges) {
    std::vector<Vec3> out;
    out.reserve(edges.size());
    for (const Edge& e : edges) out.push_back(0.5 * (m.positions[e.v0] + m.positions[e.v1]));
    return out;
}

BoundaryMetrics compare_boundaries(std::span<const Vec3> predicted, std::span<const Vec3> truth) {
    BoundaryMetrics r;
    r.predicted = predicted.size();
    r.truth = truth.size();
    if (predicted.empty() || truth.empty()) {
        r.assd = r.hd95 = r.hausdorff = std::numeric_limits<double>::infinity();
        return r;
    }
    const auto pt = nearest_distances(predicted, truth, 0.5);
    const auto tp = nearest_distances(truth, predicted, 0.5);
    double sum = 0.0;
    std::size_t p_in025 = 0, p_in050 = 0, t_in025 = 0, t_in050 = 0;
    for (double d : pt) sum += d, p_in025 += d <= 0.25, p_in050 += d <= 0.5;
    for (double d : tp) sum += d, t_in025 += d <= 0.25, t_in050 += d <= 0.5;
    r.assd = sum / static_cast<double>(pt.size() + tp.size());
    auto p95 = [](std::vector<double> v) { std::sort(v.begin(), v.end()); return v[static_cast<std::size_t>(0.95 * double(v.size() - 1))]; };
    r.hd95 = std::max(p95(pt), p95(tp));
    r.hausdorff = std::max(*std::max_element(pt.begin(), pt.end()), *std::max_element(tp.begin(), tp.end()));
    auto f1 = [&](std::size_t pin, std::size_t tin) {
        const double prec = double(pin) / double(pt.size()), rec = double(tin) / double(tp.size());
        return prec + rec > 0.0 ? 2.0 * prec * rec / (prec + rec) : 0.0;
    };
    r.f1_025 = f1(p_in025, t_in025);
    r.f1_050 = f1(p_in050, t_in050);
    return r;
}

double region_iou(std::span<const std::uint8_t> predicted, std::span<const std::uint8_t> truth, std::span<const double> area) {
    double inter = 0.0, uni = 0.0;
    for (std::size_t v = 0; v < predicted.size(); ++v) {
        if (predicted[v] && truth[v]) inter += area[v];
        if (predicted[v] || truth[v]) uni += area[v];
    }
    return uni > 0.0 ? inter / uni : 1.0;
}

MarginInputs margin_inputs(const HalfEdgeMesh& m, bool with_cusps) {
    const std::size_t nv = m.positions.size();
    MarginInputs in;
    in.arch = largest_component_mask(m);
    const Vec3 axis = occlusal_axis(m, in.arch);
    in.height.resize(nv);
    for (std::size_t v = 0; v < nv; ++v) {
        in.height[v] = dot(m.positions[v], axis);
        if (in.arch[v]) in.sorted_arch_height.push_back(in.height[v]);
    }
    std::sort(in.sorted_arch_height.begin(), in.sorted_arch_height.end());
    for (std::uint32_t h = 0; h < m.origin.size(); ++h) {
        if (m.twin[h] == kInvalid && in.arch[m.origin[h]]) in.cut_vertices.push_back(m.origin[h]);
    }
    const CurvatureField curv = compute_curvature(m);
    in.kmin.resize(nv);
    for (std::size_t v = 0; v < nv; ++v) in.kmin[v] = std::isfinite(curv.k2[v]) ? curv.k2[v] : 0.0;
    if (with_cusps) {
        const CuspDetection cusps = detect_cusps(m, cusp_operating_point());
        in.cusp_tips = cusps.vertices;
        // Seed features (D76).
        const DecOperators ops = build_dec(m);
        std::vector<double> mean(nv), gauss(nv);
        for (std::size_t v = 0; v < nv; ++v) {
            mean[v] = std::isfinite(curv.mean[v]) ? curv.mean[v] : 0.0;
            gauss[v] = std::isfinite(curv.gaussian[v]) ? curv.gaussian[v] : 0.0;
        }
        const auto smooth_mean = diffuse(ops, mean, 0.5, 1e-4), smooth_gauss = diffuse(ops, gauss, 0.5, 1e-4);
        std::vector<Vec3> normal(nv, Vec3{});  // area-weighted vertex normals
        for (std::size_t f = 0; f < m.origin.size() / 3; ++f) {
            const std::uint32_t a = m.origin[3 * f], b = m.origin[3 * f + 1], c = m.origin[3 * f + 2];
            const Vec3 n = cross(m.positions[b] - m.positions[a], m.positions[c] - m.positions[a]);
            normal[a] += n, normal[b] += n, normal[c] += n;
        }
        std::vector<Vec3> tips, cut;
        for (std::uint32_t v : in.cusp_tips) tips.push_back(m.positions[v]);
        for (std::uint32_t v : in.cut_vertices) cut.push_back(m.positions[v]);
        const auto to_cut = nearest_distances(tips, cut, 1.0);
        const auto& hs = in.sorted_arch_height;
        for (std::size_t i = 0; i < in.cusp_tips.size(); ++i) {
            const std::uint32_t v = in.cusp_tips[i];
            const double nn = norm(normal[v]);
            const double quant = hs.empty() ? 0.0
                                            : double(std::lower_bound(hs.begin(), hs.end(), in.height[v]) - hs.begin()) / double(hs.size());
            in.seed_features.insert(in.seed_features.end(),
                                    {cusps.scores[i], quant, nn > 0.0 ? dot(normal[v], axis) / nn : 0.0,
                                     std::isfinite(to_cut[i]) ? to_cut[i] : 100.0, smooth_mean[v], smooth_gauss[v]});
        }
    }
    return in;
}

std::vector<double> valley_strength(const DecOperators& ops, std::span<const double> kmin, double sigma) {
    std::vector<double> s = diffuse(ops, kmin, sigma, 1e-4);  // tolerance as in D70
    for (double& x : s) x = std::max(0.0, -x);
    return s;
}

std::vector<std::uint8_t> margin_labels(const HalfEdgeMesh& m, const MarginInputs& in, std::span<const double> valley,
                                        const MarginParams& p) {
    const std::size_t nv = m.positions.size();
    const auto& hs = in.sorted_arch_height;
    auto q = [&](double f) { return hs.empty() ? 0.0 : hs[static_cast<std::size_t>(f * double(hs.size() - 1))]; };
    std::vector<std::uint8_t> tooth(nv, 0);
    if (p.method == MarginParams::Method::HeightPlane) {
        const double cut = q(p.plane_quantile);
        for (std::size_t v = 0; v < nv; ++v) tooth[v] = (in.arch[v] && in.height[v] >= cut) ? 1 : 0;
        return tooth;
    }
    // Seeds: 1 tooth, 0 gingiva, 2 unlabelled.
    std::vector<std::uint8_t> label(nv, 2);
    const double low = q(p.gingiva_quantile), high = q(p.tooth_quantile);
    for (std::size_t v = 0; v < nv; ++v) {
        if (!in.arch[v]) continue;
        if (in.height[v] <= low) label[v] = 0;
        else if (in.height[v] >= high) label[v] = 1;
    }
    for (std::uint32_t v : in.cut_vertices) label[v] = 0;
    if (p.cusp_seeds) {
        const double seed_gate = q(p.cusp_seed_quantile);
        for (std::size_t i = 0; i < in.cusp_tips.size(); ++i) {
            const std::uint32_t v = in.cusp_tips[i];
            if (p.cusp_seed_quantile > 0.0 && in.height[v] < seed_gate) continue;
            if (p.seed_model && in.seed_features.size() >= (i + 1) * kSeedFeatureCount &&
                p.seed_model->probability(std::span<const double>(in.seed_features).subspan(i * kSeedFeatureCount, kSeedFeatureCount)) <
                    p.seed_threshold) {
                continue;  // the classifier says this seed is probably on gingiva (D76)
            }
            label[v] = 1;
        }
    }
    // Two-label multi-source Dijkstra: each vertex takes the label of the front that arrives first.
    std::vector<double> dist(nv, std::numeric_limits<double>::infinity());
    using Item = std::pair<double, std::uint32_t>;
    std::priority_queue<Item, std::vector<Item>, std::greater<>> heap;
    for (std::uint32_t v = 0; v < nv; ++v) {
        if (label[v] != 2) dist[v] = 0.0, heap.push({0.0, v});
    }
    const bool weighted = p.valley_weight > 0.0 && !valley.empty();
    while (!heap.empty()) {
        const auto [d, v] = heap.top();
        heap.pop();
        if (d > dist[v]) continue;
        for (std::uint32_t w : one_ring(m, v)) {
            if (!in.arch[w]) continue;
            double cost = norm(m.positions[w] - m.positions[v]);
            if (weighted) cost *= 1.0 + p.valley_weight * 0.5 * (valley[v] + valley[w]);
            if (d + cost < dist[w]) {
                dist[w] = d + cost;
                label[w] = label[v];
                heap.push({dist[w], w});
            }
        }
    }
    for (std::size_t v = 0; v < nv; ++v) tooth[v] = (in.arch[v] && label[v] == 1) ? 1 : 0;
    return tooth;
}

MarginParams margin_operating_point() {
    // Chosen on 60 training scans (D73): ASSD 0.647 mm. sigma = 0 means raw kappa_min (the cervical
    // crease is narrow; smoothing blurs it); alpha in its converged regime (cost ~ integrated concavity).
    MarginParams p;
    p.method = MarginParams::Method::GeodesicVoronoi;
    p.valley_weight = 2560.0, p.curvature_scale = 0.0, p.gingiva_quantile = 0.15, p.tooth_quantile = 0.9, p.cusp_seeds = true;
    return p;
}

MarginResult detect_margin(const HalfEdgeMesh& m, const MarginParams& p) {
    const bool voronoi = p.method == MarginParams::Method::GeodesicVoronoi;
    const MarginInputs in = margin_inputs(m, voronoi && p.cusp_seeds);
    std::vector<double> valley;
    if (voronoi && p.valley_weight > 0.0) valley = valley_strength(build_dec(m), in.kmin, p.curvature_scale);
    MarginResult out;
    out.tooth = margin_labels(m, in, valley, p);
    out.margin = label_boundary_edges(m, out.tooth);
    return out;
}

}  // namespace dmw
