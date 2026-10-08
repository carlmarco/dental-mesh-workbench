#include "core/margin.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <functional>
#include <limits>
#include <queue>
#include <unordered_map>
#include <vector>

#include "core/curvature.h"
#include "core/dec.h"
#include "core/cusps.h"
#include "core/maxflow.h"
#include "core/multilabel.h"
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

// Island removal (D80): connected tooth regions (over mesh edges) with lumped area below `min_area`
// become gingiva. A real crown is tens of mm^2; specks around false seeds are far smaller.
void remove_small_tooth_regions(const HalfEdgeMesh& m, std::vector<std::uint8_t>& tooth, double min_area) {
    if (min_area <= 0.0) return;
    const std::size_t nv = m.positions.size();
    std::vector<double> area(nv, 0.0);
    for (std::size_t f = 0; f < m.origin.size() / 3; ++f) {
        const std::uint32_t a = m.origin[3 * f], b = m.origin[3 * f + 1], c = m.origin[3 * f + 2];
        const double t = 0.5 * norm(cross(m.positions[b] - m.positions[a], m.positions[c] - m.positions[a])) / 3.0;
        area[a] += t, area[b] += t, area[c] += t;
    }
    std::vector<std::uint8_t> seen(nv, 0);
    std::vector<std::uint32_t> region;
    for (std::uint32_t s = 0; s < nv; ++s) {
        if (!tooth[s] || seen[s]) continue;
        region.assign(1, s);
        seen[s] = 1;
        double total = 0.0;
        for (std::size_t k = 0; k < region.size(); ++k) {
            total += area[region[k]];
            for (std::uint32_t w : one_ring(m, region[k]))
                if (tooth[w] && !seen[w]) seen[w] = 1, region.push_back(w);
        }
        if (total < min_area)
            for (std::uint32_t v : region) tooth[v] = 0;
    }
}

}  // namespace

const std::vector<std::string>& seed_feature_names() {
    static const std::vector<std::string> names{"prominence_mm", "height_quantile", "normal_dot_axis",
                                                "distance_to_cut_mm", "smoothed_mean_curvature", "smoothed_gaussian_curvature"};
    return names;
}

std::vector<double> nearest_point_distances(std::span<const Vec3> query, std::span<const Vec3> target, double cell) {
    return nearest_distances(query, target, cell);
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
    using Clock = std::chrono::steady_clock;
    auto t_last = Clock::now();
    auto lap = [&](double MarginTimings::*field) {
        if (!p.timings) return;
        const auto now = Clock::now();
        p.timings->*field += std::chrono::duration<double, std::milli>(now - t_last).count();
        t_last = now;
    };
    if (p.timings) ++p.timings->calls;
    // Seeds: 1 tooth, 0 gingiva, 2 unlabelled.
    std::vector<std::uint8_t> label(nv, 2);
    const double low = q(p.gingiva_quantile), high = q(p.tooth_quantile);
    for (std::size_t v = 0; v < nv; ++v) {
        if (!in.arch[v]) continue;
        if (in.height[v] <= low) label[v] = 0;
        else if (in.height[v] >= high) label[v] = 1;
    }
    for (std::uint32_t v : in.cut_vertices) label[v] = 0;
    auto vetoed = [&](std::uint32_t v) { return p.tooth_seed_veto && (*p.tooth_seed_veto)[v]; };
    for (std::uint32_t v = 0; v < nv; ++v) {
        if (label[v] == 1 && vetoed(v)) label[v] = 2;  // oracle: height seed on true gingiva removed
    }
    std::vector<std::uint32_t> seed_tips;  // cusp tips that became tooth seeds (for PerToothCut grouping)
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
            if (vetoed(v)) continue;
            label[v] = 1;
            seed_tips.push_back(v);
        }
    }
    const bool weighted = p.valley_weight > 0.0 && !valley.empty();
    auto edge_cost = [&](std::uint32_t v, std::uint32_t w) {
        double cost = norm(m.positions[w] - m.positions[v]);
        if (weighted) cost *= 1.0 + p.valley_weight * 0.5 * (valley[v] + valley[w]);
        return cost;
    };
    // Multi-source Dijkstra over the arch from the vertices with label[v] in `sources`; `owner` records
    // which seed label each vertex was reached from (first arrival).
    // `within(w)`: an optional spatial limit on which vertices the search may visit (PerToothCut, D86).
    auto arrival = [&](auto is_source, std::vector<double>& dist, std::vector<std::uint8_t>* owner,
                       const std::function<bool(std::uint32_t)>& within = {}, std::vector<std::uint32_t>* parent = nullptr,
                       bool plain = false) {
        dist.assign(nv, std::numeric_limits<double>::infinity());
        if (parent) parent->assign(nv, kInvalid);
        using Item = std::pair<double, std::uint32_t>;
        std::priority_queue<Item, std::vector<Item>, std::greater<>> heap;
        for (std::uint32_t v = 0; v < nv; ++v) {
            if (is_source(label[v])) {
                dist[v] = 0.0, heap.push({0.0, v});
                if (parent) (*parent)[v] = v;
            }
        }
        while (!heap.empty()) {
            const auto [d, v] = heap.top();
            heap.pop();
            if (d > dist[v]) continue;
            for (std::uint32_t w : one_ring(m, v)) {
                if (!in.arch[w] || (within && !within(w))) continue;
                const double nd = d + (plain ? norm(m.positions[w] - m.positions[v]) : edge_cost(v, w));
                if (nd < dist[w]) {
                    dist[w] = nd;
                    if (owner) (*owner)[w] = (*owner)[v];
                    if (parent) (*parent)[w] = v;
                    heap.push({nd, w});
                }
            }
        }
    };
    if (p.method == MarginParams::Method::GeodesicVoronoi) {
        // Two-label multi-source Dijkstra: each vertex takes the label of the front that arrives first.
        std::vector<double> dist;
        arrival([](std::uint8_t l) { return l != 2; }, dist, &label);
        for (std::size_t v = 0; v < nv; ++v) tooth[v] = (in.arch[v] && label[v] == 1) ? 1 : 0;
        remove_small_tooth_regions(m, tooth, p.min_tooth_region);
        return tooth;
    }

    // GraphCut (D78). Unary terms from the two arrival distances; nodes are the arch vertices.
    lap(&MarginTimings::seeds);
    std::vector<double> dt, dg;
    arrival([](std::uint8_t l) { return l == 1; }, dt, nullptr);
    arrival([](std::uint8_t l) { return l == 0; }, dg, nullptr);
    lap(&MarginTimings::binary_dijkstra);
    std::vector<double> area(nv, 0.0), dual(m.origin.size(), 0.0);  // lumped areas; dual length per half-edge
    for (std::uint32_t f = 0; f < m.origin.size() / 3; ++f) {
        const std::array<std::uint32_t, 3> c{m.origin[3 * f], m.origin[3 * f + 1], m.origin[3 * f + 2]};
        const double a = 0.5 * norm(cross(m.positions[c[1]] - m.positions[c[0]], m.positions[c[2]] - m.positions[c[0]]));
        for (std::uint32_t k = 0; k < 3; ++k) {
            area[c[k]] += a / 3.0;
            // Half-edge 3f + k runs c[k] -> c[k+1]; its opposite corner is c[k+2]. Contribution to the
            // dual edge length: cot(opposite angle) / 2 * |e| (the *1 entry of D45 times the length).
            const Vec3 u = m.positions[c[k]] - m.positions[c[(k + 2) % 3]];
            const Vec3 w = m.positions[c[(k + 1) % 3]] - m.positions[c[(k + 2) % 3]];
            const double sin2 = norm(cross(u, w));
            const double cot = sin2 > 0.0 ? dot(u, w) / sin2 : 0.0;
            dual[3 * f + k] = 0.5 * cot * norm(m.positions[c[(k + 1) % 3]] - m.positions[c[k]]);
        }
    }
    std::vector<std::uint32_t> node(nv, kInvalid);
    std::uint32_t count = 0;
    for (std::uint32_t v = 0; v < nv; ++v) {
        if (in.arch[v]) node[v] = count++;
    }
    MaxFlow g(count);
    constexpr double kHard = 1e12, kEps = 1e-6;
    for (std::uint32_t v = 0; v < nv; ++v) {
        if (!in.arch[v]) continue;
        if (label[v] == 1) {
            g.add_terminal(node[v], kHard, 0.0);
        } else if (label[v] == 0) {
            g.add_terminal(node[v], 0.0, kHard);
        } else {
            // Tooth likelihood p = d_G / (d_T + d_G): 1 at tooth seeds, 1/2 where the fronts meet.
            const double t = dt[v], s = dg[v];
            const double pt = !std::isfinite(t) ? 0.0 : !std::isfinite(s) ? 1.0 : s / (s + t);
            // Source side = tooth pays the sink capacity: D(tooth) = -A log p; D(gingiva) = -A log(1 - p).
            g.add_terminal(node[v], -area[v] * std::log(std::max(1.0 - pt, kEps)), -area[v] * std::log(std::max(pt, kEps)));
        }
    }
    if (p.cut_smoothness > 0.0) {
        for (std::uint32_t h = 0; h < m.origin.size(); ++h) {
            const std::uint32_t t = m.twin[h];
            if (t == kInvalid || t < h) continue;  // interior edges once; boundary edges cut nothing
            const std::uint32_t a = m.origin[h], b = dest(m, h);
            if (!in.arch[a] || !in.arch[b]) continue;
            // Obtuse pairs can make the cotan dual length negative; clamp (a cut cost must be >= 0).
            const double len = std::max(0.0, dual[h] + dual[t]);
            const double crease = valley.empty() ? 0.0 : 0.5 * (valley[a] + valley[b]);
            const double w = p.cut_smoothness * len / (1.0 + p.cut_crease * crease);
            if (w > 0.0) g.add_edge(node[a], node[b], w, w);
        }
    }
    g.solve();
    for (std::uint32_t v = 0; v < nv; ++v) tooth[v] = (in.arch[v] && g.source_side(node[v])) ? 1 : 0;
    lap(&MarginTimings::binary_cut);
    if (p.method != MarginParams::Method::PerToothCut) {
        remove_small_tooth_regions(m, tooth, p.min_tooth_region);
        return tooth;
    }

    // --- PerToothCut (D85): group the seed tips into teeth on the binary tooth region, then alpha-expansion.
    std::vector<int> grp;
    if (p.group_crease >= 1e8) {  // no grouping: every seed tip in the tooth region is its own label
        for (std::uint32_t v : seed_tips) grp.push_back(tooth[v] ? static_cast<int>(grp.size()) : -1);
        int next = 0;
        for (int& x : grp) x = x < 0 ? -1 : next++;
    } else {
        grp = group_cusps(m, seed_tips, tooth, in.kmin, p.group_crease);
    }
    const int groups = grp.empty() ? 0 : *std::max_element(grp.begin(), grp.end()) + 1;
    lap(&MarginTimings::grouping);
    if (groups < 2) return tooth;  // nothing to separate
    const auto labels = static_cast<std::size_t>(groups) + 1;  // 0 = gingiva, 1..groups = teeth
    std::vector<std::vector<double>> dist(labels);
    std::vector<std::vector<std::uint32_t>> tree(labels);  // per-label shortest-path parents (vertex ids; star prior)
    dist[0] = dg;
    std::vector<int> fixed(nv, -1);  // hard constraints
    for (std::uint32_t v = 0; v < nv; ++v)
        if (in.arch[v] && label[v] == 0) fixed[v] = 0;
    for (int l = 1; l <= groups; ++l) {
        std::vector<std::uint8_t> src(nv, 0);
        for (std::size_t i = 0; i < seed_tips.size(); ++i)
            if (grp[i] == l - 1) src[seed_tips[i]] = 1, fixed[seed_tips[i]] = l;
        std::vector<std::uint8_t> keep = label;
        for (std::uint32_t v = 0; v < nv; ++v) label[v] = src[v] ? 1 : 2;
        std::function<bool(std::uint32_t)> near_seeds;
        if (p.label_radius > 0.0) {
            std::vector<Vec3> tips_l;
            for (std::size_t i = 0; i < seed_tips.size(); ++i)
                if (grp[i] == l - 1) tips_l.push_back(m.positions[seed_tips[i]]);
            const double r2 = p.label_radius * p.label_radius;
            near_seeds = [&m, tips_l, r2](std::uint32_t w) {
                for (const Vec3& t : tips_l) {
                    const Vec3 d = m.positions[w] - t;
                    if (dot(d, d) <= r2) return true;
                }
                return false;
            };
        }
        arrival([](std::uint8_t x) { return x == 1; }, dist[static_cast<std::size_t>(l)], nullptr, near_seeds,
                p.star_prior && !p.star_plain ? &tree[static_cast<std::size_t>(l)] : nullptr);
        if (p.star_prior && p.star_plain) {  // straight-ray trees: a separate search over plain edge lengths
            std::vector<double> plain_dist;
            arrival([](std::uint8_t x) { return x == 1; }, plain_dist, nullptr, near_seeds, &tree[static_cast<std::size_t>(l)], true);
        }
        label = std::move(keep);
    }
    lap(&MarginTimings::label_dijkstra);
    if (p.timings) p.timings->labels += labels;
    // Unary D_l(v) = -A log p_l, p_l = (1/d_l) / sum_k (1/d_k); hard constraints as kHard.
    std::vector<double> unary(static_cast<std::size_t>(count) * labels, 0.0);
    for (std::uint32_t v = 0; v < nv; ++v) {
        if (!in.arch[v]) continue;
        double* u = &unary[static_cast<std::size_t>(node[v]) * labels];
        if (fixed[v] >= 0) {
            for (std::size_t l = 0; l < labels; ++l) u[l] = static_cast<int>(l) == fixed[v] ? 0.0 : kHard;
            continue;
        }
        double total = 0.0;
        std::vector<double> inv(labels);
        for (std::size_t l = 0; l < labels; ++l) {
            const double d = dist[l][v];
            inv[l] = std::isfinite(d) ? 1.0 / std::max(d, 1e-12) : 0.0;
            total += inv[l];
        }
        for (std::size_t l = 0; l < labels; ++l) u[l] = -area[v] * std::log(std::max(total > 0.0 ? inv[l] / total : 0.0, kEps));
    }
    // Pairwise Potts weights (same as the binary cut).
    std::vector<PottsEdge> pairs;
    if (p.cut_smoothness > 0.0) {
        for (std::uint32_t h = 0; h < m.origin.size(); ++h) {
            const std::uint32_t t = m.twin[h];
            if (t == kInvalid || t < h) continue;
            const std::uint32_t a = m.origin[h], b = dest(m, h);
            if (!in.arch[a] || !in.arch[b]) continue;
            const double len = std::max(0.0, dual[h] + dual[t]);
            const double crease = valley.empty() ? 0.0 : 0.5 * (valley[a] + valley[b]);
            const double w = p.cut_smoothness * len / (1.0 + p.cut_crease * crease);
            if (w > 0.0) pairs.push_back({node[a], node[b], w});
        }
    }
    // Initial labelling: gingiva where the binary cut said gingiva, else the most likely tooth label.
    std::vector<std::uint32_t> f(count, 0);
    for (std::uint32_t v = 0; v < nv; ++v) {
        if (!in.arch[v]) continue;
        const double* u = &unary[static_cast<std::size_t>(node[v]) * labels];
        std::size_t best = 0;
        if (tooth[v])
            for (std::size_t l = 1; l < labels; ++l)
                if (best == 0 || u[l] < u[best]) best = l;
        if (fixed[v] >= 0) best = static_cast<std::size_t>(fixed[v]);
        f[node[v]] = static_cast<std::uint32_t>(best);
    }
    lap(&MarginTimings::unary);
    ExpansionCandidates near_label;
    if (p.expansion_radius > 0.0) {
        // Dense grid of cells of size R over the arch: a vertex may take label l if its cell or one of the 26
        // around it holds a vertex labelled l (occupancy dilated once per move, then one lookup per vertex).
        const double r = p.expansion_radius;
        Vec3 lo{1e300, 1e300, 1e300};
        for (std::uint32_t v = 0; v < nv; ++v)
            if (in.arch[v]) lo = {std::min(lo.x, m.positions[v].x), std::min(lo.y, m.positions[v].y), std::min(lo.z, m.positions[v].z)};
        std::array<std::size_t, 3> dim{1, 1, 1};
        std::vector<std::array<std::size_t, 3>> ijk(count);
        for (std::uint32_t v = 0; v < nv; ++v) {
            if (!in.arch[v]) continue;
            const Vec3 d = m.positions[v] - lo;
            const std::array<std::size_t, 3> c{static_cast<std::size_t>(d.x / r) + 1, static_cast<std::size_t>(d.y / r) + 1,
                                               static_cast<std::size_t>(d.z / r) + 1};  // +1: a margin layer for dilation
            ijk[node[v]] = c;
            for (std::size_t k = 0; k < 3; ++k) dim[k] = std::max(dim[k], c[k] + 2);
        }
        std::vector<std::size_t> cell(count);
        for (std::uint32_t i = 0; i < count; ++i) cell[i] = (ijk[i][2] * dim[1] + ijk[i][1]) * dim[0] + ijk[i][0];
        near_label = [cell, count, dim](std::uint32_t alpha, const std::vector<std::uint32_t>& lab, std::vector<std::uint8_t>& mask) {
            const std::size_t cells = dim[0] * dim[1] * dim[2];
            std::vector<std::uint8_t> occupied(cells, 0), near(cells, 0);
            for (std::uint32_t i = 0; i < count; ++i)
                if (lab[i] == alpha) occupied[cell[i]] = 1;
            const auto sx = static_cast<std::ptrdiff_t>(1), sy = static_cast<std::ptrdiff_t>(dim[0]),
                       sz = static_cast<std::ptrdiff_t>(dim[0] * dim[1]);
            for (std::size_t c = 0; c < cells; ++c) {
                if (!occupied[c]) continue;
                for (std::ptrdiff_t dz = -1; dz <= 1; ++dz)
                    for (std::ptrdiff_t dy = -1; dy <= 1; ++dy)
                        for (std::ptrdiff_t dx = -1; dx <= 1; ++dx) {
                            const std::ptrdiff_t n = static_cast<std::ptrdiff_t>(c) + dx * sx + dy * sy + dz * sz;
                            if (n >= 0 && n < static_cast<std::ptrdiff_t>(cells)) near[static_cast<std::size_t>(n)] = 1;
                        }
            }
            for (std::uint32_t i = 0; i < count; ++i) mask[i] = near[cell[i]];
        };
    }
    StarParents star;
    if (p.star_prior) {
        // Node-indexed trees for the tooth labels (gingiva unconstrained); unreached nodes may not take the label.
        star.resize(labels);
        for (std::size_t l = 1; l < labels; ++l) {
            star[l].assign(count, kNoParent);
            for (std::uint32_t v = 0; v < nv; ++v) {
                if (!in.arch[v] || tree[l][v] == kInvalid) continue;
                star[l][node[v]] = node[tree[l][v]];
            }
        }
        // Repair the initial labelling: process nodes by distance in their own label's tree (parents first).
        std::vector<std::uint32_t> order(count), vertex_of(count);
        for (std::uint32_t v = 0; v < nv; ++v)
            if (in.arch[v]) vertex_of[node[v]] = v;
        for (std::uint32_t i = 0; i < count; ++i) order[i] = i;
        auto key = [&](std::uint32_t i) { return f[i] == 0 ? 0.0 : dist[f[i]][vertex_of[i]]; };
        std::sort(order.begin(), order.end(), [&](std::uint32_t a, std::uint32_t b) { return key(a) < key(b); });
        for (int pass = 0; pass < 4 && !star_feasible(star, f); ++pass) star_repair(star, order, 0, f);
    }
    alpha_expansion(static_cast<std::uint32_t>(labels), unary, pairs, f, p.expansion_sweeps, near_label, star);
    lap(&MarginTimings::expansion);
    for (std::uint32_t v = 0; v < nv; ++v) tooth[v] = (in.arch[v] && f[node[v]] != 0) ? 1 : 0;
    // Interdental strip: at tooth|tooth edges along a valley, the deeper vertex becomes gingiva.
    std::vector<std::uint32_t> carve;
    for (std::uint32_t h = 0; h < m.origin.size(); ++h) {
        const std::uint32_t t = m.twin[h];
        if (t == kInvalid || t < h) continue;
        const std::uint32_t a = m.origin[h], b = dest(m, h);
        if (!in.arch[a] || !in.arch[b]) continue;
        const std::uint32_t la = f[node[a]], lb = f[node[b]];
        if (la == 0 || lb == 0 || la == lb) continue;
        const std::uint32_t deeper = in.kmin[a] <= in.kmin[b] ? a : b;
        if (in.kmin[deeper] < p.strip_kmin && fixed[deeper] < 0) carve.push_back(deeper);
    }
    for (std::uint32_t v : carve) tooth[v] = 0;
    remove_small_tooth_regions(m, tooth, p.min_tooth_region);
    return tooth;
}

std::vector<int> group_cusps(const HalfEdgeMesh& m, std::span<const std::uint32_t> tips, std::span<const std::uint8_t> region,
                             std::span<const double> kmin, double tau) {
    const std::size_t nv = m.positions.size();
    auto open = [&](std::uint32_t v) { return region[v] && kmin[v] > tau; };
    std::vector<int> component(nv, -1);
    std::vector<int> group(tips.size(), -1);
    std::vector<int> group_of_component;  // component id -> group id (-1 until a tip lands in it)
    int components = 0;
    std::vector<std::uint32_t> queue;
    for (std::size_t i = 0; i < tips.size(); ++i) {
        const std::uint32_t t = tips[i];
        if (!open(t)) continue;
        if (component[t] < 0) {  // flood the superlevel-set component containing this tip
            const int c = components++;
            group_of_component.push_back(-1);
            queue.assign(1, t);
            component[t] = c;
            for (std::size_t k = 0; k < queue.size(); ++k)
                for (std::uint32_t w : one_ring(m, queue[k]))
                    if (component[w] < 0 && open(w)) component[w] = c, queue.push_back(w);
        }
        int& g = group_of_component[static_cast<std::size_t>(component[t])];
        if (g < 0) g = *std::max_element(group.begin(), group.end()) + 1;
        group[i] = g;
    }
    return group;
}

MarginParams margin_operating_point() {
    // Chosen on 60 training scans (D73): ASSD 0.647 mm. sigma = 0 means raw kappa_min (the cervical
    // crease is narrow; smoothing blurs it); alpha in its converged regime (cost ~ integrated concavity).
    // tooth_quantile 1.0 (D77): height-band tooth seeds leaked tooth fronts across crease gaps onto flat
    // gingiva; only the highest arch vertex remains a height seed, so the cusp tips seed the teeth.
    // GraphCut (D78): mu 1000, beta 300, the centre of the plateau of the third sweep on the 60 sweep
    // scans (ASSD 0.323 vs 0.558 for Voronoi); validation 0.518 -> 0.344 mm.
    MarginParams p;
    p.valley_weight = 2560.0, p.curvature_scale = 0.0, p.gingiva_quantile = 0.15, p.tooth_quantile = 1.0, p.cusp_seeds = true;
    p.cut_smoothness = 1000.0, p.cut_crease = 300.0;
    // PerToothCut (D85): tau -1.25 (sweeps on the 60 sweep scans: 0.323 -> 0.272), no strip carving, expansion
    // radius 3 mm (identical metrics to the full expansion, 3x faster); validation 0.344 -> 0.283 mm.
    p.method = MarginParams::Method::PerToothCut;
    p.group_crease = -1.25;
    return p;
}

MarginResult detect_margin(const HalfEdgeMesh& m, const MarginParams& p) {
    const bool voronoi = p.method != MarginParams::Method::HeightPlane;  // Voronoi or GraphCut: same seeds
    const MarginInputs in = margin_inputs(m, voronoi && p.cusp_seeds);
    std::vector<double> valley;
    if (voronoi && p.valley_weight > 0.0) valley = valley_strength(build_dec(m), in.kmin, p.curvature_scale);
    MarginResult out;
    out.tooth = margin_labels(m, in, valley, p);
    out.margin = label_boundary_edges(m, out.tooth);
    return out;
}

}  // namespace dmw
