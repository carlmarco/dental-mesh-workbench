// Tooth-gingiva margin evaluation against Teeth3DS per-vertex labels (D71, D72, D73).
//   margin_eval sweep <scan-dir> [stride]   parameter sweep on TRAINING scans (every stride-th scan)
//   margin_eval test  <scan-dir>            fixed operating points on TEST scans (run once)
//   margin_eval validate <scan-dir> <stride> [min_remainder] [model.json]
//                                            VALIDATION: scans with index % stride >= min_remainder (default
//                                            1, i.e. not in the sweep); optional seed classifier (D74-D76)
// Each <name>.obj needs its <name>.json (labels: FDI per vertex, 0 = gingiva). Per-scan inputs are
// computed once and reused across configurations; test mode re-checks the library path (detect_margin).
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "core/io.h"
#include "core/learn.h"
#include "core/margin.h"
#include "core/topology.h"
#include "dataset.h"

using namespace dmw;
using Method = MarginParams::Method;

namespace {

struct Config {
    MarginParams p;
    std::string label() const {
        char b[160];
        if (p.method == Method::HeightPlane) std::snprintf(b, sizeof b, "plane cut q=%.2f", p.plane_quantile);
        else if (p.cusp_seed_quantile < 0.0) std::snprintf(b, sizeof b, "ORACLE: cusp seeds on true gingiva removed");
        else if (p.seed_model) std::snprintf(b, sizeof b, "Voronoi + seed classifier, drop P(tooth) < %.2f", p.seed_threshold);
        else std::snprintf(b, sizeof b, "Voronoi alpha=%.0f sigma=%.2f gq=%.2f tq=%.2f cusps=%d seed_q=%.1f", p.valley_weight,
                           p.curvature_scale, p.gingiva_quantile, p.tooth_quantile, int(p.cusp_seeds), p.cusp_seed_quantile);
        return b;
    }
};

MarginParams voronoi(double alpha, double sigma, double gq, double tq, bool cusps) {
    MarginParams p;
    p.method = Method::GeodesicVoronoi, p.valley_weight = alpha, p.curvature_scale = sigma;
    p.gingiva_quantile = gq, p.tooth_quantile = tq, p.cusp_seeds = cusps;
    return p;
}
MarginParams plane(double q) {
    MarginParams p;
    p.method = Method::HeightPlane, p.plane_quantile = q;
    return p;
}

std::vector<Config> sweep_configs() {
    // Widened until optima are interior (D73). Sweep 1 (61 configs): best at alpha 40 (largest) and
    // sigma 0.2 (smallest); gq 0.15 interior; cusp seeds essential. Sweep 2: best at alpha 320 (largest)
    // and sigma 0.1 (smallest), tq 0.9 > 0.8, ASSD 0.908 -> 0.716 mm. Sweep 3 (this): sigma down to 0
    // (raw kappa_min, the natural end of the range), alpha up to 2560.
    std::vector<Config> c;
    c.push_back({plane(0.5)});
    for (double a : {320.0, 640.0, 1280.0, 2560.0})
        for (double s : {0.0, 0.05, 0.1}) c.push_back({voronoi(a, s, 0.15, 0.9, true)});
    return c;
}

const LogisticModel* g_seed_model = nullptr;  // set from the command line (D76)

// Fixed operating points for the test set, from the TRAINING sweeps (60 scans of part 1, 2026-10-06).
// Edit only from training results, never from test results.
//   plane cut q 0.5 (ASSD 1.743); Voronoi without concavity weighting (alpha 0) to isolate its effect;
//   chosen: alpha 2560, sigma 0 (raw kappa_min; the natural end of the range), gq 0.15, tq 0.9, cusp
//   seeds (ASSD 0.647, median 0.461, F1@0.5 0.799, IoU 0.820; alpha gains converging: 1280 -> 2560 = -0.004).
// Second test evaluation (2026-10-06, disclosed in D76): adds the seed classifier at threshold 0.5,
// chosen on 120 validation scans for reliability (t = -4.1, worse on 2/120) over the best mean (0.85:
// -0.009 mm better than 0.5, statistically indistinguishable, worse on 34/120).
std::vector<Config> test_configs() {
    std::vector<Config> c{{plane(0.5)}, {voronoi(0.0, 0.0, 0.15, 0.9, true)}, {voronoi(2560.0, 0.0, 0.15, 0.9, true)}};
    if (g_seed_model) {
        MarginParams p = voronoi(2560.0, 0.0, 0.15, 0.9, true);
        p.seed_model = g_seed_model, p.seed_threshold = 0.5;
        c.push_back({p});
    }
    return c;
}

// Validation experiments (D75): seed filtering by a stricter height gate, plus an ORACLE that drops cusp
// seeds lying on true gingiva (uses labels: an upper bound on what seed filtering can gain, not a method).
constexpr double kOracle = -1.0;
std::vector<Config> validate_configs() {
    std::vector<Config> c;
    c.push_back({voronoi(2560.0, 0.0, 0.15, 0.9, true)});  // the operating point
    if (g_seed_model) {
        for (double t : {0.5, 0.8, 0.85, 0.9, 0.95}) {
            MarginParams p = voronoi(2560.0, 0.0, 0.15, 0.9, true);
            p.seed_model = g_seed_model, p.seed_threshold = t;
            c.push_back({p});
        }
    } else {
        for (double sq : {0.6, 0.7, 0.8, 0.9}) {
            MarginParams p = voronoi(2560.0, 0.0, 0.15, 0.9, true);
            p.cusp_seed_quantile = sq;
            c.push_back({p});
        }
    }
    MarginParams oracle = voronoi(2560.0, 0.0, 0.15, 0.9, true);
    oracle.cusp_seed_quantile = kOracle;
    c.push_back({oracle});
    return c;
}

struct Totals {
    std::vector<double> assd, hd95, f1_025, f1_050, iou;
    void add(const BoundaryMetrics& b, double i) {
        assd.push_back(b.assd), hd95.push_back(b.hd95), f1_025.push_back(b.f1_025), f1_050.push_back(b.f1_050), iou.push_back(i);
    }
};
double mean(const std::vector<double>& v) {
    double s = 0.0;
    for (double x : v) s += x;
    return v.empty() ? 0.0 : s / double(v.size());
}
double median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: margin_eval sweep|test <scan-dir> [stride]\n");
        return 2;
    }
    const std::string mode = argv[1];
    const std::size_t stride = argc > 3 ? std::stoul(argv[3]) : 1;
    const std::size_t min_remainder = argc > 4 ? std::stoul(argv[4]) : 1;
    LogisticModel seed_model;
    if (argc > 5) {
        std::string err;
        seed_model = parse_logistic_model(dataset::read_text(argv[5]), err);
        if (!err.empty()) {
            std::fprintf(stderr, "model: %s\n", err.c_str());
            return 2;
        }
    }
    g_seed_model = seed_model.weights.empty() ? nullptr : &seed_model;
    const auto objs = dataset::index_files({argv[2]}, ".obj");
    const auto labels = dataset::index_files({argv[2]}, ".json");
    const std::vector<Config> configs = mode == "sweep" ? sweep_configs() : mode == "validate" ? validate_configs() : test_configs();
    (void)seed_model;
    std::vector<Totals> totals(configs.size()), library(configs.size());

    std::size_t index = 0, used = 0, skipped = 0;
    double input_ms = 0.0;
    for (const auto& [stem, obj] : objs) {
        const std::size_t rem = index++ % stride;
        if (mode == "validate" ? rem < min_remainder : rem != 0) continue;
        const auto lab = labels.find(stem);
        if (lab == labels.end()) continue;
        const JsonResult j = parse_json(dataset::read_text(lab->second));
        const LoadResult r = parse_obj(dataset::read_text(obj));
        const JsonValue* arr = j.ok() ? j.value.find("labels") : nullptr;
        if (!r.ok() || !arr || arr->array.size() != r.mesh.positions.size()) { ++skipped; continue; }
        const AnalysisMesh analysis = manifold_analysis_mesh(r.mesh);
        if (!analysis.manifold) { ++skipped; continue; }
        const HalfEdgeMesh& m = analysis.halfedge;
        std::vector<std::uint8_t> truth(m.positions.size());
        for (std::size_t v = 0; v < truth.size(); ++v) truth[v] = arr->array[v].number != 0.0 ? 1 : 0;
        const auto gt_line = edge_midpoints(m, label_boundary_edges(m, truth));
        ++used;

        const auto t0 = std::chrono::steady_clock::now();
        const MarginInputs in = margin_inputs(m, true);
        const DecOperators ops = build_dec(m);
        input_ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        std::map<double, std::vector<double>> valley;  // by sigma
        for (std::size_t c = 0; c < configs.size(); ++c) {
            const MarginParams& p = configs[c].p;
            std::span<const double> val;
            if (p.method == Method::GeodesicVoronoi && p.valley_weight > 0.0) {
                auto it = valley.find(p.curvature_scale);
                if (it == valley.end()) it = valley.emplace(p.curvature_scale, valley_strength(ops, in.kmin, p.curvature_scale)).first;
                val = it->second;
            }
            std::vector<std::uint8_t> tooth;
            if (p.cusp_seed_quantile == kOracle) {  // oracle: drop seeds on true gingiva (diagnostic only)
                MarginInputs filtered = in;
                std::erase_if(filtered.cusp_tips, [&](std::uint32_t v) { return truth[v] == 0; });
                MarginParams pp = p;
                pp.cusp_seed_quantile = 0.0;
                tooth = margin_labels(m, filtered, val, pp);
            } else {
                tooth = margin_labels(m, in, val, p);
            }
            totals[c].add(compare_boundaries(edge_midpoints(m, label_boundary_edges(m, tooth)), gt_line),
                          region_iou(tooth, truth, ops.star0));
            if (mode == "test") {  // the library path, end to end
                const MarginResult lib = detect_margin(m, p);
                library[c].add(compare_boundaries(edge_midpoints(m, lib.margin), gt_line), region_iou(lib.tooth, truth, ops.star0));
            }
        }
        std::fprintf(stderr, "\r%zu scans", used);
    }
    std::fprintf(stderr, "\n");
    std::printf("scans evaluated: %zu (skipped: %zu), mean input time %.0f ms/scan\n\n", used, skipped, input_ms / double(std::max<std::size_t>(used, 1)));
    std::vector<std::size_t> order(configs.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    if (mode == "sweep") std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return mean(totals[a].assd) < mean(totals[b].assd); });
    std::printf("| configuration | ASSD mean (mm) | ASSD median | HD95 mean (mm) | boundary F1@0.25 | F1@0.5 | tooth IoU |\n|---|---:|---:|---:|---:|---:|---:|\n");
    int shown_plane = 0, shown_vor = 0;
    for (std::size_t i : order) {
        const bool is_plane = configs[i].p.method == Method::HeightPlane;
        if (mode == "sweep" && (is_plane ? shown_plane++ >= 3 : shown_vor++ >= 8)) continue;
        const Totals& t = totals[i];
        std::printf("| %s | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |\n", configs[i].label().c_str(), mean(t.assd), median(t.assd),
                    mean(t.hd95), mean(t.f1_025), mean(t.f1_050), mean(t.iou));
    }
    if (mode == "sweep") {  // the cusp-seed ablation, wherever it ranked
        for (std::size_t i = 0; i < configs.size(); ++i)
            if (!configs[i].p.cusp_seeds && configs[i].p.method == Method::GeodesicVoronoi)
                std::printf("| (ablation) %s | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f |\n", configs[i].label().c_str(), mean(totals[i].assd),
                            median(totals[i].assd), mean(totals[i].hd95), mean(totals[i].f1_025), mean(totals[i].f1_050), mean(totals[i].iou));
    }
    if (mode == "validate" || (mode == "test" && g_seed_model)) {
        // Paired comparison against the operating point on the same scans: is a gain real?
        const std::size_t base_index = mode == "validate" ? 0 : 2;
        std::printf("\npaired vs operating point (per-scan ASSD difference; negative = better):\n");
        const auto& base = totals[base_index].assd;
        for (std::size_t c = base_index + 1; c < configs.size(); ++c) {
            const auto& t = totals[c].assd;
            double sum = 0.0, sq = 0.0;
            std::size_t wins = 0, losses = 0;
            for (std::size_t i = 0; i < base.size(); ++i) {
                const double d = t[i] - base[i];
                sum += d, sq += d * d;
                wins += d < -1e-9, losses += d > 1e-9;
            }
            const double n = double(base.size()), mean_d = sum / n;
            const double se = std::sqrt(std::max(sq / n - mean_d * mean_d, 0.0) / (n - 1.0));
            std::printf("  %-55s mean %+.4f mm, SE %.4f (t = %+.1f), better on %zu, worse on %zu of %zu scans\n",
                        configs[c].label().c_str(), mean_d, se, se > 0 ? mean_d / se : 0.0, wins, losses, base.size());
        }
    }
    if (mode == "test") {
        std::printf("\nlibrary path check (detect_margin, mean ASSD):");
        for (const auto& t : library) std::printf(" %.3f", mean(t.assd));
        std::printf("\n");
    }
    return 0;
}
