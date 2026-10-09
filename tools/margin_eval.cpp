// Tooth-gingiva margin evaluation against Teeth3DS per-vertex labels (D71, D72, D73).
//   margin_eval sweep <scan-dir> [stride]   parameter sweep on TRAINING scans (every stride-th scan)
//   margin_eval test  <scan-dir>            fixed operating points on TEST scans (run once)
//   margin_eval score <scan-dir> <pred-dir> [label]
//                                            EXTERNAL predictions (D81): <pred-dir>/<name>.json in the Teeth3DS
//                                            label format (e.g. ToothGroupNetwork output), any label != 0 = tooth;
//                                            scored with the same metrics, paired against the public graph cut
//   margin_eval validate <scan-dir> <stride> [min_remainder] [model.json]
//                                            VALIDATION: scans with index % stride >= min_remainder (default
//                                            1, i.e. not in the sweep); optional seed classifier (D74-D76)
// Each <name>.obj needs its <name>.json (labels: FDI per vertex, 0 = gingiva). Per-scan inputs are
// computed once and reused across configurations; test mode re-checks the library path (detect_margin).
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <limits>
#include <map>
#include <set>
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

// Sentinel values of cusp_seed_quantile that select label-using ORACLE diagnostics (D75, D77).
constexpr double kOracle = -1.0, kHeightOracle = -2.0, kBothOracles = -3.0;
constexpr double kExternal = -4.0;  // score mode: labels read from a prediction file, not computed (D81)
std::string g_external_label = "external predictions";

struct Config {
    MarginParams p;
    std::string label() const {
        char b[160];
        if (p.cusp_seed_quantile == kExternal) std::snprintf(b, sizeof b, "%s", g_external_label.c_str());
        else if (p.method == Method::HeightPlane) std::snprintf(b, sizeof b, "plane cut q=%.2f", p.plane_quantile);
        else if (p.method == Method::GraphCut)
            std::snprintf(b, sizeof b, "GraphCut mu=%.2f beta=%.0f alpha=%.0f tq=%.2f islands<%.0f%s", p.cut_smoothness, p.cut_crease,
                          p.valley_weight, p.tooth_quantile, p.min_tooth_region, p.seed_model ? " + classifier" : "");
        else if (p.method == Method::PerToothCut)
            std::snprintf(b, sizeof b, "PerToothCut tau=%s strip<%s R=%.0f label_r=%.0f%s%s", p.group_crease >= 1e8 ? "none" : std::to_string(p.group_crease).substr(0, 5).c_str(),
                          p.strip_kmin < -1e8 ? "off" : std::to_string(p.strip_kmin).substr(0, 4).c_str(), p.expansion_radius,
                          p.label_radius, p.star_prior ? (p.star_plain ? " star(plain)" : " star") : "",
                          p.vertex_model ? (" + unary w=" + std::to_string(p.vertex_weight).substr(0, 4)).c_str() : (p.seed_model ? " + classifier" : ""));
        else if (p.cusp_seed_quantile == kHeightOracle) std::snprintf(b, sizeof b, "ORACLE: height seeds on true gingiva removed (+clf)");
        else if (p.cusp_seed_quantile == kBothOracles) std::snprintf(b, sizeof b, "ORACLE: all tooth seeds on true gingiva removed");
        else if (p.cusp_seed_quantile < 0.0) std::snprintf(b, sizeof b, "ORACLE: cusp seeds on true gingiva removed");
        else if (p.seed_model) std::snprintf(b, sizeof b, "Voronoi + seed classifier %.2f, tq=%.2f", p.seed_threshold, p.tooth_quantile);
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
MarginTimings g_timings;  // D86 profile experiment
const LogisticModel* g_vertex_model = nullptr;  // D88: from $DMW_VERTEX_MODEL (learned per-vertex data term)

// Fixed operating points for the test set, from the TRAINING sweeps (60 scans of part 1, 2026-10-06).
// Edit only from training results, never from test results.
//   plane cut q 0.5 (ASSD 1.743); Voronoi without concavity weighting (alpha 0) to isolate its effect;
//   chosen: alpha 2560, sigma 0 (raw kappa_min; the natural end of the range), gq 0.15, tq 0.9, cusp
//   seeds (ASSD 0.647, median 0.461, F1@0.5 0.799, IoU 0.820; alpha gains converging: 1280 -> 2560 = -0.004).
// Test evaluations (disclosed): #1 (D73) plane / alpha-0 / chosen Voronoi; #2 (D76) + seed classifier 0.5;
// #3 (D77) no height seeds (tooth_quantile 1.0); #4 (D78, 2026-10-06) graph-cut labelling, mu 1000, beta 300,
// tuned on the 60 sweep scans (3 sweeps) and validated on 120 scans before this run. Row 0 is the baseline.
MarginParams graph_cut(bool clf) {
    MarginParams p = voronoi(2560.0, 0.0, 0.15, 1.0, true);
    p.method = Method::GraphCut, p.cut_smoothness = 1000.0, p.cut_crease = 300.0;
    if (clf) p.seed_model = g_seed_model, p.seed_threshold = 0.5;
    return p;
}
std::vector<Config> score_configs() {
    MarginParams ext;
    ext.cusp_seed_quantile = kExternal;
    return {{graph_cut(false)}, {ext}};
}
std::vector<Config> test_configs() {
    // #5 (D85, 2026-10-08): per-tooth cut, tau -1.25, chosen on sweep scans and validated before the run. On part 6
    // this is its second look (the first, D81, evaluated the binary graph cut). Row 0 = the D78/D81 public point.
    MarginParams pt = graph_cut(false), ptc;
    pt.method = Method::PerToothCut, pt.group_crease = -1.25;
    ptc = pt;
    std::vector<Config> c{{graph_cut(false)}, {pt}};
    if (g_seed_model) {
        ptc.seed_model = g_seed_model, ptc.seed_threshold = 0.5;
        c.push_back({ptc});
    }
    return c;
}

// Validation experiments (D75): seed filtering by a stricter height gate, plus an ORACLE that drops cusp
// seeds lying on true gingiva (uses labels: an upper bound on what seed filtering can gain, not a method).
std::string g_experiment = "thresholds";      // validate mode: which experiment (argv[6])
std::vector<Config> validate_configs() {
    std::vector<Config> c;
    if (g_experiment == "graphcut") {
        // D78: tune the graph cut on TRAINING scans (the 60 sweep scans), paired vs the public point.
        // Sweep 1: mu {0.3, 1, 3, 10} x beta {0, 100, 1000}: best mu 10 (grid edge), beta 100 (interior);
        // beta 0 is far worse at every mu. Sweep 2: mu {10..300} x beta {30, 100, 300}: best mu 300 (edge),
        // beta 100: ASSD 0.315. Sweep 3 (this): mu up to 1e6, where the unary term is negligible (pure
        // crease-weighted minimum cut between the seeds).
        c.push_back({voronoi(2560.0, 0.0, 0.15, 1.0, true)});
        for (double mu : {300.0, 1000.0, 3000.0, 10000.0, 1e6})
            for (double beta : {30.0, 100.0, 300.0}) {
                MarginParams p = voronoi(2560.0, 0.0, 0.15, 1.0, true);
                p.method = Method::GraphCut, p.cut_smoothness = mu, p.cut_crease = beta;
                c.push_back({p});
            }
        return c;
    }
    if (g_experiment == "errors") {
        c.push_back({graph_cut(false)});  // D80: where does the per-vertex error of the public method sit?
        return c;
    }
    if (g_experiment == "unary") {  // D88: learned data term weight w on the per-tooth cut (sweep scans)
        MarginParams p = graph_cut(false);
        p.method = Method::PerToothCut, p.group_crease = -1.25;
        c.push_back({p});
        for (double w : {0.25, 0.5, 1.0, 2.0, 4.0}) {
            MarginParams q = p;
            q.vertex_model = g_vertex_model, q.vertex_weight = w;
            c.push_back({q});
        }
        return c;
    }
    if (g_experiment == "star") {  // D87: geodesic star-convexity prior on the per-tooth cut (sweep scans)
        MarginParams p = graph_cut(false);
        p.method = Method::PerToothCut, p.group_crease = -1.25;
        MarginParams st = p;
        st.star_prior = true;
        // 12-scan slice: weighted trees -0.001 mm (t -1.0); straight-ray trees +0.023 (2 scans worse by > 0.1 mm).
        c.push_back({p});
        c.push_back({st});
        return c;
    }
    if (g_experiment == "speed_check") {  // D86: the faster per-tooth cut must reproduce sweep 3 (tau -1.25: 0.272)
        MarginParams p = graph_cut(false);
        p.method = Method::PerToothCut, p.group_crease = -1.25;
        c.push_back({graph_cut(false)});
        c.push_back({p});
        return c;
    }
    if (g_experiment == "cutoff") {  // D86: per-label Dijkstra cutoff; metrics must not move
        // A weighted-distance cutoff (5..50) was catastrophic: crease weighting makes distances inside one crown
        // reach thousands. Now: a straight-line radius around each label's seed tips.
        for (double rad : {0.0, 20.0, 15.0, 12.0, 9.0}) {
            MarginParams p = graph_cut(false);
            p.method = Method::PerToothCut, p.group_crease = -1.25, p.label_radius = rad;
            c.push_back({p});
        }
        return c;
    }
    if (g_experiment == "profile") {  // D86: per-stage timing of the public per-tooth cut
        MarginParams p = graph_cut(false);
        p.method = Method::PerToothCut, p.group_crease = -1.25;
        p.timings = &g_timings;
        c.push_back({p});
        return c;
    }
    if (g_experiment == "pertooth_validate" || g_experiment == "errors_pertooth") {  // D85: chosen tau -1.25, R 3
        MarginParams p = graph_cut(false);
        p.method = Method::PerToothCut, p.group_crease = -1.25;
        if (g_experiment == "pertooth_validate") {
            c.push_back({graph_cut(false)});
            c.push_back({p});
            MarginParams pc = p;
            pc.seed_model = g_seed_model, pc.seed_threshold = 0.5;
            c.push_back({pc});
        } else {
            c.push_back({p});
        }
        return c;
    }
    if (g_experiment == "pertooth") {  // D85: per-tooth cut vs the public binary cut (sweep scans)
        c.push_back({graph_cut(false)});
        auto pt = [](double tau, double strip) {
            MarginParams p = graph_cut(false);
            p.method = Method::PerToothCut, p.group_crease = tau, p.strip_kmin = strip;
            return p;
        };
        // Sweep 1: tau -1.5 without carving 0.285 (t -3.0, worse 3/60); every carving variant (strip < 0 or -2, tau
        // -1.25/-1.5/-1.75) was worse than the binary cut (+0.040..+0.080). Sweep 2 (this): tau, no carving.
        // Sweep 2: tau {-1.25..-2}, no carving: best -1.25 (0.272, t -3.2), the grid edge. Sweep 3 (this): toward
        // more splitting, up to one label per tip (1e9); plus -1.25 without the expansion radius (equivalence check).
        for (double tau : {-0.75, -1.0, -1.25, 1e9}) c.push_back({pt(tau, -1e9)});
        MarginParams full = pt(-1.25, -1e9);
        full.expansion_radius = 0.0;
        c.push_back({full});
        return c;
    }
    if (g_experiment == "groups") {  // D85: cusp grouping on the public cut's tooth region (sweep scans)
        c.push_back({graph_cut(false)});
        return c;
    }
    if (g_experiment == "errors_voronoi") {  // D82: the same decomposition for first-arrival Voronoi
        c.push_back({voronoi(2560.0, 0.0, 0.15, 1.0, true)});
        return c;
    }
    if (g_experiment == "islands_validate") {
        MarginParams p = graph_cut(false);
        p.min_tooth_region = 20.0;
        c.push_back({graph_cut(false)});
        c.push_back({p});
        return c;
    }
    if (g_experiment == "islands") {
        // D80: island removal (as in ToothGroupNetwork's post-processing); tuned on the 60 sweep scans.
        // (Seed discs of radius 0.5-2 mm were tried first and hurt: +0.042..+0.092 mm.)
        c.push_back({graph_cut(false)});
        for (double a : {20.0, 40.0, 80.0}) {  // sweep 1: 2, 5, 10, 20 (best 20, edge); sweep 2: 20 best, 80 deletes teeth
            MarginParams p = graph_cut(false);
            p.min_tooth_region = a;
            c.push_back({p});
        }
        return c;
    }
    if (g_experiment == "graphcut_validate") {
        // D78: the chosen cut (mu 1000, beta 300; centre of the sweep-3 plateau) vs Voronoi, with and
        // without the seed classifier, on validation scans.
        auto cut = [](bool clf) {
            MarginParams p = voronoi(2560.0, 0.0, 0.15, 1.0, true);
            p.method = Method::GraphCut, p.cut_smoothness = 1000.0, p.cut_crease = 300.0;
            if (clf) p.seed_model = g_seed_model, p.seed_threshold = 0.5;
            return p;
        };
        MarginParams vclf = voronoi(2560.0, 0.0, 0.15, 1.0, true);
        vclf.seed_model = g_seed_model, vclf.seed_threshold = 0.5;
        c.push_back({voronoi(2560.0, 0.0, 0.15, 1.0, true)});
        c.push_back({cut(false)});
        c.push_back({vclf});
        c.push_back({cut(true)});
        return c;
    }
    if (g_experiment == "confirm") {
        // D77: confirm the height-seed effect on independent scans (the 60 sweep scans).
        MarginParams base = voronoi(2560.0, 0.0, 0.15, 0.9, true);
        base.seed_model = g_seed_model, base.seed_threshold = 0.5;
        c.push_back({base});
        for (double tq : {0.98, 1.0}) {
            MarginParams p = base;
            p.tooth_quantile = tq;
            c.push_back({p});
        }
        MarginParams no_clf = voronoi(2560.0, 0.0, 0.15, 1.0, true);  // does the classifier still matter?
        c.push_back({no_clf});
        return c;
    }
    if (g_experiment == "height") {
        // D77: height-seed experiments on top of the classifier (threshold 0.5) operating point.
        auto with_clf = [](double tq) {
            MarginParams p = voronoi(2560.0, 0.0, 0.15, tq, true);
            p.seed_model = g_seed_model, p.seed_threshold = 0.5;
            return p;
        };
        c.push_back({with_clf(0.9)});  // current method
        for (double tq : {0.95, 0.98, 1.0}) c.push_back({with_clf(tq)});
        MarginParams ho = with_clf(0.9), both = with_clf(0.9);
        ho.cusp_seed_quantile = kHeightOracle;
        both.cusp_seed_quantile = kBothOracles;
        c.push_back({ho});
        c.push_back({both});
        return c;
    }
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
    double label_ms = 0.0;  // labelling step only (inputs are shared across configurations)
    // Per-vertex (unweighted) binary metrics, the convention of the deep-learning papers (comparison only).
    std::vector<double> vacc, viou_tooth, viou_gingiva;
    void add_vertex(std::span<const std::uint8_t> pred, std::span<const std::uint8_t> truth) {
        std::size_t ok = 0, tt = 0, tu = 0, gi = 0, gu = 0;
        for (std::size_t v = 0; v < pred.size(); ++v) {
            ok += pred[v] == truth[v];
            tt += pred[v] && truth[v], tu += pred[v] || truth[v];
            gi += !pred[v] && !truth[v], gu += !pred[v] || !truth[v];
        }
        vacc.push_back(double(ok) / double(pred.size()));
        viou_tooth.push_back(tu ? double(tt) / double(tu) : 1.0);
        viou_gingiva.push_back(gu ? double(gi) / double(gu) : 1.0);
    }
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

// D80 error breakdown (area-weighted, summed over scans): false gingiva (tooth vertices labelled gingiva)
// split into whole missed teeth (< 50% of the tooth labelled tooth) vs partially covered teeth; false
// tooth; plus per-tooth-type miss rates and whether a missed tooth had a seed.
struct ErrorStats {
    double false_gingiva_missed = 0, false_gingiva_partial = 0, false_tooth = 0, false_tooth_near = 0, total_tooth = 0, total_area = 0;
    std::map<int, std::pair<int, int>> by_type;  // FDI unit digit -> (missed, total)
    int missed_with_seed = 0, missed_total = 0;
    // D82 F1@0.5 decomposition, pooled over scans (boundary samples = edge midpoints, as in compare_boundaries).
    // Recall misses (true samples > 0.5 mm from the predicted line), by the true tooth beside them:
    std::size_t gt_samples = 0, r_missed_tooth = 0, r_near = 0, r_far = 0;  // near: 0.5-1 mm; far: > 1 mm
    // Far recall misses: is the true-gingiva vertex beside them labelled tooth (over-extension), and does it lie
    // within 1.5 mm of two different teeth (interdental papilla)?
    std::size_t r_far_overext = 0, r_far_interdental = 0, r_far_overext_interdental = 0;
    // D83: crease strength along the TRUE boundary. kappa_min at each true boundary edge (mean of endpoints) and
    // the strongest crease within the one-rings of its endpoints, split interdental / other and matched / missed
    // by the cut; plus all true-gingiva arch vertices as a flat baseline.
    std::vector<double> k_edge[2][2], k_ring[2][2];  // [interdental][missed]
    std::vector<double> k_gingiva;
    std::vector<double> kp_ring[2][2];  // PREDICTED boundary edges: [interdental][> 1 mm from the true line]
    // Precision misses (predicted samples > 0.5 mm from the true line), by the predicted-tooth vertex beside them:
    std::size_t pr_samples = 0, p_fake = 0, p_near = 0, p_far = 0;  // fake: that vertex is gingiva > 0.5 mm from any tooth
    // D84: per-scan worst cases for visual inspection: missed interdental true samples, and the densest one.
    struct Case {
        std::size_t missed;
        std::string obj, labels;
        Vec3 focus;
    };
    std::vector<Case> cases;
    // D84: width of the labelled gingiva strip at missed interdental samples: distance from the sample to the
    // nearest vertex of a tooth OTHER than the one beside it (the neighbouring crown across the gap).
    std::vector<double> strip_width;
};
std::string g_scan_obj, g_scan_labels;  // the scan being evaluated (for the inspection list)
ErrorStats g_errors;
void error_breakdown(const HalfEdgeMesh& m, const MarginInputs& in, std::span<const std::uint8_t> tooth,
                     std::span<const std::uint8_t> truth, const std::vector<int>& fdi, std::span<const double> area) {
    std::map<int, double> covered, total;
    for (std::size_t v = 0; v < tooth.size(); ++v) {
        g_errors.total_area += area[v];
        if (fdi[v] != 0) total[fdi[v]] += area[v], covered[fdi[v]] += tooth[v] ? area[v] : 0.0;
        if (tooth[v] && !truth[v]) g_errors.false_tooth += area[v];
    }
    // False tooth within 0.5 mm of a true tooth vertex = boundary offset; farther = leak / fake region.
    std::vector<Vec3> ft, tt;
    std::vector<double> ft_area;
    for (std::size_t v = 0; v < tooth.size(); ++v) {
        if (truth[v]) tt.push_back(m.positions[v]);
        if (tooth[v] && !truth[v]) ft.push_back(m.positions[v]), ft_area.push_back(area[v]);
    }
    const auto dn = nearest_point_distances(ft, tt, 0.5);
    for (std::size_t k = 0; k < ft.size(); ++k) g_errors.false_tooth_near += dn[k] <= 0.5 ? ft_area[k] : 0.0;
    std::map<int, bool> seeded;
    for (std::uint32_t v : in.cusp_tips)
        if (fdi[v] != 0) seeded[fdi[v]] = true;
    for (const auto& [t, a] : total) {
        const bool missed = covered[t] < 0.5 * a;
        g_errors.total_tooth += a;
        (missed ? g_errors.false_gingiva_missed : g_errors.false_gingiva_partial) += a - covered[t];
        auto& bt = g_errors.by_type[t % 10];
        bt.second++, bt.first += missed;
        if (missed) g_errors.missed_total++, g_errors.missed_with_seed += seeded[t];
    }

    // --- F1@0.5 decomposition (D82).
    const auto gt_edges = label_boundary_edges(m, truth), pr_edges = label_boundary_edges(m, tooth);
    const auto gt_pts = edge_midpoints(m, gt_edges), pr_pts = edge_midpoints(m, pr_edges);
    const auto gt_to_pr = nearest_point_distances(gt_pts, pr_pts, 0.5), pr_to_gt = nearest_point_distances(pr_pts, gt_pts, 0.5);
    g_errors.gt_samples += gt_pts.size();
    g_errors.pr_samples += pr_pts.size();
    for (std::size_t k = 0; k < gt_pts.size(); ++k) {
        if (gt_to_pr[k] <= 0.5) continue;
        const std::uint32_t tv = truth[gt_edges[k].v0] ? gt_edges[k].v0 : gt_edges[k].v1;  // the true-tooth side
        const int t = fdi[tv];
        const bool missed = t != 0 && covered[t] < 0.5 * total[t];
        if (missed) {
            ++g_errors.r_missed_tooth;
        } else if (gt_to_pr[k] <= 1.0) {
            ++g_errors.r_near;
        } else {
            ++g_errors.r_far;
            const std::uint32_t gv = tv == gt_edges[k].v0 ? gt_edges[k].v1 : gt_edges[k].v0;  // true-gingiva side
            const bool over = tooth[gv] != 0;
            // Distinct teeth within 1.5 mm of the gingiva vertex (brute force: few far misses per scan).
            int first = 0;
            bool two = false;
            for (std::size_t v = 0; v < fdi.size() && !two; ++v) {
                if (fdi[v] == 0 || fdi[v] == first) continue;
                const Vec3 d{m.positions[v].x - m.positions[gv].x, m.positions[v].y - m.positions[gv].y, m.positions[v].z - m.positions[gv].z};
                if (d.x * d.x + d.y * d.y + d.z * d.z > 2.25) continue;
                if (first == 0) first = fdi[v];
                else two = true;
            }
            g_errors.r_far_overext += over, g_errors.r_far_interdental += two, g_errors.r_far_overext_interdental += over && two;
        }
    }
    // --- D83 crease measurement along the true boundary.
    {
        // Distinct teeth within 1.5 mm of each true-gingiva-side vertex: fixed-radius lookup on a 1.5 mm grid
        // (27 cells), O(n). (A per-tooth nearest-point query expands its search across the arch: too slow.)
        constexpr double kR = 1.5;
        using Key = std::array<std::int64_t, 3>;
        auto key = [](const Vec3& p) {
            return Key{static_cast<std::int64_t>(std::floor(p.x / kR)), static_cast<std::int64_t>(std::floor(p.y / kR)),
                       static_cast<std::int64_t>(std::floor(p.z / kR))};
        };
        std::map<Key, std::vector<std::uint32_t>> grid;
        for (std::uint32_t v = 0; v < fdi.size(); ++v)
            if (fdi[v] != 0) grid[key(m.positions[v])].push_back(v);
        auto teeth_near = [&](const Vec3& g) {  // 0, 1 or 2 (= at least two distinct teeth within 1.5 mm)
            const Key c = key(g);
            int first = 0, n = 0;
            for (std::int64_t dx = -1; dx <= 1 && n < 2; ++dx)
                for (std::int64_t dy = -1; dy <= 1 && n < 2; ++dy)
                    for (std::int64_t dz = -1; dz <= 1 && n < 2; ++dz) {
                        const auto it = grid.find({c[0] + dx, c[1] + dy, c[2] + dz});
                        if (it == grid.end()) continue;
                        for (std::uint32_t v : it->second) {
                            const Vec3 d{m.positions[v].x - g.x, m.positions[v].y - g.y, m.positions[v].z - g.z};
                            if (d.x * d.x + d.y * d.y + d.z * d.z > kR * kR || fdi[v] == first) continue;
                            if (first == 0) first = fdi[v], n = 1;
                            else { n = 2; break; }
                        }
                    }
            return n;
        };
        auto ring_min = [&](std::uint32_t a, std::uint32_t b) {
            double r = std::min(in.kmin[a], in.kmin[b]);
            for (std::uint32_t w : one_ring(m, a)) r = std::min(r, in.kmin[w]);
            for (std::uint32_t w : one_ring(m, b)) r = std::min(r, in.kmin[w]);
            return r;
        };
        std::vector<int> near_teeth(gt_edges.size(), 0);
        for (std::size_t k = 0; k < gt_edges.size(); ++k)
            near_teeth[k] = teeth_near(m.positions[truth[gt_edges[k].v0] ? gt_edges[k].v1 : gt_edges[k].v0]);
        for (std::size_t k = 0; k < pr_edges.size(); ++k) {
            const int inter = teeth_near(pr_pts[k]) >= 2, far = pr_to_gt[k] > 1.0;
            g_errors.kp_ring[inter][far].push_back(ring_min(pr_edges[k].v0, pr_edges[k].v1));
        }
        std::vector<Vec3> missed_inter;
        for (std::size_t k = 0; k < gt_edges.size(); ++k) {
            if (near_teeth[k] < 2 || gt_to_pr[k] <= 0.5) continue;
            missed_inter.push_back(gt_pts[k]);
            const int own = fdi[truth[gt_edges[k].v0] ? gt_edges[k].v0 : gt_edges[k].v1];
            const Key c = key(gt_pts[k]);
            double best = std::numeric_limits<double>::infinity();
            for (std::int64_t dx = -1; dx <= 1; ++dx)
                for (std::int64_t dy = -1; dy <= 1; ++dy)
                    for (std::int64_t dz = -1; dz <= 1; ++dz) {
                        const auto it = grid.find({c[0] + dx, c[1] + dy, c[2] + dz});
                        if (it == grid.end()) continue;
                        for (std::uint32_t v : it->second) {
                            if (fdi[v] == own) continue;
                            const Vec3 d{m.positions[v].x - gt_pts[k].x, m.positions[v].y - gt_pts[k].y, m.positions[v].z - gt_pts[k].z};
                            best = std::min(best, std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z));
                        }
                    }
            g_errors.strip_width.push_back(best);  // <= 1.5 mm by construction (interdental test)
        }
        if (!missed_inter.empty()) {
            // Focus = the missed interdental sample with the most others within 2 mm (brute force; small sets).
            std::size_t best = 0, best_n = 0;
            for (std::size_t i = 0; i < missed_inter.size(); i += 4) {
                std::size_t n = 0;
                for (const Vec3& q : missed_inter) {
                    const Vec3 d{q.x - missed_inter[i].x, q.y - missed_inter[i].y, q.z - missed_inter[i].z};
                    n += d.x * d.x + d.y * d.y + d.z * d.z < 4.0;
                }
                if (n > best_n) best = i, best_n = n;
            }
            g_errors.cases.push_back({missed_inter.size(), g_scan_obj, g_scan_labels, missed_inter[best]});
        }
        for (std::size_t k = 0; k < gt_edges.size(); ++k) {
            const std::uint32_t a = gt_edges[k].v0, b = gt_edges[k].v1;
            const double ring = ring_min(a, b);
            const int inter = near_teeth[k] >= 2, miss = gt_to_pr[k] > 0.5;
            g_errors.k_edge[inter][miss].push_back(0.5 * (in.kmin[a] + in.kmin[b]));
            g_errors.k_ring[inter][miss].push_back(ring);
        }
        for (std::size_t v = 0; v < truth.size(); v += 7)  // subsample: baseline only
            if (in.arch[v] && !truth[v]) g_errors.k_gingiva.push_back(in.kmin[v]);
    }
    // Distance from each predicted-tooth-side vertex to the nearest true tooth vertex (fake-region test).
    std::vector<Vec3> side;
    for (const Edge& e : pr_edges) side.push_back(m.positions[tooth[e.v0] ? e.v0 : e.v1]);
    const auto side_to_tooth = nearest_point_distances(side, tt, 0.5);
    for (std::size_t k = 0; k < pr_pts.size(); ++k) {
        if (pr_to_gt[k] <= 0.5) continue;
        const std::uint32_t pv = tooth[pr_edges[k].v0] ? pr_edges[k].v0 : pr_edges[k].v1;
        if (!truth[pv] && side_to_tooth[k] > 0.5) ++g_errors.p_fake;
        else if (pr_to_gt[k] <= 1.0) ++g_errors.p_near;
        else ++g_errors.p_far;
    }
}

// D85: how well does group_cusps recover teeth? Tips on true teeth only (FDI at the tip vertex); pairs of
// tips: same tooth & same group = kept together; different teeth & same group = merge; same tooth & different
// groups = split. A tooth is "clean" if all its tips share one group that holds no other tooth's tips.
const std::vector<double> kGroupTaus{-0.75, -1.0, -1.25, -1.5, -1.75, -2.0};  // sweep 1: -0.5..-6, best -1 (clean 82%)
struct GroupStats {
    std::size_t together = 0, merged = 0, split = 0, teeth = 0, clean = 0, tips_on_gingiva = 0, tips = 0, ungrouped = 0;
};
std::vector<GroupStats> g_groups(kGroupTaus.size());
void group_eval(const HalfEdgeMesh& m, const MarginInputs& in, std::span<const std::uint8_t> tooth, const std::vector<int>& fdi) {
    for (std::size_t ti = 0; ti < kGroupTaus.size(); ++ti) {
        const auto g = group_cusps(m, in.cusp_tips, tooth, in.kmin, kGroupTaus[ti]);
        GroupStats& st = g_groups[ti];
        std::map<int, std::vector<int>> groups_of_tooth;  // FDI -> groups of its tips
        std::map<int, std::set<int>> teeth_of_group;
        for (std::size_t i = 0; i < g.size(); ++i) {
            const int t = fdi[in.cusp_tips[i]];
            ++st.tips;
            if (t == 0) { ++st.tips_on_gingiva; continue; }
            if (g[i] < 0) { ++st.ungrouped; continue; }
            groups_of_tooth[t].push_back(g[i]);
            teeth_of_group[g[i]].insert(t);
        }
        std::vector<std::pair<int, int>> tg;  // (tooth, group) per grouped tip on a tooth
        for (const auto& [t, gs] : groups_of_tooth)
            for (int x : gs) tg.push_back({t, x});
        for (std::size_t i = 0; i < tg.size(); ++i)
            for (std::size_t j = i + 1; j < tg.size(); ++j) {
                const bool same_tooth = tg[i].first == tg[j].first, same_group = tg[i].second == tg[j].second;
                st.together += same_tooth && same_group, st.merged += !same_tooth && same_group, st.split += same_tooth && !same_group;
            }
        for (const auto& [t, gs] : groups_of_tooth) {
            ++st.teeth;
            const bool one = std::all_of(gs.begin(), gs.end(), [&](int x) { return x == gs[0]; });
            st.clean += one && teeth_of_group[gs[0]].size() == 1;
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: margin_eval sweep|test|validate|score <scan-dir> ... (see the header comment)\n");
        return 2;
    }
    const std::string mode = argv[1];
    const bool score = mode == "score";
    if (score && argc < 4) {
        std::fprintf(stderr, "usage: margin_eval score <scan-dir> <pred-dir> [label]\n");
        return 2;
    }
    if (score && argc > 4) g_external_label = argv[4];
    const std::size_t stride = !score && argc > 3 ? std::stoul(argv[3]) : 1;
    const std::size_t min_remainder = !score && argc > 4 ? std::stoul(argv[4]) : 1;
    LogisticModel seed_model;
    if (!score && argc > 5) {
        std::string err;
        seed_model = parse_logistic_model(dataset::read_text(argv[5]), err);
        if (!err.empty()) {
            std::fprintf(stderr, "model: %s\n", err.c_str());
            return 2;
        }
    }
    g_seed_model = seed_model.weights.empty() ? nullptr : &seed_model;
    LogisticModel vertex_model;
    if (const char* vm = std::getenv("DMW_VERTEX_MODEL")) {
        std::string err;
        vertex_model = parse_logistic_model(dataset::read_text(vm), err);
        if (!err.empty()) {
            std::fprintf(stderr, "vertex model: %s\n", err.c_str());
            return 2;
        }
        g_vertex_model = &vertex_model;
    }
    if (!score && argc > 6) g_experiment = argv[6];
    const std::size_t max_remainder = !score && argc > 7 ? std::stoul(argv[7]) : stride - 1;
    const auto objs = dataset::index_files({argv[2]}, ".obj");
    const auto labels = dataset::index_files({argv[2]}, ".json");
    const auto predictions = score ? dataset::index_files({argv[3]}, ".json") : decltype(labels){};
    const std::vector<Config> configs = mode == "sweep"      ? sweep_configs()
                                        : mode == "validate" ? validate_configs()
                                        : score              ? score_configs()
                                                             : test_configs();
    std::size_t no_prediction = 0;
    (void)seed_model;
    std::vector<Totals> totals(configs.size()), library(configs.size());

    std::size_t index = 0, used = 0, skipped = 0;
    double input_ms = 0.0;
    for (const auto& [stem, obj] : objs) {
        const std::size_t rem = index++ % stride;
        if (mode == "validate" ? (rem < min_remainder || rem > max_remainder) : rem != 0) continue;
        const auto lab = labels.find(stem);
        if (lab == labels.end()) continue;
        std::vector<std::uint8_t> external;  // score mode: the prediction, binarized (any tooth label = tooth)
        if (score) {
            const auto pr = predictions.find(stem);
            if (pr == predictions.end()) { ++no_prediction; continue; }
            const JsonResult pj = parse_json(dataset::read_text(pr->second));
            const JsonValue* parr = pj.ok() ? pj.value.find("labels") : nullptr;
            if (!parr) { ++no_prediction; continue; }
            for (const JsonValue& x : parr->array) external.push_back(x.number != 0.0 ? 1 : 0);
        }
        g_scan_obj = std::filesystem::absolute(obj).string(), g_scan_labels = std::filesystem::absolute(lab->second).string();
        const JsonResult j = parse_json(dataset::read_text(lab->second));
        const LoadResult r = parse_obj(dataset::read_text(obj));
        const JsonValue* arr = j.ok() ? j.value.find("labels") : nullptr;
        if (!r.ok() || !arr || arr->array.size() != r.mesh.positions.size()) { ++skipped; continue; }
        const AnalysisMesh analysis = manifold_analysis_mesh(r.mesh);
        if (!analysis.manifold) { ++skipped; continue; }
        const HalfEdgeMesh& m = analysis.halfedge;
        if (score && external.size() != m.positions.size()) { ++no_prediction; continue; }  // wrong vertex count
        std::vector<std::uint8_t> truth(m.positions.size());
        for (std::size_t v = 0; v < truth.size(); ++v) truth[v] = arr->array[v].number != 0.0 ? 1 : 0;
        std::vector<int> fdi(m.positions.size());
        for (std::size_t v = 0; v < fdi.size(); ++v) fdi[v] = static_cast<int>(arr->array[v].number);
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
            if (p.method != Method::HeightPlane && p.valley_weight > 0.0) {  // Voronoi and GraphCut
                auto it = valley.find(p.curvature_scale);
                if (it == valley.end()) it = valley.emplace(p.curvature_scale, valley_strength(ops, in.kmin, p.curvature_scale)).first;
                val = it->second;
            }
            std::vector<std::uint8_t> tooth;
            std::vector<std::uint8_t> veto;
            if (p.cusp_seed_quantile == kExternal) {
                tooth = external;
            } else if (p.cusp_seed_quantile == kHeightOracle || p.cusp_seed_quantile == kBothOracles) {
                // Oracles (diagnostics): veto tooth seeds on true gingiva. Height-only keeps cusp seeds
                // as the classifier left them; "both" vetoes every gingival tooth seed.
                MarginParams pp = p;
                pp.cusp_seed_quantile = 0.0;
                veto.assign(truth.size(), 0);
                for (std::size_t v = 0; v < truth.size(); ++v) veto[v] = truth[v] == 0;
                MarginInputs filtered = in;
                if (p.cusp_seed_quantile == kHeightOracle) {
                    // keep classifier-surviving cusp seeds even on gingiva: veto applies to height seeds only
                    std::vector<std::uint8_t> cusp_mark(truth.size(), 0);
                    for (std::uint32_t v : in.cusp_tips) cusp_mark[v] = 1;
                    for (std::size_t v = 0; v < truth.size(); ++v) if (cusp_mark[v]) veto[v] = 0;
                }
                pp.tooth_seed_veto = &veto;
                tooth = margin_labels(m, filtered, val, pp);
            } else if (p.cusp_seed_quantile == kOracle) {  // oracle: drop seeds on true gingiva (diagnostic only)
                MarginInputs filtered = in;
                std::erase_if(filtered.cusp_tips, [&](std::uint32_t v) { return truth[v] == 0; });
                MarginParams pp = p;
                pp.cusp_seed_quantile = 0.0;
                tooth = margin_labels(m, filtered, val, pp);
            } else {
                const auto tl = std::chrono::steady_clock::now();
                tooth = margin_labels(m, in, val, p);
                totals[c].label_ms += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - tl).count();
            }
            totals[c].add(compare_boundaries(edge_midpoints(m, label_boundary_edges(m, tooth)), gt_line),
                          region_iou(tooth, truth, ops.star0));
            totals[c].add_vertex(tooth, truth);
            if (g_experiment == "groups") group_eval(m, in, tooth, fdi);
            if (g_experiment.rfind("errors", 0) == 0) error_breakdown(m, in, tooth, truth, fdi, ops.star0);
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
    std::printf("\nper-vertex binary metrics (unweighted; for comparison with published per-point results):\n");
    for (std::size_t i = 0; i < configs.size(); ++i)
        std::printf("  %-60s accuracy %.4f, tooth IoU %.4f, gingiva IoU %.4f, mean of the two %.4f\n", configs[i].label().c_str(),
                    mean(totals[i].vacc), mean(totals[i].viou_tooth), mean(totals[i].viou_gingiva),
                    0.5 * (mean(totals[i].viou_tooth) + mean(totals[i].viou_gingiva)));
    if (g_experiment == "groups") {
        std::printf("\ncusp grouping (tips on true teeth; pairs: together / merged across teeth / split within a tooth):\n");
        for (std::size_t ti = 0; ti < kGroupTaus.size(); ++ti) {
            const GroupStats& st = g_groups[ti];
            const double pairs_same = double(st.together + st.split), pairs_grouped = double(st.together + st.merged);
            std::printf("  tau %5.2f: pair precision %.3f, pair recall %.3f, clean teeth %.1f%% of %zu, ungrouped tips %.1f%% (tips on gingiva %.1f%%)\n",
                        kGroupTaus[ti], pairs_grouped > 0 ? double(st.together) / pairs_grouped : 0.0,
                        pairs_same > 0 ? double(st.together) / pairs_same : 0.0, 100.0 * double(st.clean) / double(std::max<std::size_t>(st.teeth, 1)),
                        st.teeth, 100.0 * double(st.ungrouped) / double(std::max<std::size_t>(st.tips, 1)),
                        100.0 * double(st.tips_on_gingiva) / double(std::max<std::size_t>(st.tips, 1)));
        }
    }
    if (g_experiment.rfind("errors", 0) == 0) {
        const auto& e = g_errors;
        const double err = e.false_gingiva_missed + e.false_gingiva_partial + e.false_tooth;
        std::printf("\nerror breakdown (area-weighted, all scans): total error %.1f%% of arch area\n", 100.0 * err / e.total_area);
        std::printf("  false gingiva in MISSED teeth (<50%% covered): %.1f%% of error\n", 100.0 * e.false_gingiva_missed / err);
        std::printf("  false gingiva in partially covered teeth:      %.1f%% of error\n", 100.0 * e.false_gingiva_partial / err);
        std::printf("  false tooth (gingiva labelled tooth):          %.1f%% of error (%.1f%% within 0.5 mm of a true tooth)\n",
                    100.0 * e.false_tooth / err, 100.0 * e.false_tooth_near / err);
        std::printf("  missed teeth: %d (%d had a cusp seed on them)\n  miss rate by tooth type (FDI unit digit):",
                    e.missed_total, e.missed_with_seed);
        for (const auto& [t, mt] : e.by_type) std::printf(" %d: %d/%d", t, mt.first, mt.second);
        std::printf("\n");
        const double rg = double(e.gt_samples), pg = double(e.pr_samples);
        const double r_miss = double(e.r_missed_tooth + e.r_near + e.r_far), p_miss = double(e.p_fake + e.p_near + e.p_far);
        const double recall = 1.0 - r_miss / rg, precision = 1.0 - p_miss / pg;
        std::printf("\nF1@0.5 decomposition (pooled boundary samples): precision %.3f, recall %.3f, pooled F1 %.3f\n", precision,
                    recall, 2 * precision * recall / (precision + recall));
        std::printf("  recall misses (%.1f%% of true samples): on missed teeth %.1f%%, 0.5-1 mm off %.1f%%, > 1 mm off %.1f%%\n",
                    100 * r_miss / rg, 100 * double(e.r_missed_tooth) / rg, 100 * double(e.r_near) / rg, 100 * double(e.r_far) / rg);
        std::printf("    of the > 1 mm recall misses: gingiva side labelled tooth (over-extension) %.1f%%, interdental (two teeth"
                    " within 1.5 mm) %.1f%%, both %.1f%%\n", 100.0 * double(e.r_far_overext) / double(e.r_far),
                    100.0 * double(e.r_far_interdental) / double(e.r_far), 100.0 * double(e.r_far_overext_interdental) / double(e.r_far));
        auto q = [](std::vector<double> v, double f) {
            if (v.empty()) return 0.0;
            std::sort(v.begin(), v.end());
            return v[static_cast<std::size_t>(f * double(v.size() - 1))];
        };
        auto frac_below = [](const std::vector<double>& v, double t) {
            std::size_t n = 0;
            for (double x : v) n += x < t;
            return v.empty() ? 0.0 : double(n) / double(v.size());
        };
        std::printf("\ncrease strength along the TRUE boundary (kappa_min, 1/mm; more negative = stronger crease):\n");
        std::printf("  %-34s %8s %10s %10s %12s %14s\n", "", "edges", "edge med", "ring med", "ring q25", "ring < -1 /mm");
        const char* names[2][2] = {{"cheek/tongue side, matched", "cheek/tongue side, missed"}, {"interdental, matched", "interdental, missed"}};
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                std::printf("  %-34s %8zu %10.2f %10.2f %12.2f %13.1f%%\n", names[i][j], e.k_edge[i][j].size(), q(e.k_edge[i][j], 0.5),
                            q(e.k_ring[i][j], 0.5), q(e.k_ring[i][j], 0.25), 100.0 * frac_below(e.k_ring[i][j], -1.0));
        std::printf("  %-34s %8zu %10.2f %10s %12s %13.1f%%\n", "baseline: true gingiva vertices", e.k_gingiva.size(), q(e.k_gingiva, 0.5),
                    "-", "-", 100.0 * frac_below(e.k_gingiva, -1.0));
        std::printf("\nmissed interdental samples: distance to the neighbouring crown (gingiva strip width): q25 %.2f, median %.2f,"
                    " q75 %.2f mm; < 0.5 mm: %.1f%%\n", q(e.strip_width, 0.25), q(e.strip_width, 0.5), q(e.strip_width, 0.75),
                    100.0 * frac_below(e.strip_width, 0.5));
        auto cases = e.cases;
        std::sort(cases.begin(), cases.end(), [](const auto& x, const auto& y) { return x.missed > y.missed; });
        std::printf("\nworst scans by missed interdental boundary samples (viewer: npm run dev, then open the URL):\n");
        for (std::size_t i = 0; i < std::min<std::size_t>(cases.size(), 5); ++i)
            std::printf("  %zu missed  http://localhost:5173/?scan=%s&labels=%s&focus=%.2f,%.2f,%.2f&dist=10\n", cases[i].missed,
                        cases[i].obj.c_str(), cases[i].labels.c_str(), cases[i].focus.x, cases[i].focus.y, cases[i].focus.z);
        std::printf("  PREDICTED boundary (ring min kappa_min):\n");
        const char* pnames[2][2] = {{"cheek/tongue side, within 1 mm", "cheek/tongue side, > 1 mm off"}, {"interdental, within 1 mm", "interdental, > 1 mm off"}};
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                std::printf("  %-34s %8zu %10s %10.2f %12.2f %13.1f%%\n", pnames[i][j], e.kp_ring[i][j].size(), "-",
                            q(e.kp_ring[i][j], 0.5), q(e.kp_ring[i][j], 0.25), 100.0 * frac_below(e.kp_ring[i][j], -1.0));
        std::printf("  precision misses (%.1f%% of predicted samples): fake regions on gingiva %.1f%%, 0.5-1 mm off %.1f%%,"
                    " > 1 mm off (not fake) %.1f%%\n", 100 * p_miss / pg, 100 * double(e.p_fake) / pg, 100 * double(e.p_near) / pg,
                    100 * double(e.p_far) / pg);
    }
    if (g_experiment == "profile" && g_timings.calls) {
        const double n = double(g_timings.calls);
        std::printf("\nper-tooth cut, mean ms per scan: seeds %.0f, binary Dijkstras %.0f, binary cut %.0f, grouping %.0f, per-label"
                    " Dijkstras %.0f, unary %.0f, alpha-expansion %.0f; labels per scan %.1f\n", g_timings.seeds / n,
                    g_timings.binary_dijkstra / n, g_timings.binary_cut / n, g_timings.grouping / n, g_timings.label_dijkstra / n,
                    g_timings.unary / n, g_timings.expansion / n, double(g_timings.labels) / n);
    }
    std::printf("\nlabelling time per scan (ms, after shared inputs; native Release):");
    for (const auto& t : totals) std::printf(" %.0f", t.label_ms / double(std::max<std::size_t>(used, 1)));
    std::printf("\n");
    if (score) std::printf("scans without a usable prediction (skipped, not scored): %zu\n", no_prediction);
    if (mode == "validate" || score || (mode == "test" && g_seed_model)) {
        // Paired comparison on the same scans: is a gain real? Validate pairs every row with row 0; test
        // pairs every two rows (the disclosed runs compare old and new points with and without the classifier).
        auto paired = [&](std::size_t a, std::size_t b) {
            const auto &base = totals[a].assd, &t = totals[b].assd;
            double sum = 0.0, sq = 0.0;
            std::size_t wins = 0, losses = 0, big_losses = 0;
            double worst = 0.0;
            for (std::size_t i = 0; i < base.size(); ++i) {
                const double d = t[i] - base[i];
                sum += d, sq += d * d;
                wins += d < -1e-9, losses += d > 1e-9, big_losses += d > 0.1;
                worst = std::max(worst, d);
            }
            const double n = double(base.size()), mean_d = sum / n;
            const double se = std::sqrt(std::max(sq / n - mean_d * mean_d, 0.0) / (n - 1.0));
            std::printf("  [%zu] -> [%zu] %-50s mean %+.4f mm, SE %.4f (t = %+.1f), better on %zu, worse on %zu of %zu scans"
                        " (worse by > 0.1 mm: %zu; largest loss %+.3f mm)\n", a, b, configs[b].label().c_str(), mean_d, se,
                        se > 0 ? mean_d / se : 0.0, wins, losses, base.size(), big_losses, worst);
        };
        std::printf("\npaired per-scan ASSD difference (row index -> row index; negative = better):\n");
        for (std::size_t a = 0; a < configs.size(); ++a) {
            for (std::size_t b = a + 1; b < configs.size(); ++b) {
                if (mode == "test" || score || a == 0) paired(a, b);
            }
        }
    }
    if (mode == "test") {
        std::printf("\nlibrary path check (detect_margin, mean ASSD):");
        for (const auto& t : library) std::printf(" %.3f", mean(t.assd));
        std::printf("\n");
    }
    return 0;
}
