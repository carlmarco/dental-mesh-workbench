// Cusp detection evaluation against 3DTeethLand Cusp landmarks (D67, D68).
//   cusp_eval sweep <landmark-dir> <scan-dir>...   parameter sweep (use the TRAINING landmarks)
//   cusp_eval test  <landmark-dir> <scan-dir>...   fixed operating points (use the TEST landmarks)
// Score fields are computed once per scan and reused across configurations; the chosen points are
// re-run through detect_cusps() itself to confirm the fast path matches the library.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "core/curvature.h"
#include "core/cusps.h"
#include "core/io.h"
#include "core/metrics.h"
#include "core/topology.h"
#include "dataset.h"

using namespace dmw;
using Clock = std::chrono::steady_clock;
using Method = CuspParams::Method;

namespace {

const double kTolerances[3] = {0.5, 1.0, 2.0};

struct Config {
    CuspParams p;
    std::string label() const {
        char b[160];
        const char* m = p.method == Method::RawCurvature ? "A raw H" : p.method == Method::SmoothedCurvature ? "B smooth H" : "C prominence";
        std::snprintf(b, sizeof b, "%-12s sigma=%.1f R=%.1f q=%.1f thr=%.2f nms=%.1f", m, p.curvature_scale,
                      p.prominence_scale, p.height_quantile, p.threshold, p.nms_radius);
        return b;
    }
};

std::vector<Config> sweep_configs() {
    // Ranges widened until each method's best F1 lies strictly inside the grid (D68). History:
    // sweep 1: optimum at the grid edge for all methods. Sweep 2: interior for A (thr 2.0, nms 2.0),
    // edge for B (sigma 1.8, nms 3.0) and C (thr 0.8, nms 3.0). Sweep 3: C interior in thr/nms/q but
    // R = 4.0 at the edge; B thr 0.2 and nms 4.0 at the edge. Sweep 4 (this): extends those.
    std::vector<Config> c;
    for (double thr : {1.2, 2.0, 3.0})
        for (double nms : {1.5, 2.0, 3.0}) c.push_back({{Method::RawCurvature, 0.5, 2.0, 0.5, nms, thr}});
    for (double sg : {1.8, 2.5, 3.5})
        for (double thr : {0.1, 0.15, 0.2, 0.3})
            for (double nms : {3.5, 4.0, 4.5, 5.0}) c.push_back({{Method::SmoothedCurvature, sg, 2.0, 0.5, nms, thr}});
    for (double R : {2.5, 4.0, 6.0, 8.0})
        for (double thr : {1.0, 1.3, 1.6, 2.0})
            for (double nms : {3.0, 3.5, 4.0}) c.push_back({{Method::OcclusalProminence, 0.5, R, 0.5, nms, thr}});
    return c;
}

// Fixed operating points for the test set: the best F1 at 1 mm per method on the TRAINING sweep
// (67 scans, sweep 4, 2026-10-05). Edit only from training results, never from test results.
//   A: F1 0.101 (thr 2.0, nms 2.0; interior)
//   B: F1 0.596 (sigma 3.5 at the edge of its range, +0.007 over sigma 2.5: diminishing; thr 0.15, nms 4.0)
//   C: F1 0.663 (sigma 0.5, R 4.0, q 0.5, thr 1.3, nms 3.5; interior in every swept parameter)
std::vector<Config> test_configs() {
    return {
        {{Method::RawCurvature, 0.5, 2.0, 0.5, 2.0, 2.0}},
        {{Method::SmoothedCurvature, 3.5, 2.0, 0.5, 4.0, 0.15}},
        {{Method::OcclusalProminence, 0.5, 4.0, 0.5, 3.5, 1.3}},
    };
}

// Per-scan score fields, computed once (keyed by the parameters that affect them).
struct Fields {
    std::vector<std::uint8_t> arch;
    std::vector<double> raw_h, height, smooth_h_05;
    std::map<double, std::vector<double>> smooth_h;  // by sigma
    std::map<double, std::vector<double>> base;      // prominence by R (D70)
    std::vector<double> sorted_arch_height;
};

std::vector<double> score_for(const HalfEdgeMesh& m, Fields& f, const DecOperators& ops, const CuspParams& p) {
    const std::size_t nv = m.positions.size();
    const double neg = -std::numeric_limits<double>::infinity();
    std::vector<double> score(nv, neg);
    auto smooth = [&](double sg) -> const std::vector<double>& {
        auto it = f.smooth_h.find(sg);
        if (it == f.smooth_h.end()) it = f.smooth_h.emplace(sg, diffuse(ops, f.raw_h, sg, 1e-4)).first;  // as detect_cusps
        return it->second;
    };
    if (p.method == Method::RawCurvature) {
        for (std::size_t v = 0; v < nv; ++v)
            if (f.arch[v]) score[v] = f.raw_h[v];
    } else if (p.method == Method::SmoothedCurvature) {
        const auto& s = smooth(p.curvature_scale);
        for (std::size_t v = 0; v < nv; ++v)
            if (f.arch[v]) score[v] = s[v];
    } else {
        const auto& s = smooth(p.curvature_scale);
        auto it = f.base.find(p.prominence_scale);
        if (it == f.base.end()) it = f.base.emplace(p.prominence_scale, prominence(ops, f.height, p.prominence_scale, 1e-4)).first;
        const auto& hs = f.sorted_arch_height;
        const double gate = hs[static_cast<std::size_t>(p.height_quantile * double(hs.size() - 1))];
        for (std::size_t v = 0; v < nv; ++v)
            if (f.arch[v] && s[v] > 0.0 && f.height[v] >= gate) score[v] = it->second[v];
    }
    return score;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        std::fprintf(stderr, "usage: cusp_eval sweep|test <landmark-dir> <scan-dir>...\n");
        return 2;
    }
    const std::string mode = argv[1];
    const auto landmarks = dataset::index_files({argv[2]}, "__kpt.json");
    std::vector<std::string> scan_dirs(argv + 3, argv + argc);
    const auto scans = dataset::index_files(scan_dirs, ".obj");
    const std::vector<Config> configs = mode == "sweep" ? sweep_configs() : test_configs();

    std::vector<std::array<MatchResult, 3>> totals(configs.size());
    std::vector<MatchResult> library_check(configs.size());  // detect_cusps() at tolerance 1 mm
    int used = 0, skipped = 0;
    double feature_ms = 0.0;
    for (const auto& [stem, lm_path] : landmarks) {
        const auto it = scans.find(stem);
        if (it == scans.end()) continue;
        const std::vector<Vec3> gt = dataset::read_landmarks(lm_path, "Cusp");
        const LoadResult r = parse_obj(dataset::read_text(it->second));
        if (!r.ok()) { ++skipped; continue; }
        const AnalysisMesh analysis = manifold_analysis_mesh(r.mesh);  // D69
        if (!analysis.manifold) { ++skipped; continue; }
        const HalfEdgeMesh& m = analysis.halfedge;
        ++used;

        const auto t0 = Clock::now();
        Fields f;
        f.arch = largest_component_mask(m);
        const Vec3 axis = occlusal_axis(m, f.arch);
        const CurvatureField curv = compute_curvature(m);
        f.raw_h.resize(m.positions.size());
        f.height.resize(m.positions.size());
        for (std::size_t v = 0; v < m.positions.size(); ++v) {
            f.raw_h[v] = std::isfinite(curv.mean[v]) ? curv.mean[v] : 0.0;
            f.height[v] = m.positions[v].x * axis.x + m.positions[v].y * axis.y + m.positions[v].z * axis.z;
            if (f.arch[v]) f.sorted_arch_height.push_back(f.height[v]);
        }
        std::sort(f.sorted_arch_height.begin(), f.sorted_arch_height.end());
        const DecOperators ops = build_dec(m);
        feature_ms += std::chrono::duration<double, std::milli>(Clock::now() - t0).count();

        for (std::size_t c = 0; c < configs.size(); ++c) {
            const auto score = score_for(m, f, ops, configs[c].p);
            const auto cand = local_maxima(m, score, configs[c].p.threshold);
            const auto kept = non_max_suppression(m.positions, score, cand, configs[c].p.nms_radius);
            std::vector<Vec3> det;
            for (std::uint32_t v : kept) det.push_back(m.positions[v]);
            for (std::size_t t = 0; t < 3; ++t) totals[c][t].add(match_points(det, gt, kTolerances[t]));
            if (mode == "test") {  // the library path, end to end
                std::vector<Vec3> lib;
                for (std::uint32_t v : detect_cusps(m, configs[c].p).vertices) lib.push_back(m.positions[v]);
                library_check[c].add(match_points(lib, gt, 1.0));
            }
        }
        std::fprintf(stderr, "\r%d scans", used);
    }
    std::fprintf(stderr, "\n");
    std::printf("scans evaluated: %d (skipped: %d), mean feature time %.0f ms/scan\n\n", used, skipped, feature_ms / std::max(used, 1));
    std::vector<std::size_t> order(configs.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    if (mode == "sweep") {
        std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return totals[a][1].f1() > totals[b][1].f1(); });
    }
    std::printf("| configuration | P@1mm | R@1mm | F1@0.5 | F1@1 | F1@2 | median err (mm, @1) | detections/scan |\n|---|---:|---:|---:|---:|---:|---:|---:|\n");
    std::map<Method, int> shown;
    for (std::size_t i : order) {
        const auto& t = totals[i];
        if (mode == "sweep" && shown[configs[i].p.method]++ >= 3) continue;  // top 3 per method
        auto d = t[1].matched_distances;
        std::sort(d.begin(), d.end());
        std::printf("| %s | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.1f |\n", configs[i].label().c_str(), t[1].precision(),
                    t[1].recall(), t[0].f1(), t[1].f1(), t[2].f1(), d.empty() ? 0.0 : d[d.size() / 2],
                    double(t[1].detections) / std::max(used, 1));
    }
    if (mode == "test") {
        std::printf("\nlibrary path check (detect_cusps, F1@1mm):");
        for (std::size_t c = 0; c < configs.size(); ++c) std::printf(" %.3f", library_check[c].f1());
        std::printf("\n");
    }
    return 0;
}
