// Margin failure analysis (D74). Runs the margin operating point on each scan and records where and
// how it fails. Use VALIDATION scans (training scans not used in the sweeps) to guide changes; the
// held-out test set is only ever described, never used for tuning.
//   margin_analyze <scan-dir> <out.csv> [stride] [min_remainder] [model.json]
//   (keeps scans with index % stride >= min_remainder; optional seed classifier at threshold 0.5, D76)
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "core/cusps.h"
#include "core/dec.h"
#include "core/io.h"
#include "core/learn.h"
#include "core/margin.h"
#include "core/topology.h"
#include "dataset.h"

using namespace dmw;

namespace {

const char* kCategory[4] = {"incisor", "canine", "premolar", "molar"};
std::size_t category(int fdi) {
    const int t = fdi % 10;
    return t <= 2 ? 0 : t == 3 ? 1 : t <= 5 ? 2 : 3;
}

// Brute-force nearest distances (margins have a few thousand points; this is analysis, not runtime).
std::vector<double> nearest(const std::vector<Vec3>& from, const std::vector<Vec3>& to) {
    std::vector<double> d(from.size(), 1e9);
    for (std::size_t i = 0; i < from.size(); ++i) {
        for (const Vec3& q : to) {
            const double dx = from[i].x - q.x, dy = from[i].y - q.y, dz = from[i].z - q.z;
            d[i] = std::min(d[i], dx * dx + dy * dy + dz * dz);
        }
        d[i] = std::sqrt(d[i]);
    }
    return d;
}

double mean(const std::vector<double>& v) {
    double s = 0.0;
    for (double x : v) s += x;
    return v.empty() ? 0.0 : s / double(v.size());
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: margin_analyze <scan-dir> <out.csv> [stride] [offset]\n");
        return 2;
    }
    const std::size_t stride = argc > 3 ? std::stoul(argv[3]) : 0, min_rem = argc > 4 ? std::stoul(argv[4]) : 0;
    LogisticModel seed_model;
    if (argc > 5) {
        std::string err;
        seed_model = parse_logistic_model(dataset::read_text(argv[5]), err);
        if (!err.empty()) return std::fprintf(stderr, "model: %s\n", err.c_str()), 2;
    }
    MarginParams params = margin_operating_point();
    if (!seed_model.weights.empty()) params.seed_model = &seed_model, params.seed_threshold = 0.5;
    const auto objs = dataset::index_files({argv[1]}, ".obj");
    const auto labels = dataset::index_files({argv[1]}, ".json");
    std::ofstream csv(argv[2]);
    csv << "scan,vertices,teeth,assd,hd95,f1_050,iou,false_margin_mm,missed_margin_mm,over_seg,under_seg,"
           "excluded_faces,components,cusp_seeds,worst_tooth,worst_tooth_mm,seeds_on_gingiva,mean_edge_mm,genus,"
           "false_far_frac,false_near_frac,height_seeds_on_gingiva\n";

    std::array<std::vector<double>, 4> cat_err;  // per-point missed-margin distances by tooth category
    std::size_t index = 0, n = 0;
    for (const auto& [stem, obj] : objs) {
        const std::size_t i = index++;
        if (stride && i % stride < min_rem) continue;
        const auto lab = labels.find(stem);
        if (lab == labels.end()) continue;
        const JsonResult j = parse_json(dataset::read_text(lab->second));
        const LoadResult r = parse_obj(dataset::read_text(obj));
        const JsonValue* arr = j.ok() ? j.value.find("labels") : nullptr;
        if (!r.ok() || !arr || arr->array.size() != r.mesh.positions.size()) continue;
        const AnalysisMesh a = manifold_analysis_mesh(r.mesh);
        if (!a.manifold) continue;
        const HalfEdgeMesh& m = a.halfedge;
        const std::size_t nv = m.positions.size();
        std::vector<int> fdi(nv);
        std::vector<std::uint8_t> truth(nv);
        std::map<int, int> teeth;
        for (std::size_t v = 0; v < nv; ++v) {
            fdi[v] = static_cast<int>(arr->array[v].number);
            truth[v] = fdi[v] != 0;
            if (fdi[v]) teeth[fdi[v]]++;
        }
        const MarginResult pred = detect_margin(m, params);
        const auto gt_edges = label_boundary_edges(m, truth);
        const auto gt_pts = edge_midpoints(m, gt_edges), pr_pts = edge_midpoints(m, pred.margin);
        const BoundaryMetrics b = compare_boundaries(pr_pts, gt_pts);
        const auto area = build_dec(m).star0;
        double over = 0.0, under = 0.0, tooth_area = 0.0;
        for (std::size_t v = 0; v < nv; ++v) {
            if (truth[v]) tooth_area += area[v];
            if (pred.tooth[v] && !truth[v]) over += area[v];
            if (!pred.tooth[v] && truth[v]) under += area[v];
        }
        // Missed margin per tooth: each GT margin point belongs to the tooth on its tooth side.
        const auto missed = nearest(gt_pts, pr_pts);
        const auto falsem = nearest(pr_pts, gt_pts);
        std::map<int, std::vector<double>> per_tooth;
        for (std::size_t e = 0; e < gt_edges.size(); ++e) {
            const int t = fdi[gt_edges[e].v0] ? fdi[gt_edges[e].v0] : fdi[gt_edges[e].v1];
            per_tooth[t].push_back(missed[e]);
            cat_err[category(t)].push_back(missed[e]);
        }
        int worst = 0;
        double worst_mm = 0.0;
        for (const auto& [t, d] : per_tooth) {
            if (mean(d) > worst_mm) worst_mm = mean(d), worst = t;
        }
        const TopologyReport topo = analyze_topology(a.mesh);
        const std::size_t comps = topo.components.size();
        std::uint32_t genus = 0;
        for (const auto& c : topo.components) genus += c.genus.value_or(0);
        // H1: cusp seeds that land on true gingiva start tooth fronts inside the gingiva.
        const auto seeds = detect_cusps(m, cusp_operating_point()).vertices;
        const std::size_t cusps = seeds.size();
        std::size_t on_gingiva = 0;
        for (std::uint32_t v : seeds) on_gingiva += truth[v] == 0;
        // D77: split false-margin samples (> 1 mm from the true margin) by distance to the nearest true
        // TOOTH vertex: > 2 mm = a spurious tooth region inside gingiva (e.g. distal loops), else an offset.
        std::vector<Vec3> tooth_pts;
        for (std::size_t v = 0; v < nv; ++v)
            if (truth[v]) tooth_pts.push_back(m.positions[v]);
        const auto to_tooth = nearest_point_distances(pr_pts, tooth_pts, 1.0);
        std::size_t far = 0, near = 0;
        for (std::size_t k = 0; k < pr_pts.size(); ++k) {
            if (falsem[k] <= 1.0) continue;
            (to_tooth[k] > 2.0 ? far : near)++;
        }
        // Height seeds (top tooth_quantile of arch height) that land on true gingiva.
        const MarginInputs mi = margin_inputs(m, false);
        const auto& hs = mi.sorted_arch_height;
        const double high = hs.empty() ? 0.0 : hs[static_cast<std::size_t>(params.tooth_quantile * double(hs.size() - 1))];
        std::size_t hseeds = 0, hseeds_gingiva = 0;
        for (std::size_t v = 0; v < nv; ++v) {
            if (!mi.arch[v] || mi.height[v] < high) continue;
            ++hseeds;
            hseeds_gingiva += truth[v] == 0;
        }
        // H2: resolution (mean edge length).
        double edge_sum = 0.0;
        std::size_t edge_n = 0;
        for (std::uint32_t h = 0; h < m.origin.size(); ++h) {
            const Vec3 &p = m.positions[m.origin[h]], &q = m.positions[dest(m, h)];
            edge_sum += std::sqrt((p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y) + (p.z - q.z) * (p.z - q.z));
            ++edge_n;
        }
        csv << stem << ',' << nv << ',' << teeth.size() << ',' << b.assd << ',' << b.hd95 << ',' << b.f1_050 << ','
            << region_iou(pred.tooth, truth, area) << ',' << mean(falsem) << ',' << mean(missed) << ','
            << over / tooth_area << ',' << under / tooth_area << ',' << a.excluded_faces << ',' << comps << ','
            << cusps << ',' << worst << ',' << worst_mm << ',' << (cusps ? double(on_gingiva) / double(cusps) : 0.0) << ','
            << edge_sum / double(std::max<std::size_t>(edge_n, 1)) << ',' << genus << ','
            << double(far) / double(std::max<std::size_t>(pr_pts.size(), 1)) << ',' << double(near) / double(std::max<std::size_t>(pr_pts.size(), 1))
            << ',' << (hseeds ? double(hseeds_gingiva) / double(hseeds) : 0.0) << '\n';
        ++n;
        std::fprintf(stderr, "\r%zu scans", n);
    }
    std::fprintf(stderr, "\n");
    std::printf("scans analysed: %zu\nmissed-margin distance by tooth category (all margin points pooled):\n", n);
    for (std::size_t c = 0; c < 4; ++c) {
        auto d = cat_err[c];
        std::sort(d.begin(), d.end());
        if (d.empty()) continue;
        std::printf("  %-9s mean %.3f mm, median %.3f, p95 %.3f, points %zu\n", kCategory[c], mean(d), d[d.size() / 2],
                    d[static_cast<std::size_t>(0.95 * double(d.size() - 1))], d.size());
    }
    return 0;
}
