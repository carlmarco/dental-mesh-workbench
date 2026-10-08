// Trains the learned per-vertex data term (D88): P(tooth | local, label-free features) by logistic regression.
// Split of part-1 scans by index % 5: {1,2} = training, {0} = kept for tuning the weight in margin_eval,
// {3,4} = validation (AUC reported here). Per scan: 1,500 random arch vertices + 1,500 within 1.5 mm of the true
// boundary (where the cut actually decides). Writes the model JSON (local, git-ignored, like D76).
//   vertex_train <scan-dir> <model-out.json>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <limits>
#include <queue>
#include <random>
#include <string>
#include <vector>

#include "core/io.h"
#include "core/learn.h"
#include "core/margin.h"
#include "core/topology.h"
#include "dataset.h"

using namespace dmw;

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: vertex_train <scan-dir> <model-out.json>\n");
        return 2;
    }
    constexpr std::size_t K = kVertexFeatureCount, kRandom = 1500, kNear = 1500;
    const auto objs = dataset::index_files({argv[1]}, ".obj");
    const auto labels = dataset::index_files({argv[1]}, ".json");
    std::vector<double> xtr, xva;
    std::vector<int> ytr, yva;
    std::mt19937 rng(88);
    std::size_t index = 0, ntr = 0, nva = 0;
    for (const auto& [stem, obj] : objs) {
        const std::size_t split = index++ % 5;
        if (split == 0) continue;  // tuning scans: not used for fitting
        const bool train = split <= 2;
        const auto lab = labels.find(stem);
        if (lab == labels.end()) continue;
        const JsonResult j = parse_json(dataset::read_text(lab->second));
        const LoadResult r = parse_obj(dataset::read_text(obj));
        const JsonValue* arr = j.ok() ? j.value.find("labels") : nullptr;
        if (!r.ok() || !arr || arr->array.size() != r.mesh.positions.size()) continue;
        const AnalysisMesh a = manifold_analysis_mesh(r.mesh);
        if (!a.manifold) continue;
        const HalfEdgeMesh& m = a.halfedge;
        const MarginInputs in = margin_inputs(m, true);
        MarginParams p = margin_operating_point();
        p.method = MarginParams::Method::GraphCut;  // features come from the binary stage; skip the per-tooth cost
        std::vector<double> feat;
        p.vertex_features_out = &feat;
        margin_labels(m, in, valley_strength(build_dec(m), in.kmin, p.curvature_scale), p);
        const std::size_t nv = m.positions.size();
        std::vector<std::uint8_t> truth(nv);
        for (std::size_t v = 0; v < nv; ++v) truth[v] = arr->array[v].number != 0.0 ? 1 : 0;
        // Vertices within 1.5 mm (geodesic, edge lengths) of the true boundary.
        std::vector<double> d(nv, std::numeric_limits<double>::infinity());
        using Item = std::pair<double, std::uint32_t>;
        std::priority_queue<Item, std::vector<Item>, std::greater<>> heap;
        for (const Edge& e : label_boundary_edges(m, truth)) {
            for (std::uint32_t v : {e.v0, e.v1})
                if (d[v] > 0.0) d[v] = 0.0, heap.push({0.0, v});
        }
        while (!heap.empty()) {
            const auto [dv, v] = heap.top();
            heap.pop();
            if (dv > d[v] || dv > 1.5) continue;
            for (std::uint32_t w : one_ring(m, v)) {
                const Vec3 e{m.positions[w].x - m.positions[v].x, m.positions[w].y - m.positions[v].y, m.positions[w].z - m.positions[v].z};
                const double nd = dv + std::sqrt(e.x * e.x + e.y * e.y + e.z * e.z);
                if (nd < d[w]) d[w] = nd, heap.push({nd, w});
            }
        }
        std::vector<std::uint32_t> all, near;
        for (std::uint32_t v = 0; v < nv; ++v) {
            if (!in.arch[v]) continue;
            all.push_back(v);
            if (d[v] <= 1.5) near.push_back(v);
        }
        std::shuffle(all.begin(), all.end(), rng);
        std::shuffle(near.begin(), near.end(), rng);
        all.resize(std::min(all.size(), kRandom));
        near.resize(std::min(near.size(), kNear));
        auto& x = train ? xtr : xva;
        auto& y = train ? ytr : yva;
        for (const auto* set : {&all, &near})
            for (std::uint32_t v : *set) {
                x.insert(x.end(), feat.begin() + static_cast<std::ptrdiff_t>(v * K), feat.begin() + static_cast<std::ptrdiff_t>((v + 1) * K));
                y.push_back(truth[v]);
            }
        (train ? ntr : nva)++;
        std::fprintf(stderr, "\r%zu train / %zu validation scans", ntr, nva);
    }
    std::fprintf(stderr, "\n");
    auto tooth_rate = [](const std::vector<int>& y) {
        std::size_t t = 0;
        for (int v : y) t += v != 0;
        return double(t) / double(std::max<std::size_t>(y.size(), 1));
    };
    std::printf("samples: %zu training (%zu scans, %.1f%% tooth), %zu validation (%zu scans, %.1f%% tooth)\n", ytr.size(), ntr,
                100 * tooth_rate(ytr), yva.size(), nva, 100 * tooth_rate(yva));
    LogisticModel model = train_logistic(xtr, ytr, K, 1e-3);
    model.features = vertex_feature_names();
    auto scores = [&](const std::vector<double>& x) {
        std::vector<double> s;
        for (std::size_t i = 0; i * K < x.size(); ++i) s.push_back(model.probability(std::span<const double>(x).subspan(i * K, K)));
        return s;
    };
    std::printf("ROC AUC (tooth vs gingiva vertex): training %.4f, validation %.4f\n", roc_auc(scores(xtr), ytr), roc_auc(scores(xva), yva));
    for (std::size_t f = 0; f < K; ++f) {
        std::vector<double> col;
        for (std::size_t i = f; i < xva.size(); i += K) col.push_back(xva[i]);
        const double auc = roc_auc(col, yva);
        std::printf("  %-26s weight %+8.3f (standardized) | single-feature validation AUC %.3f\n", model.features[f].c_str(),
                    model.weights[f], std::max(auc, 1.0 - auc));
    }
    std::ofstream(argv[2]) << model.to_json() << '\n';
    std::printf("model written to %s\n", argv[2]);
    return 0;
}
