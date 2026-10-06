// Trains the cusp-seed classifier (D76): is a detected cusp seed on a tooth (1) or on gingiva (0)?
// Split of part-1 scans by index % 5: {0,1,2} = classifier training, {3,4} = validation.
// Writes the model JSON (local, git-ignored: data/models; publication pending the author's licence decision).
//   seed_train <scan-dir> <model-out.json>
#include <cstdio>
#include <fstream>
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
        std::fprintf(stderr, "usage: seed_train <scan-dir> <model-out.json>\n");
        return 2;
    }
    const auto objs = dataset::index_files({argv[1]}, ".obj");
    const auto labels = dataset::index_files({argv[1]}, ".json");
    std::vector<double> xtr, xva;
    std::vector<int> ytr, yva;
    std::size_t index = 0, ntr = 0, nva = 0;
    for (const auto& [stem, obj] : objs) {
        const std::size_t split = index++ % 5;
        const auto lab = labels.find(stem);
        if (lab == labels.end()) continue;
        const JsonResult j = parse_json(dataset::read_text(lab->second));
        const LoadResult r = parse_obj(dataset::read_text(obj));
        const JsonValue* arr = j.ok() ? j.value.find("labels") : nullptr;
        if (!r.ok() || !arr || arr->array.size() != r.mesh.positions.size()) continue;
        const AnalysisMesh a = manifold_analysis_mesh(r.mesh);
        if (!a.manifold) continue;
        const MarginInputs in = margin_inputs(a.halfedge, true);
        auto& x = split <= 2 ? xtr : xva;
        auto& y = split <= 2 ? ytr : yva;
        for (std::size_t i = 0; i < in.cusp_tips.size(); ++i) {
            x.insert(x.end(), in.seed_features.begin() + static_cast<std::ptrdiff_t>(i * kSeedFeatureCount),
                     in.seed_features.begin() + static_cast<std::ptrdiff_t>((i + 1) * kSeedFeatureCount));
            y.push_back(arr->array[in.cusp_tips[i]].number != 0.0 ? 1 : 0);
        }
        (split <= 2 ? ntr : nva)++;
        std::fprintf(stderr, "\r%zu train / %zu validation scans", ntr, nva);
    }
    std::fprintf(stderr, "\n");
    auto rate = [](const std::vector<int>& y) {
        std::size_t g = 0;
        for (int v : y) g += v == 0;
        return double(g) / double(y.size());
    };
    std::printf("seeds: %zu training (%zu scans, %.1f%% on gingiva), %zu validation (%zu scans, %.1f%% on gingiva)\n", ytr.size(),
                ntr, 100 * rate(ytr), yva.size(), nva, 100 * rate(yva));
    LogisticModel m = train_logistic(xtr, ytr, kSeedFeatureCount, 1e-3);
    m.features = seed_feature_names();
    auto scores = [&](const std::vector<double>& x) {
        std::vector<double> s;
        for (std::size_t i = 0; i * kSeedFeatureCount < x.size(); ++i)
            s.push_back(m.probability(std::span<const double>(x).subspan(i * kSeedFeatureCount, kSeedFeatureCount)));
        return s;
    };
    std::printf("ROC AUC (tooth vs gingiva seed): training %.3f, validation %.3f\n", roc_auc(scores(xtr), ytr), roc_auc(scores(xva), yva));
    // Single-feature AUCs on validation: which signal carries the discrimination?
    for (std::size_t f = 0; f < kSeedFeatureCount; ++f) {
        std::vector<double> col;
        for (std::size_t i = f; i < xva.size(); i += kSeedFeatureCount) col.push_back(xva[i]);
        const double auc = roc_auc(col, yva);
        std::printf("  %-28s weight %+7.3f (standardized) | single-feature validation AUC %.3f\n", m.features[f].c_str(), m.weights[f],
                    std::max(auc, 1.0 - auc));
    }
    std::ofstream(argv[2]) << m.to_json() << '\n';
    std::printf("model written to %s\n", argv[2]);
    return 0;
}
