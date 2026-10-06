#pragma once

#include <cstddef>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace dmw {

// Binary logistic regression (D76): P(y = 1 | x) = sigmoid(w . z + b), z = standardized features.
// Trained by Newton's method (iteratively reweighted least squares) with an L2 penalty on w.
// Small, dependency-free, and explainable: one weight per feature.
struct LogisticModel {
    std::vector<std::string> features;  // names, for reporting and file round-trips
    std::vector<double> mean, scale;    // standardization: z = (x - mean) / scale
    std::vector<double> weights;
    double bias = 0.0;

    double probability(std::span<const double> x) const;
    std::string to_json() const;
};

// rows: n samples, each `dim` features (row-major in `x`); labels 0/1.
LogisticModel train_logistic(std::span<const double> x, std::span<const int> y, std::size_t dim, double l2 = 1e-2,
                             int iterations = 25);

// Parses to_json() output; on failure returns a model with no weights and sets `error`.
LogisticModel parse_logistic_model(std::string_view json, std::string& error);

// Area under the ROC curve (rank statistic; ties count half).
double roc_auc(std::span<const double> score, std::span<const int> label);

}  // namespace dmw
