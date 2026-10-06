#include <cmath>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/learn.h"

using namespace dmw;
using Catch::Matchers::WithinAbs;

namespace {
double hash01(double k) {
    const double s = std::sin(k * 12.9898) * 43758.5453;
    return s - std::floor(s);
}
}  // namespace

TEST_CASE("learn: recovers a known 1D logistic relationship", "[learn]") {
    // Labels drawn from P(y=1|x) = sigmoid(3x - 1) with a deterministic hash: the fit should
    // recover slope and intercept in raw units (after undoing the standardization).
    std::vector<double> x;
    std::vector<int> y;
    for (int i = 0; i < 20000; ++i) {
        const double xi = 4.0 * hash01(i) - 2.0;
        const double p = 1.0 / (1.0 + std::exp(-(3.0 * xi - 1.0)));
        x.push_back(xi);
        y.push_back(hash01(i + 0.5) < p ? 1 : 0);
    }
    const auto m = train_logistic(x, y, 1, 1e-6);
    const double slope = m.weights[0] / m.scale[0];
    const double intercept = m.bias - m.weights[0] * m.mean[0] / m.scale[0];
    CHECK_THAT(slope, WithinAbs(3.0, 0.15));
    CHECK_THAT(intercept, WithinAbs(-1.0, 0.1));
}

TEST_CASE("learn: a 2D problem where only one feature matters", "[learn]") {
    std::vector<double> x;
    std::vector<int> y;
    for (int i = 0; i < 4000; ++i) {
        const double a = hash01(2 * i) - 0.5, noise = hash01(2 * i + 1) - 0.5;
        x.insert(x.end(), {a, noise});
        y.push_back(a > 0.1 ? 1 : 0);
    }
    const auto m = train_logistic(x, y, 2, 1e-3);
    CHECK(std::abs(m.weights[0]) > 10.0 * std::abs(m.weights[1]));
    CHECK(m.probability(std::vector<double>{0.4, 0.0}) > 0.95);
    CHECK(m.probability(std::vector<double>{-0.4, 0.0}) < 0.05);
}

TEST_CASE("learn: JSON round trip preserves predictions", "[learn]") {
    std::vector<double> x{0, 1, 2, 3, 4, 5};
    std::vector<int> y{0, 0, 0, 1, 1, 1};
    auto m = train_logistic(x, y, 1, 1.0);
    m.features = {"height"};
    std::string err;
    const auto back = parse_logistic_model(m.to_json(), err);
    REQUIRE(err.empty());
    CHECK(back.features == m.features);
    CHECK_THAT(back.probability(std::vector<double>{2.5}), WithinAbs(m.probability(std::vector<double>{2.5}), 1e-12));
    std::string err2;
    CHECK(parse_logistic_model("{\"weights\":[1]}", err2).weights.empty());  // inconsistent model rejected
    CHECK_FALSE(err2.empty());
}

TEST_CASE("learn: ROC AUC on perfect, inverted, and tied scores", "[learn]") {
    const std::vector<int> lab{0, 0, 1, 1};
    CHECK(roc_auc(std::vector<double>{0.1, 0.2, 0.8, 0.9}, lab) == 1.0);
    CHECK(roc_auc(std::vector<double>{0.9, 0.8, 0.2, 0.1}, lab) == 0.0);
    CHECK(roc_auc(std::vector<double>{0.5, 0.5, 0.5, 0.5}, lab) == 0.5);
}
