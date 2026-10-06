#include "core/learn.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>

#include "core/json.h"

namespace dmw {
namespace {

double sigmoid(double t) { return t >= 0 ? 1.0 / (1.0 + std::exp(-t)) : std::exp(t) / (1.0 + std::exp(t)); }

// Solves the small dense SPD system A u = r in place (Cholesky); A is n x n row-major.
bool solve_spd(std::vector<double>& a, std::vector<double>& r, std::size_t n) {
    for (std::size_t j = 0; j < n; ++j) {
        double d = a[j * n + j];
        for (std::size_t k = 0; k < j; ++k) d -= a[j * n + k] * a[j * n + k];
        if (!(d > 0.0)) return false;
        d = std::sqrt(d);
        a[j * n + j] = d;
        for (std::size_t i = j + 1; i < n; ++i) {
            double s = a[i * n + j];
            for (std::size_t k = 0; k < j; ++k) s -= a[i * n + k] * a[j * n + k];
            a[i * n + j] = s / d;
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        double s = r[i];
        for (std::size_t k = 0; k < i; ++k) s -= a[i * n + k] * r[k];
        r[i] = s / a[i * n + i];
    }
    for (std::size_t i = n; i-- > 0;) {
        double s = r[i];
        for (std::size_t k = i + 1; k < n; ++k) s -= a[k * n + i] * r[k];
        r[i] = s / a[i * n + i];
    }
    return true;
}

}  // namespace

double LogisticModel::probability(std::span<const double> x) const {
    double t = bias;
    for (std::size_t j = 0; j < weights.size(); ++j) t += weights[j] * (x[j] - mean[j]) / scale[j];
    return sigmoid(t);
}

std::string LogisticModel::to_json() const {
    std::ostringstream o;
    o.precision(17);
    auto list = [&](const std::vector<double>& v) {
        o << '[';
        for (std::size_t i = 0; i < v.size(); ++i) o << (i ? "," : "") << v[i];
        o << ']';
    };
    o << "{\"model\":\"logistic\",\"features\":[";
    for (std::size_t i = 0; i < features.size(); ++i) o << (i ? "," : "") << '"' << features[i] << '"';
    o << "],\"mean\":";
    list(mean);
    o << ",\"scale\":";
    list(scale);
    o << ",\"weights\":";
    list(weights);
    o << ",\"bias\":" << bias << "}";
    return o.str();
}

LogisticModel train_logistic(std::span<const double> x, std::span<const int> y, std::size_t dim, double l2, int iterations) {
    const std::size_t n = y.size();
    LogisticModel m;
    m.mean.assign(dim, 0.0);
    m.scale.assign(dim, 1.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < dim; ++j) m.mean[j] += x[i * dim + j];
    for (double& v : m.mean) v /= double(std::max<std::size_t>(n, 1));
    for (std::size_t j = 0; j < dim; ++j) {
        double var = 0.0;
        for (std::size_t i = 0; i < n; ++i) var += (x[i * dim + j] - m.mean[j]) * (x[i * dim + j] - m.mean[j]);
        const double sd = std::sqrt(var / double(std::max<std::size_t>(n, 1)));
        m.scale[j] = sd > 0.0 ? sd : 1.0;
    }
    // Parameters theta = (w_1..w_dim, b). Newton step: (X^T S X + L) d = X^T (y - p) - L theta,
    // S = diag(p (1 - p)), L = l2 * n on the weights (not the bias).
    const std::size_t k = dim + 1;
    std::vector<double> theta(k, 0.0), z(k);
    for (int it = 0; it < iterations; ++it) {
        std::vector<double> hess(k * k, 0.0), grad(k, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < dim; ++j) z[j] = (x[i * dim + j] - m.mean[j]) / m.scale[j];
            z[dim] = 1.0;
            double t = 0.0;
            for (std::size_t j = 0; j < k; ++j) t += theta[j] * z[j];
            const double p = sigmoid(t), s = std::max(p * (1.0 - p), 1e-12);
            for (std::size_t a = 0; a < k; ++a) {
                grad[a] += (double(y[i]) - p) * z[a];
                for (std::size_t b = 0; b <= a; ++b) hess[a * k + b] += s * z[a] * z[b];
            }
        }
        for (std::size_t a = 0; a < k; ++a)
            for (std::size_t b = 0; b < a; ++b) hess[b * k + a] = hess[a * k + b];
        for (std::size_t j = 0; j < dim; ++j) {
            hess[j * k + j] += l2 * double(n);
            grad[j] -= l2 * double(n) * theta[j];
        }
        if (!solve_spd(hess, grad, k)) break;
        double step = 0.0;
        for (std::size_t j = 0; j < k; ++j) theta[j] += grad[j], step = std::max(step, std::abs(grad[j]));
        if (step < 1e-10) break;
    }
    m.weights.assign(theta.begin(), theta.begin() + static_cast<std::ptrdiff_t>(dim));
    m.bias = theta[dim];
    return m;
}

LogisticModel parse_logistic_model(std::string_view json, std::string& error) {
    LogisticModel m;
    const JsonResult r = parse_json(json);
    if (!r.ok()) {
        error = r.error;
        return {};
    }
    auto nums = [&](const char* key, std::vector<double>& out) {
        const JsonValue* v = r.value.find(key);
        if (!v) return false;
        for (const JsonValue& e : v->array) out.push_back(e.number);
        return true;
    };
    const JsonValue* bias = r.value.find("bias");
    if (!nums("mean", m.mean) || !nums("scale", m.scale) || !nums("weights", m.weights) || !bias ||
        m.mean.size() != m.weights.size() || m.scale.size() != m.weights.size()) {
        error = "incomplete or inconsistent logistic model";
        return {};
    }
    m.bias = bias->number;
    if (const JsonValue* f = r.value.find("features")) {
        for (const JsonValue& e : f->array) m.features.push_back(e.string);
    }
    return m;
}

double roc_auc(std::span<const double> score, std::span<const int> label) {
    // Mann-Whitney U: rank all scores (average ranks on ties), sum the positives' ranks.
    std::vector<std::size_t> idx(score.size());
    std::iota(idx.begin(), idx.end(), std::size_t{0});
    std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) { return score[a] < score[b]; });
    std::vector<double> rank(score.size());
    for (std::size_t i = 0; i < idx.size();) {
        std::size_t j = i;
        while (j + 1 < idx.size() && score[idx[j + 1]] == score[idx[i]]) ++j;
        const double avg = 0.5 * double(i + j) + 1.0;
        for (std::size_t t = i; t <= j; ++t) rank[idx[t]] = avg;
        i = j + 1;
    }
    double pos = 0.0, neg = 0.0, sum = 0.0;
    for (std::size_t i = 0; i < label.size(); ++i) {
        if (label[i]) pos += 1.0, sum += rank[i];
        else neg += 1.0;
    }
    if (pos == 0.0 || neg == 0.0) return 0.5;
    return (sum - pos * (pos + 1.0) / 2.0) / (pos * neg);
}

}  // namespace dmw
