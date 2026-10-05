#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/cholesky.h"
#include "core/dec.h"
#include "core/generate.h"

using namespace dmw;
using Catch::Matchers::WithinAbs;

namespace {

// Path graph Laplacian with Dirichlet ends, nodes labelled in a scrambled order, so the
// natural numbering has a large bandwidth that RCM must undo.
SparseMatrix scrambled_path(std::uint32_t n, std::vector<std::uint32_t>& label) {
    label.resize(n);
    for (std::uint32_t i = 0; i < n; ++i) label[i] = (i * 37) % n;  // 37 coprime to n
    std::vector<Triplet> t;
    for (std::uint32_t i = 0; i < n; ++i) {
        t.push_back({label[i], label[i], 2.0});
        if (i + 1 < n) {
            t.push_back({label[i], label[i + 1], -1.0});
            t.push_back({label[i + 1], label[i], -1.0});
        }
    }
    return SparseMatrix::from_triplets(n, n, std::move(t));
}

std::uint32_t bandwidth(const SparseMatrix& a, const std::vector<std::uint32_t>& order) {
    std::vector<std::uint32_t> pos(order.size());
    for (std::uint32_t k = 0; k < order.size(); ++k) pos[order[k]] = k;
    std::uint32_t bw = 0;
    for (std::uint32_t r = 0; r < a.rows; ++r) {
        for (std::uint32_t k = a.row_start[r]; k < a.row_start[r + 1]; ++k) {
            const std::uint32_t i = pos[r], j = pos[a.col[k]];
            bw = std::max(bw, i > j ? i - j : j - i);
        }
    }
    return bw;
}

double max_residual(const SparseMatrix& a, const std::vector<double>& x, const std::vector<double>& b) {
    const auto ax = a.multiply(x);
    double r = 0.0;
    for (std::size_t i = 0; i < b.size(); ++i) r = std::max(r, std::abs(ax[i] - b[i]));
    return r;
}

}  // namespace

TEST_CASE("rcm: returns a permutation and restores a scrambled path to bandwidth 1", "[cholesky]") {
    std::vector<std::uint32_t> label;
    const auto a = scrambled_path(101, label);
    std::vector<std::uint32_t> identity(101);
    std::iota(identity.begin(), identity.end(), 0u);
    CHECK(bandwidth(a, identity) > 30);
    const auto order = reverse_cuthill_mckee(a);
    auto sorted = order;
    std::sort(sorted.begin(), sorted.end());
    CHECK(sorted == identity);  // each old index exactly once
    CHECK(bandwidth(a, order) == 1);
}

TEST_CASE("ldlt: solves an SPD system to roundoff", "[cholesky]") {
    std::vector<std::uint32_t> label;
    const auto a = scrambled_path(300, label);
    std::vector<double> b(300);
    for (std::size_t i = 0; i < b.size(); ++i) b[i] = std::cos(0.1 * static_cast<double>(i));
    const EnvelopeLdlt f(a);
    REQUIRE(f.ok());
    CHECK(max_residual(a, f.solve(b), b) < 1e-10);
}

TEST_CASE("ldlt: heat matrix *0 - t L of a mesh, and its solution is positive everywhere", "[cholesky]") {
    // The M6 regression (D47): at t = h^2 the heat solution far from the source is tiny but
    // positive. A direct solve keeps every entry positive; CG at 1e-10 zeroed a third of them.
    auto r = build_halfedge(make_icosphere(4));
    REQUIRE(r.ok());
    const auto ops = build_dec(r.mesh);
    const double t = ops.mean_edge_length * ops.mean_edge_length;
    const auto heat = scaled_plus_diagonal(ops.laplacian, -t, ops.star0);
    const EnvelopeLdlt f(heat);
    REQUIRE(f.ok());
    std::vector<double> delta(heat.rows, 0.0);
    delta[0] = 1.0;
    const auto u = f.solve(delta);
    CHECK(max_residual(heat, u, delta) < 1e-12);
    CHECK(*std::min_element(u.begin(), u.end()) > 0.0);
}

TEST_CASE("ldlt: rejects a matrix that is not positive definite", "[cholesky]") {
    const auto a = SparseMatrix::from_triplets(2, 2, {{0, 0, 1.0}, {0, 1, 2.0}, {1, 0, 2.0}, {1, 1, 1.0}});
    const EnvelopeLdlt f(a);  // eigenvalues 3 and -1
    CHECK_FALSE(f.ok());
}
