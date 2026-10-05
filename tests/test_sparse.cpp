#include <cmath>
#include <cstdint>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/sparse.h"

using namespace dmw;
using Catch::Matchers::WithinAbs;

namespace {

// 1D Laplacian on a path of n nodes: 2 on the diagonal, -1 off-diagonal (Dirichlet ends):
// symmetric positive definite.
SparseMatrix path_laplacian(std::uint32_t n, bool neumann) {
    std::vector<Triplet> t;
    for (std::uint32_t i = 0; i < n; ++i) {
        const double deg = neumann ? double((i > 0) + (i + 1 < n)) : 2.0;
        t.push_back({i, i, deg});
        if (i + 1 < n) {
            t.push_back({i, i + 1, -1.0});
            t.push_back({i + 1, i, -1.0});
        }
    }
    return SparseMatrix::from_triplets(n, n, std::move(t));
}

double residual(const SparseMatrix& a, const std::vector<double>& x, const std::vector<double>& b) {
    const auto ax = a.multiply(x);
    double r = 0.0;
    for (std::size_t i = 0; i < b.size(); ++i) r = std::max(r, std::abs(ax[i] - b[i]));
    return r;
}

}  // namespace

TEST_CASE("sparse: triplets are sorted and duplicates summed", "[sparse]") {
    const auto a = SparseMatrix::from_triplets(2, 3, {{1, 2, 4.0}, {0, 1, 1.0}, {1, 0, 2.0}, {0, 1, 0.5}});
    CHECK(a.nonzeros() == 3);
    CHECK(a.at(0, 1) == 1.5);
    CHECK(a.at(1, 0) == 2.0);
    CHECK(a.at(1, 2) == 4.0);
    CHECK(a.at(0, 0) == 0.0);
    CHECK(a.row_start == std::vector<std::uint32_t>{0, 1, 3});
    CHECK(a.col == std::vector<std::uint32_t>{1, 0, 2});  // ascending within each row
}

TEST_CASE("sparse: multiply and transpose match the dense definitions", "[sparse]") {
    // [[1, 0, 2], [0, 3, 0]]
    const auto a = SparseMatrix::from_triplets(2, 3, {{0, 0, 1}, {0, 2, 2}, {1, 1, 3}});
    CHECK(a.multiply(std::vector<double>{1, 2, 3}) == std::vector<double>{7, 6});
    const auto at = a.transpose();
    CHECK(at.rows == 3);
    CHECK(at.at(2, 0) == 2.0);
    CHECK(at.multiply(std::vector<double>{1, 1}) == std::vector<double>{1, 3, 2});
}

TEST_CASE("sparse: weighted Gram D^T W D matches the dense product", "[sparse]") {
    // D = [[-1, 1, 0], [0, -1, 1]] (incidence of a path), W = diag(2, 5)
    // D^T W D = [[2, -2, 0], [-2, 7, -5], [0, -5, 5]]
    const auto d = SparseMatrix::from_triplets(2, 3, {{0, 0, -1}, {0, 1, 1}, {1, 1, -1}, {1, 2, 1}});
    const auto g = weighted_gram(d, std::vector<double>{2, 5});
    const double expected[3][3] = {{2, -2, 0}, {-2, 7, -5}, {0, -5, 5}};
    for (std::uint32_t r = 0; r < 3; ++r) {
        for (std::uint32_t c = 0; c < 3; ++c) CHECK(g.at(r, c) == expected[r][c]);
    }
}

TEST_CASE("sparse: scaled_plus_diagonal", "[sparse]") {
    const auto a = path_laplacian(3, false);
    const auto m = scaled_plus_diagonal(a, 0.5, std::vector<double>{10, 20, 30});
    CHECK(m.at(0, 0) == 11.0);
    CHECK(m.at(1, 1) == 21.0);
    CHECK(m.at(1, 2) == -0.5);
}

TEST_CASE("cg: solves an SPD system to tolerance", "[sparse][cg]") {
    const std::uint32_t n = 200;
    const auto a = path_laplacian(n, false);
    std::vector<double> x_true(n), b;
    for (std::uint32_t i = 0; i < n; ++i) x_true[i] = std::sin(0.05 * i) + 0.01 * i;
    b = a.multiply(x_true);
    const auto r = conjugate_gradient(a, b, 1e-12);
    REQUIRE(r.converged);
    CHECK(r.iterations <= n);  // exact arithmetic: at most n steps
    for (std::uint32_t i = 0; i < n; ++i) CHECK_THAT(r.x[i], WithinAbs(x_true[i], 1e-8));
}

TEST_CASE("cg: singular but consistent system (Neumann Laplacian, zero-sum rhs)", "[sparse][cg]") {
    // Kernel = constants. b sums to zero, so a solution exists (unique up to a constant).
    const std::uint32_t n = 50;
    const auto a = path_laplacian(n, true);
    std::vector<double> b(n, 0.0);
    b[3] = 1.0;
    b[40] = -1.0;
    const auto r = conjugate_gradient(a, b, 1e-12);
    REQUIRE(r.converged);
    CHECK(residual(a, r.x, b) < 1e-9);
}

TEST_CASE("cg: zero right-hand side gives zero immediately", "[sparse][cg]") {
    const auto r = conjugate_gradient(path_laplacian(10, false), std::vector<double>(10, 0.0));
    CHECK(r.converged);
    CHECK(r.iterations == 0);
    CHECK(r.x == std::vector<double>(10, 0.0));
}
