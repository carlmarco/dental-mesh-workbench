#include <algorithm>
#include <cstdint>
#include <random>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "core/multilabel.h"

using namespace dmw;

namespace {

struct Problem {
    std::uint32_t n = 0, labels = 0;
    std::vector<double> unary;
    std::vector<PottsEdge> edges;
};

Problem random_problem(std::mt19937& rng, std::uint32_t n, std::uint32_t labels) {
    std::uniform_real_distribution<double> u(0.0, 10.0), w(0.0, 6.0);
    std::uniform_int_distribution<int> coin(0, 2);
    Problem p{n, labels, {}, {}};
    for (std::uint32_t i = 0; i < n * labels; ++i) p.unary.push_back(u(rng));
    for (std::uint32_t a = 0; a < n; ++a)
        for (std::uint32_t b = a + 1; b < n; ++b)
            if (coin(rng)) p.edges.push_back({a, b, w(rng)});
    return p;
}

// Exhaustive minimum over all labels^n labellings.
double brute_force(const Problem& p) {
    std::vector<std::uint32_t> f(p.n, 0);
    double best = 1e300;
    for (;;) {
        best = std::min(best, potts_energy(p.labels, p.unary, p.edges, f));
        std::uint32_t i = 0;
        while (i < p.n && ++f[i] == p.labels) f[i++] = 0;
        if (i == p.n) return best;
    }
}

// Best energy reachable from f by ONE expansion move on alpha (exhaustive over the 2^n move subsets).
double best_expansion(const Problem& p, const std::vector<std::uint32_t>& f, std::uint32_t alpha) {
    double best = 1e300;
    for (std::uint32_t mask = 0; mask < (1u << p.n); ++mask) {
        std::vector<std::uint32_t> g = f;
        for (std::uint32_t i = 0; i < p.n; ++i)
            if ((mask >> i) & 1u) g[i] = alpha;
        best = std::min(best, potts_energy(p.labels, p.unary, p.edges, g));
    }
    return best;
}

}  // namespace

TEST_CASE("multilabel: alpha-expansion never increases energy, reaches an expansion-move local minimum, within 2x optimum", "[multilabel]") {
    std::mt19937 rng(2024);
    for (int trial = 0; trial < 300; ++trial) {
        const std::uint32_t n = 2 + static_cast<std::uint32_t>(trial % 6), labels = 2 + static_cast<std::uint32_t>(trial % 3);
        const Problem p = random_problem(rng, n, labels);
        std::vector<std::uint32_t> f(n);
        for (auto& x : f) x = std::uniform_int_distribution<std::uint32_t>(0, labels - 1)(rng);
        const double e0 = potts_energy(labels, p.unary, p.edges, f);
        const double e = alpha_expansion(labels, p.unary, p.edges, f, 50);
        INFO("trial " << trial);
        CHECK(e <= e0 + 1e-9);
        for (std::uint32_t alpha = 0; alpha < labels; ++alpha) CHECK(best_expansion(p, f, alpha) >= e - 1e-9);  // local min
        CHECK(e <= 2.0 * brute_force(p) + 1e-9);  // BVZ bound for Potts
    }
}

TEST_CASE("multilabel: with two labels one expansion pass is the exact binary min cut", "[multilabel]") {
    std::mt19937 rng(99);
    for (int trial = 0; trial < 200; ++trial) {
        const Problem p = random_problem(rng, 2 + static_cast<std::uint32_t>(trial % 7), 2);
        std::vector<std::uint32_t> f(p.n, 0);
        const double e = alpha_expansion(2, p.unary, p.edges, f, 1);  // the alpha = 1 move from all-0 is the full binary problem
        INFO("trial " << trial);
        CHECK(e <= brute_force(p) + 1e-9);
    }
}

TEST_CASE("multilabel: an all-true candidate filter gives the same result; a restrictive one never raises energy", "[multilabel]") {
    std::mt19937 rng(7);
    for (int trial = 0; trial < 200; ++trial) {
        const std::uint32_t n = 3 + static_cast<std::uint32_t>(trial % 5), labels = 2 + static_cast<std::uint32_t>(trial % 3);
        const Problem p = random_problem(rng, n, labels);
        std::vector<std::uint32_t> init(n);
        for (auto& x : init) x = std::uniform_int_distribution<std::uint32_t>(0, labels - 1)(rng);
        std::vector<std::uint32_t> a = init, b = init, c = init;
        alpha_expansion(labels, p.unary, p.edges, a, 50);
        alpha_expansion(labels, p.unary, p.edges, b, 50,
                        [](std::uint32_t, const std::vector<std::uint32_t>&, std::vector<std::uint8_t>& mask) { std::fill(mask.begin(), mask.end(), 1); });
        INFO("trial " << trial);
        CHECK(a == b);
        // Only even nodes may move: still a descent method.
        const double e0 = potts_energy(labels, p.unary, p.edges, c);
        const double e = alpha_expansion(labels, p.unary, p.edges, c, 50,
                                         [](std::uint32_t, const std::vector<std::uint32_t>&, std::vector<std::uint8_t>& mask) {
                                             for (std::size_t i = 0; i < mask.size(); ++i) mask[i] = i % 2 == 0;
                                         });
        CHECK(e <= e0 + 1e-9);
        for (std::size_t i = 1; i < n; i += 2) CHECK(c[i] == init[i]);
    }
}

TEST_CASE("multilabel: star constraints hold throughout, and the result is a local minimum over feasible expansion moves", "[multilabel]") {
    // Each non-zero label gets a random tree over the nodes (a random root, parents chosen among earlier nodes in a
    // random order, some nodes forbidden). Label 0 is unconstrained (the fallback).
    std::mt19937 rng(31);
    for (int trial = 0; trial < 300; ++trial) {
        const std::uint32_t n = 3 + static_cast<std::uint32_t>(trial % 5), labels = 2 + static_cast<std::uint32_t>(trial % 3);
        const Problem p = random_problem(rng, n, labels);
        StarParents star(labels);
        std::vector<std::uint32_t> order(n);
        for (std::uint32_t i = 0; i < n; ++i) order[i] = i;
        std::shuffle(order.begin(), order.end(), rng);  // one global order: parents always earlier (valid for repair)
        for (std::uint32_t l = 1; l < labels; ++l) {
            star[l].assign(n, kNoParent);
            for (std::size_t k = 0; k < n; ++k) {
                const std::uint32_t i = order[k];
                if (k == 0) star[l][i] = i;  // root
                else if (std::uniform_int_distribution<int>(0, 5)(rng) == 0) star[l][i] = kNoParent;
                else star[l][i] = order[std::uniform_int_distribution<std::size_t>(0, k - 1)(rng)];
            }
        }
        std::vector<std::uint32_t> f(n);
        for (auto& x : f) x = std::uniform_int_distribution<std::uint32_t>(0, labels - 1)(rng);
        star_repair(star, order, 0, f);
        INFO("trial " << trial);
        REQUIRE(star_feasible(star, f));
        const double e0 = potts_energy(labels, p.unary, p.edges, f);
        const double e = alpha_expansion(labels, p.unary, p.edges, f, 50, {}, star);
        CHECK(star_feasible(star, f));
        CHECK(e <= e0 + 1e-9);
        // No feasible single expansion move improves the result (exhaustive over move subsets).
        for (std::uint32_t alpha = 0; alpha < labels; ++alpha) {
            for (std::uint32_t mask = 0; mask < (1u << n); ++mask) {
                std::vector<std::uint32_t> g = f;
                for (std::uint32_t i = 0; i < n; ++i)
                    if ((mask >> i) & 1u) g[i] = alpha;
                if (!star_feasible(star, g)) continue;
                CHECK(potts_energy(labels, p.unary, p.edges, g) >= e - 1e-9);
            }
        }
    }
}
