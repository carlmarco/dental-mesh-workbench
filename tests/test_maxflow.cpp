#include <algorithm>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/maxflow.h"

using namespace dmw;
using Catch::Matchers::WithinAbs;

namespace {

struct Problem {
    std::uint32_t n = 0;
    std::vector<double> source, sink;  // per node
    struct E {
        std::uint32_t u, v;
        double cap, rev;
    };
    std::vector<E> edges;
};

// Value of the cut that puts node set {v : in_source[v]} with the source.
double cut_value(const Problem& p, const std::vector<bool>& in_source) {
    double c = 0.0;
    for (std::uint32_t v = 0; v < p.n; ++v) c += in_source[v] ? p.sink[v] : p.source[v];
    for (const auto& e : p.edges) {
        if (in_source[e.u] && !in_source[e.v]) c += e.cap;
        if (in_source[e.v] && !in_source[e.u]) c += e.rev;
    }
    return c;
}

double brute_force_min_cut(const Problem& p) {
    double best = 1e300;
    for (std::uint32_t mask = 0; mask < (1u << p.n); ++mask) {
        std::vector<bool> s(p.n);
        for (std::uint32_t v = 0; v < p.n; ++v) s[v] = (mask >> v) & 1u;
        best = std::min(best, cut_value(p, s));
    }
    return best;
}

MaxFlow build(const Problem& p) {
    MaxFlow g(p.n);
    for (std::uint32_t v = 0; v < p.n; ++v) g.add_terminal(v, p.source[v], p.sink[v]);
    for (const auto& e : p.edges) g.add_edge(e.u, e.v, e.cap, e.rev);
    return g;
}

// Independent reference: Edmonds-Karp (BFS augmenting paths) on a dense capacity matrix with
// super-source n and super-sink n + 1.
double edmonds_karp(const Problem& p) {
    const std::uint32_t n = p.n + 2, s = p.n, t = p.n + 1;
    std::vector<std::vector<double>> c(n, std::vector<double>(n, 0.0));
    for (std::uint32_t v = 0; v < p.n; ++v) c[s][v] += p.source[v], c[v][t] += p.sink[v];
    for (const auto& e : p.edges) c[e.u][e.v] += e.cap, c[e.v][e.u] += e.rev;
    double flow = 0.0;
    for (;;) {
        std::vector<std::uint32_t> prev(n, n);
        std::vector<std::uint32_t> queue{s};
        prev[s] = s;
        for (std::size_t q = 0; q < queue.size() && prev[t] == n; ++q)
            for (std::uint32_t w = 0; w < n; ++w)
                if (prev[w] == n && c[queue[q]][w] > 1e-12) prev[w] = queue[q], queue.push_back(w);
        if (prev[t] == n) return flow;
        double b = std::numeric_limits<double>::infinity();
        for (std::uint32_t v = t; v != s; v = prev[v]) b = std::min(b, c[prev[v]][v]);
        for (std::uint32_t v = t; v != s; v = prev[v]) c[prev[v]][v] -= b, c[v][prev[v]] += b;
        flow += b;
    }
}

}  // namespace

TEST_CASE("maxflow: textbook network (CLRS Fig. 26.1, max flow 23)", "[maxflow]") {
    // s=0 and t=5 of the CLRS network become terminal capacities on v1..v4 (nodes 0..3).
    MaxFlow g(4);
    g.add_terminal(0, 16, 0);   // s -> v1
    g.add_terminal(1, 13, 0);   // s -> v2
    g.add_terminal(2, 0, 20);   // v3 -> t
    g.add_terminal(3, 0, 4);    // v4 -> t
    g.add_edge(0, 2, 12, 0);    // v1 -> v3
    g.add_edge(1, 0, 4, 0);     // v2 -> v1
    g.add_edge(1, 3, 14, 0);    // v2 -> v4
    g.add_edge(2, 1, 9, 0);     // v3 -> v2
    g.add_edge(3, 2, 7, 0);     // v4 -> v3
    CHECK_THAT(g.solve(), WithinAbs(23.0, 1e-12));
    // Minimal source set of the min cut {s, v1, v2, v4} | {v3, t}: 12 + 7 + 4 = 23.
    CHECK(g.source_side(0));
    CHECK(g.source_side(1));
    CHECK_FALSE(g.source_side(2));
    CHECK(g.source_side(3));
}

TEST_CASE("maxflow: matches brute-force min cut on random graphs, and the labelling achieves it", "[maxflow]") {
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> cap(0.0, 10.0);
    std::uniform_int_distribution<int> coin(0, 3);
    for (int trial = 0; trial < 300; ++trial) {
        Problem p;
        p.n = 2 + static_cast<std::uint32_t>(trial % 9);  // 2..10 nodes: 2^10 cuts at most
        for (std::uint32_t v = 0; v < p.n; ++v) {
            p.source.push_back(coin(rng) == 0 ? cap(rng) : 0.0);
            p.sink.push_back(coin(rng) == 0 ? cap(rng) : 0.0);
        }
        for (std::uint32_t u = 0; u < p.n; ++u)
            for (std::uint32_t v = u + 1; v < p.n; ++v)
                if (coin(rng) < 2) p.edges.push_back({u, v, coin(rng) ? cap(rng) : 0.0, coin(rng) ? cap(rng) : 0.0});
        MaxFlow g = build(p);
        const double flow = g.solve();
        const double best = brute_force_min_cut(p);
        INFO("trial " << trial);
        CHECK_THAT(flow, WithinAbs(best, 1e-9));
        std::vector<bool> s(p.n);
        for (std::uint32_t v = 0; v < p.n; ++v) s[v] = g.source_side(v);
        CHECK_THAT(cut_value(p, s), WithinAbs(best, 1e-9));  // max-flow = min-cut, certified by the labels
    }
}

TEST_CASE("maxflow: hard constraints and isolated nodes", "[maxflow]") {
    // A chain 0 - 1 - 2 - 3 with a weak link in the middle; ends are hard-tied to the terminals.
    MaxFlow g(5);  // node 4 is isolated
    g.add_terminal(0, 1e12, 0);
    g.add_terminal(3, 0, 1e12);
    g.add_edge(0, 1, 5, 5);
    g.add_edge(1, 2, 0.5, 0.5);
    g.add_edge(2, 3, 5, 5);
    CHECK_THAT(g.solve(), WithinAbs(0.5, 1e-12));
    CHECK(g.source_side(0));
    CHECK(g.source_side(1));
    CHECK_FALSE(g.source_side(2));
    CHECK_FALSE(g.source_side(3));
    CHECK_FALSE(g.source_side(4));  // unreachable from the source: sink side by convention
}

TEST_CASE("maxflow: a long grid (orphan adoption path) agrees with the analytic cut", "[maxflow]") {
    // W x H grid, source tied to the left column, sink to the right; uniform horizontal capacity c
    // and vertical capacity 1. The min cut severs one column of H horizontal edges: H * c.
    const std::uint32_t w = 40, h = 25;
    const double c = 0.7;
    MaxFlow g(w * h);
    auto id = [w](std::uint32_t x, std::uint32_t y) { return y * w + x; };
    for (std::uint32_t y = 0; y < h; ++y) {
        g.add_terminal(id(0, y), 1e9, 0);
        g.add_terminal(id(w - 1, y), 0, 1e9);
        for (std::uint32_t x = 0; x < w; ++x) {
            if (x + 1 < w) g.add_edge(id(x, y), id(x + 1, y), c, c);
            if (y + 1 < h) g.add_edge(id(x, y), id(x, y + 1), 1.0, 1.0);
        }
    }
    CHECK_THAT(g.solve(), WithinAbs(h * c, 1e-9));
}

TEST_CASE("maxflow: matches Edmonds-Karp on larger sparse random graphs (deep trees, many orphans)", "[maxflow]") {
    // Mesh-like sparsity: a ring lattice plus random chords, terminals on a few nodes only, so the
    // search trees get deep and augmentations orphan whole subtrees (exercises adoption and freeing).
    std::mt19937 rng(777);
    std::uniform_real_distribution<double> cap(0.0, 10.0);
    for (int trial = 0; trial < 600; ++trial) {
        Problem p;
        p.n = 6 + static_cast<std::uint32_t>(trial % 30);  // small enough that freed nodes matter (kills a missing re-activation)
        std::uniform_int_distribution<std::uint32_t> node(0, p.n - 1);
        p.source.assign(p.n, 0.0);
        p.sink.assign(p.n, 0.0);
        for (int k = 0; k < 3; ++k) p.source[node(rng)] += cap(rng) * 3, p.sink[node(rng)] += cap(rng) * 3;
        for (std::uint32_t v = 0; v < p.n; ++v) {
            p.edges.push_back({v, (v + 1) % p.n, cap(rng), cap(rng)});
            p.edges.push_back({v, (v + 2) % p.n, cap(rng), cap(rng)});
        }
        for (std::uint32_t k = 0; k < p.n / 2; ++k) {
            const std::uint32_t u = node(rng), v = node(rng);
            if (u != v) p.edges.push_back({u, v, cap(rng), cap(rng)});
        }
        MaxFlow g = build(p);
        const double flow = g.solve();
        INFO("trial " << trial << ", n " << p.n);
        CHECK_THAT(flow, WithinAbs(edmonds_karp(p), 1e-7));
        std::vector<bool> s(p.n);
        for (std::uint32_t v = 0; v < p.n; ++v) s[v] = g.source_side(v);
        CHECK_THAT(cut_value(p, s), WithinAbs(flow, 1e-7));
    }
}
