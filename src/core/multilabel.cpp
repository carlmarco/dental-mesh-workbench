#include "core/multilabel.h"

#include <algorithm>

#include "core/maxflow.h"

namespace dmw {

double potts_energy(std::uint32_t labels, std::span<const double> unary, std::span<const PottsEdge> edges,
                    std::span<const std::uint32_t> f) {
    double e = 0.0;
    for (std::size_t i = 0; i < f.size(); ++i) e += unary[i * labels + f[i]];
    for (const PottsEdge& p : edges) e += f[p.a] != f[p.b] ? p.w : 0.0;
    return e;
}

bool star_feasible(const StarParents& star_parent, std::span<const std::uint32_t> f) {
    for (std::size_t i = 0; i < f.size(); ++i) {
        if (f[i] >= star_parent.size() || star_parent[f[i]].empty()) continue;
        const std::uint32_t p = star_parent[f[i]][i];
        if (p == kNoParent || (p != i && f[p] != f[i])) return false;
    }
    return true;
}

void star_repair(const StarParents& star_parent, std::span<const std::uint32_t> order, std::uint32_t fallback,
                 std::vector<std::uint32_t>& f) {
    for (std::uint32_t i : order) {
        const std::uint32_t l = f[i];
        if (l >= star_parent.size() || star_parent[l].empty()) continue;
        const std::uint32_t p = star_parent[l][i];
        if (p == kNoParent || (p != i && f[p] != l)) f[i] = fallback;  // parents come first, so cascades resolve
    }
}

double alpha_expansion(std::uint32_t labels, std::span<const double> unary, std::span<const PottsEdge> edges,
                       std::vector<std::uint32_t>& f, int max_sweeps, const ExpansionCandidates& candidates,
                       const StarParents& star_parent) {
    constexpr double kInf = 1e15;  // finite stand-in for an infinite capacity (max-flow needs finite values)
    auto parents_of = [&](std::uint32_t l) -> const std::vector<std::uint32_t>* {
        return l < star_parent.size() && !star_parent[l].empty() ? &star_parent[l] : nullptr;
    };
    const auto n = static_cast<std::uint32_t>(f.size());
    constexpr std::uint32_t kOut = 0xFFFFFFFFu;
    std::vector<std::uint8_t> mask(n, 1);
    std::vector<std::uint32_t> local(n);  // node -> index in the move graph, kOut if not in it
    std::vector<double> cost0, cost1;
    for (int sweep = 0; sweep < max_sweeps; ++sweep) {
        bool changed = false;
        for (std::uint32_t alpha = 0; alpha < labels; ++alpha) {
            // Move variable x_i: 0 = keep f_i (source side), 1 = take alpha (sink side). Only nodes not yet
            // labelled alpha (and passing the filter) are variables; the rest are fixed at their label.
            if (candidates) candidates(alpha, f, mask);
            std::uint32_t m = 0;
            for (std::uint32_t i = 0; i < n; ++i) local[i] = (f[i] != alpha && mask[i]) ? m++ : kOut;
            if (m == 0) continue;
            MaxFlow move(m);
            cost0.assign(m, 0.0), cost1.assign(m, 0.0);
            for (std::uint32_t i = 0; i < n; ++i) {
                if (local[i] == kOut) continue;
                cost0[local[i]] = unary[static_cast<std::size_t>(i) * labels + f[i]];
                cost1[local[i]] = unary[static_cast<std::size_t>(i) * labels + alpha];
            }
            // E(x_a, x_b) with A = E00, B = E01, C = E10, D = E11 (= 0 here: both alpha) equals
            // A + (C - A) x_a + (D - C) x_b + (B + C - A - D)(1 - x_a) x_b. The constant A is dropped, the
            // linear terms go to the unaries, the last term is an arc a -> b (cut iff a keeps and b takes alpha).
            // A fixed endpoint has x = 0 (it keeps its label; for a node labelled alpha both choices agree).
            for (const PottsEdge& e : edges) {
                const std::uint32_t ia = local[e.a], ib = local[e.b];
                if (ia == kOut && ib == kOut) continue;
                const std::uint32_t la = f[e.a], lb = f[e.b];
                const double A = la != lb ? e.w : 0.0, B = la != alpha ? e.w : 0.0, C = alpha != lb ? e.w : 0.0;
                if (ia != kOut && ib != kOut) {
                    cost1[ia] += C - A;
                    cost1[ib] -= C;
                    const double c = B + C - A;  // >= 0 by the triangle inequality of Potts
                    if (c > 0.0) move.add_edge(ia, ib, c, 0.0);
                } else if (ia != kOut) {
                    cost1[ia] += C - A;  // x_b = 0: E = A + (C - A) x_a
                } else {
                    cost1[ib] += B - A;  // x_a = 0: E = A + (B - A) x_b
                }
            }
            // Star constraints (D87). (1) Taking alpha: node i may switch only if its alpha-parent ends up alpha.
            if (const auto* pa = parents_of(alpha)) {
                for (std::uint32_t i = 0; i < n; ++i) {
                    if (local[i] == kOut) continue;
                    const std::uint32_t p = (*pa)[i];
                    if (p == i) continue;                                          // a root
                    if (p == kNoParent) { cost1[local[i]] += kInf; continue; }     // may not take alpha
                    if (f[p] == alpha) continue;                                   // parent already alpha (stays)
                    if (local[p] == kOut) { cost1[local[i]] += kInf; continue; }   // parent cannot become alpha
                    move.add_edge(local[p], local[i], kInf, 0.0);                  // cut iff p keeps and i takes alpha
                }
            }
            // (2) Keeping l: a node that keeps label l forbids its l-parent from switching to alpha.
            for (std::uint32_t i = 0; i < n; ++i) {
                const std::uint32_t l = f[i];
                if (l == alpha) continue;
                const auto* pl = parents_of(l);
                if (!pl) continue;
                const std::uint32_t p = (*pl)[i];
                if (p == i || p == kNoParent || local[p] == kOut) continue;        // root, or the parent cannot move
                if (local[i] == kOut) cost1[local[p]] += kInf;                     // i is fixed: p may not switch
                else move.add_edge(local[i], local[p], kInf, 0.0);                 // cut iff i keeps and p switches
            }
            for (std::uint32_t i = 0; i < m; ++i) {
                const double lo = std::min(cost0[i], cost1[i]);
                move.add_terminal(i, cost1[i] - lo, cost0[i] - lo);  // source cap = cost of x = 1, sink cap = cost of x = 0
            }
            move.solve();
            for (std::uint32_t i = 0; i < n; ++i)
                if (local[i] != kOut && !move.source_side(local[i])) f[i] = alpha, changed = true;
        }
        if (!changed) break;
    }
    return potts_energy(labels, unary, edges, f);
}

}  // namespace dmw
