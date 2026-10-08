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

double alpha_expansion(std::uint32_t labels, std::span<const double> unary, std::span<const PottsEdge> edges,
                       std::vector<std::uint32_t>& f, int max_sweeps, const ExpansionCandidates& candidates) {
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
