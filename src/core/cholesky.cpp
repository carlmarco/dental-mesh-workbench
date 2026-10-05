#include "core/cholesky.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace dmw {

std::vector<std::uint32_t> reverse_cuthill_mckee(const SparseMatrix& a) {
    const std::uint32_t n = a.rows;
    std::vector<std::uint32_t> degree(n);
    for (std::uint32_t r = 0; r < n; ++r) degree[r] = a.row_start[r + 1] - a.row_start[r];
    std::vector<std::uint32_t> by_degree(n);  // start candidates: lowest degree first
    for (std::uint32_t i = 0; i < n; ++i) by_degree[i] = i;
    std::stable_sort(by_degree.begin(), by_degree.end(),
                     [&](std::uint32_t x, std::uint32_t y) { return degree[x] < degree[y]; });

    std::vector<std::uint32_t> order;
    order.reserve(n);
    std::vector<bool> placed(n, false);
    std::vector<std::uint32_t> neighbors;
    for (std::uint32_t start : by_degree) {  // one breadth-first sweep per component
        if (placed[start]) continue;
        placed[start] = true;
        std::size_t head = order.size();  // `order` doubles as the breadth-first queue
        order.push_back(start);
        while (head < order.size()) {
            const std::uint32_t v = order[head++];
            neighbors.clear();
            for (std::uint32_t k = a.row_start[v]; k < a.row_start[v + 1]; ++k) {
                const std::uint32_t w = a.col[k];
                if (!placed[w]) placed[w] = true, neighbors.push_back(w);
            }
            std::sort(neighbors.begin(), neighbors.end(),
                      [&](std::uint32_t x, std::uint32_t y) { return degree[x] < degree[y]; });
            order.insert(order.end(), neighbors.begin(), neighbors.end());
        }
    }
    std::reverse(order.begin(), order.end());  // the reversal shrinks the envelope further
    return order;
}

EnvelopeLdlt::EnvelopeLdlt(const SparseMatrix& a) : order_(reverse_cuthill_mckee(a)) {
    const std::uint32_t n = a.rows;
    std::vector<std::uint32_t> pos(n);  // pos[old] = new
    for (std::uint32_t k = 0; k < n; ++k) pos[order_[k]] = k;

    // Envelope of the permuted matrix B = P A P^T: row i spans columns first[i] .. i.
    first_.resize(n);
    for (std::uint32_t i = 0; i < n; ++i) {
        std::uint32_t f = i;
        const std::uint32_t old = order_[i];
        for (std::uint32_t k = a.row_start[old]; k < a.row_start[old + 1]; ++k) f = std::min(f, pos[a.col[k]]);
        first_[i] = f;
    }
    row_start_.resize(n + 1, 0);
    for (std::uint32_t i = 0; i < n; ++i) row_start_[i + 1] = row_start_[i] + (i - first_[i]);
    l_.assign(row_start_[n], 0.0);
    d_.assign(n, 0.0);
    // Scatter B into the envelope (lower triangle and diagonal).
    for (std::uint32_t i = 0; i < n; ++i) {
        const std::uint32_t old = order_[i];
        for (std::uint32_t k = a.row_start[old]; k < a.row_start[old + 1]; ++k) {
            const std::uint32_t j = pos[a.col[k]];
            if (j == i) d_[i] = a.value[k];
            else if (j < i) l_[row_start_[i] + (j - first_[i])] = a.value[k];
        }
    }

    // Row-by-row LDL^T. Row i of B = L D L^T reads, for each j < i (L_jj = 1):
    //   B_ij = sum_{k<j} L_ik d_k L_jk + L_ij d_j.
    // With g_k = L_ik d_k this gives  g_j = B_ij - sum_{k<j} g_k L_jk,  L_ij = g_j / d_j,
    // and d_i = B_ii - sum_{j<i} g_j L_ij. Rows j < i are finished (they hold L_jk).
    for (std::uint32_t i = 0; i < n; ++i) {
        double* li = l_.data() + row_start_[i];
        const std::uint32_t fi = first_[i];
        double di = d_[i];
        // Pass 1: slot j goes from B_ij to g_j. The inner sum needs g_k for k < j, which
        // is exactly what the earlier slots of this row hold at that moment.
        for (std::uint32_t j = fi; j < i; ++j) {
            const double* lj = l_.data() + row_start_[j];
            const std::uint32_t fj = first_[j];
            double g = li[j - fi];
            for (std::uint32_t k = std::max(fi, fj); k < j; ++k) g -= li[k - fi] * lj[k - fj];
            li[j - fi] = g;
        }
        // Pass 2: slot j goes from g_j to L_ij, and d_i accumulates.
        for (std::uint32_t j = fi; j < i; ++j) {
            const double g = li[j - fi];
            const double lij = g / d_[j];
            di -= g * lij;
            li[j - fi] = lij;
        }
        if (!(di > 0.0)) {
            error_ = "not positive definite at pivot " + std::to_string(i);
            return;
        }
        d_[i] = di;
    }
}

std::vector<double> EnvelopeLdlt::solve(std::span<const double> b) const {
    const std::size_t n = d_.size();
    std::vector<double> y(n);
    for (std::size_t i = 0; i < n; ++i) y[i] = b[order_[i]];  // permute: y = P b
    for (std::size_t i = 0; i < n; ++i) {                     // forward: L z = y (unit lower)
        const double* li = l_.data() + row_start_[i];
        double s = y[i];
        for (std::uint32_t k = first_[i]; k < i; ++k) s -= li[k - first_[i]] * y[k];
        y[i] = s;
    }
    for (std::size_t i = 0; i < n; ++i) y[i] /= d_[i];  // diagonal: D w = z
    for (std::size_t i = n; i-- > 0;) {                // backward: L^T x = w, column-wise
        const double* li = l_.data() + row_start_[i];
        for (std::uint32_t k = first_[i]; k < i; ++k) y[k] -= li[k - first_[i]] * y[i];
    }
    std::vector<double> x(n);
    for (std::size_t i = 0; i < n; ++i) x[order_[i]] = y[i];  // un-permute
    return x;
}

}  // namespace dmw
