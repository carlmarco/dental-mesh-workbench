#include "core/sparse.h"

#include <algorithm>
#include <cmath>
#include <tuple>

namespace dmw {
namespace {

double dot(std::span<const double> a, std::span<const double> b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}

}  // namespace

SparseMatrix SparseMatrix::from_triplets(std::uint32_t rows, std::uint32_t cols, std::vector<Triplet> entries) {
    std::sort(entries.begin(), entries.end(), [](const Triplet& a, const Triplet& b) {
        return std::tie(a.row, a.col) < std::tie(b.row, b.col);
    });
    SparseMatrix m;
    m.rows = rows;
    m.cols = cols;
    m.row_start.assign(rows + 1, 0);
    const Triplet* prev = nullptr;
    for (const Triplet& t : entries) {
        if (prev && prev->row == t.row && prev->col == t.col) {
            m.value.back() += t.value;  // sorted, so duplicates are adjacent: sum them
        } else {
            m.col.push_back(t.col);
            m.value.push_back(t.value);
            ++m.row_start[t.row + 1];  // count per row for now; prefix-summed below
        }
        prev = &t;
    }
    for (std::uint32_t r = 0; r < rows; ++r) m.row_start[r + 1] += m.row_start[r];  // counts -> offsets
    return m;
}

void SparseMatrix::multiply(std::span<const double> x, std::span<double> y) const {
    for (std::uint32_t r = 0; r < rows; ++r) {
        double s = 0.0;
        for (std::uint32_t k = row_start[r]; k < row_start[r + 1]; ++k) s += value[k] * x[col[k]];
        y[r] = s;
    }
}

std::vector<double> SparseMatrix::multiply(std::span<const double> x) const {
    std::vector<double> y(rows);
    multiply(x, y);
    return y;
}

double SparseMatrix::at(std::uint32_t r, std::uint32_t c) const {
    // Raw pointers: an unsigned offset needs no conversion. (Iterator + uint32_t converts to
    // the signed difference_type, which is 32-bit on wasm32 and trips -Wsign-conversion.)
    const std::uint32_t* begin = col.data() + row_start[r];
    const std::uint32_t* end = col.data() + row_start[r + 1];
    const std::uint32_t* it = std::lower_bound(begin, end, c);
    return (it != end && *it == c) ? value[static_cast<std::size_t>(it - col.data())] : 0.0;
}

std::vector<double> SparseMatrix::diagonal() const {
    std::vector<double> d(std::min(rows, cols), 0.0);
    for (std::uint32_t r = 0; r < d.size(); ++r) d[r] = at(r, r);
    return d;
}

SparseMatrix SparseMatrix::transpose() const {
    std::vector<Triplet> t;
    t.reserve(value.size());
    for (std::uint32_t r = 0; r < rows; ++r) {
        for (std::uint32_t k = row_start[r]; k < row_start[r + 1]; ++k) t.push_back({col[k], r, value[k]});
    }
    return from_triplets(cols, rows, std::move(t));
}

SparseMatrix weighted_gram(const SparseMatrix& d, std::span<const double> w) {
    // (D^T W D)_{ij} = sum_r w_r D_ri D_rj: every pair of nonzeros in a row of D contributes.
    // Rows of d0 have exactly two nonzeros (an edge's endpoints), so this is 4 triplets per edge.
    std::vector<Triplet> t;
    for (std::uint32_t r = 0; r < d.rows; ++r) {
        for (std::uint32_t a = d.row_start[r]; a < d.row_start[r + 1]; ++a) {
            for (std::uint32_t b = d.row_start[r]; b < d.row_start[r + 1]; ++b) {
                t.push_back({d.col[a], d.col[b], w[r] * d.value[a] * d.value[b]});
            }
        }
    }
    return SparseMatrix::from_triplets(d.cols, d.cols, std::move(t));
}

SparseMatrix scaled_plus_diagonal(const SparseMatrix& a, double s, std::span<const double> d) {
    std::vector<Triplet> t;
    t.reserve(a.value.size() + d.size());
    for (std::uint32_t r = 0; r < a.rows; ++r) {
        for (std::uint32_t k = a.row_start[r]; k < a.row_start[r + 1]; ++k) t.push_back({r, a.col[k], s * a.value[k]});
        t.push_back({r, r, d[r]});
    }
    return SparseMatrix::from_triplets(a.rows, a.cols, std::move(t));
}

CgResult conjugate_gradient(const SparseMatrix& a, std::span<const double> b, double tolerance,
                            std::uint32_t max_iterations) {
    const std::size_t n = a.rows;
    if (max_iterations == 0) max_iterations = static_cast<std::uint32_t>(10 * n);
    CgResult res;
    res.x.assign(n, 0.0);
    const double b_norm = std::sqrt(dot(b, b));
    if (b_norm == 0.0) {
        res.converged = true;
        return res;
    }
    // Jacobi preconditioner M = diag(A): cheap, and it equalizes the scaling of rows (vertex
    // areas vary across a mesh). Non-positive diagonals (isolated vertices) fall back to 1.
    std::vector<double> inv_diag = a.diagonal();
    for (double& d : inv_diag) d = d > 0.0 ? 1.0 / d : 1.0;

    std::vector<double> r(b.begin(), b.end()), z(n), p(n), ap(n);
    for (std::size_t i = 0; i < n; ++i) z[i] = inv_diag[i] * r[i];
    p = z;
    double rz = dot(r, z);
    for (res.iterations = 0; res.iterations < max_iterations; ++res.iterations) {
        a.multiply(p, ap);
        const double p_ap = dot(p, ap);
        if (p_ap <= 0.0) break;  // direction in the kernel (semidefinite A): nothing left to reduce
        const double alpha = rz / p_ap;
        for (std::size_t i = 0; i < n; ++i) {
            res.x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        res.relative_residual = std::sqrt(dot(r, r)) / b_norm;
        if (res.relative_residual < tolerance) {
            res.converged = true;
            ++res.iterations;
            break;
        }
        for (std::size_t i = 0; i < n; ++i) z[i] = inv_diag[i] * r[i];
        const double rz_next = dot(r, z);
        const double beta = rz_next / rz;  // keeps the new direction A-conjugate to the old ones
        rz = rz_next;
        for (std::size_t i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
    }
    return res;
}

}  // namespace dmw
