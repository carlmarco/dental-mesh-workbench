#pragma once

#include <cstdint>
#include <span>
#include <vector>

namespace dmw {

struct Triplet {
    std::uint32_t row;
    std::uint32_t col;
    double value;
};

// Compressed sparse row (CSR) matrix: the entries of row r are
// (col[k], value[k]) for k in [row_start[r], row_start[r + 1]), columns ascending.
struct SparseMatrix {
    std::uint32_t rows = 0;
    std::uint32_t cols = 0;
    std::vector<std::uint32_t> row_start;  // rows + 1 entries
    std::vector<std::uint32_t> col;
    std::vector<double> value;

    // Duplicate (row, col) entries are summed (the natural way to assemble FEM-style matrices).
    static SparseMatrix from_triplets(std::uint32_t rows, std::uint32_t cols, std::vector<Triplet> entries);

    void multiply(std::span<const double> x, std::span<double> y) const;  // y = A x
    std::vector<double> multiply(std::span<const double> x) const;
    double at(std::uint32_t r, std::uint32_t c) const;  // 0 if not stored (binary search)
    std::vector<double> diagonal() const;
    SparseMatrix transpose() const;
    std::size_t nonzeros() const { return value.size(); }
};

// D^T diag(w) D, with w.size() == D.rows. In DEC: L = -d0^T *1 d0 (D45).
SparseMatrix weighted_gram(const SparseMatrix& d, std::span<const double> w);

// s * A + diag(d), with d.size() == A.rows (A square).
SparseMatrix scaled_plus_diagonal(const SparseMatrix& a, double s, std::span<const double> d);

struct CgResult {
    std::vector<double> x;
    std::uint32_t iterations = 0;
    double relative_residual = 0.0;  // |b - A x| / |b|
    bool converged = false;
};

// Jacobi-preconditioned conjugate gradient for symmetric positive (semi)definite A (D46).
// For singular A the system must be consistent (b in the range of A), e.g. -L with b summing
// to zero on each connected component. max_iterations == 0 means 10 * rows.
CgResult conjugate_gradient(const SparseMatrix& a, std::span<const double> b, double tolerance = 1e-10,
                            std::uint32_t max_iterations = 0);

}  // namespace dmw
