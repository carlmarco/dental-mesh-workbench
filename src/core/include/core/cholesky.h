#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "core/sparse.h"

namespace dmw {

// Reverse Cuthill-McKee ordering of a symmetric matrix's adjacency graph: breadth-first from a
// low-degree vertex in each component, neighbors in increasing degree, then reversed. Keeps
// coupled unknowns close in index, so the matrix becomes banded. Returns order[new] = old.
std::vector<std::uint32_t> reverse_cuthill_mckee(const SparseMatrix& a);

// Sparse LDL^T factorization in envelope (profile/skyline) storage (D49): row i of the
// permuted matrix stores columns first[i] .. i-1. Cholesky fill-in never leaves the envelope,
// so the storage is fixed before factoring. Factor once, then solve() any number of times.
class EnvelopeLdlt {
public:
    // Factors a symmetric positive definite matrix (both triangles must be stored).
    explicit EnvelopeLdlt(const SparseMatrix& a);
    bool ok() const { return error_.empty(); }
    const std::string& error() const { return error_; }  // "not positive definite at ..." etc.
    std::vector<double> solve(std::span<const double> b) const;
    std::size_t stored_entries() const { return l_.size(); }  // strictly-lower envelope size

private:
    std::vector<std::uint32_t> order_;    // order_[new] = old
    std::vector<std::uint32_t> first_;    // first stored column of each (permuted) row
    std::vector<std::size_t> row_start_;  // offset of row i's entries in l_
    std::vector<double> l_;               // unit-lower L, strictly below the diagonal
    std::vector<double> d_;               // diagonal D
    std::string error_;
};

}  // namespace dmw
