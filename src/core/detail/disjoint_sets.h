#pragma once
// Private to dmw_core.

#include <cstdint>
#include <numeric>
#include <utility>
#include <vector>

namespace dmw::detail {

// Union-find over {0, ..., n-1}. Union by size + path halving: any sequence of m operations
// costs O(m * alpha(n)), where alpha (inverse Ackermann) is <= 4 for any practical n.
class DisjointSets {
public:
    explicit DisjointSets(std::size_t n = 0) { reset(n); }

    // Re-initializes to n singletons, reusing storage.
    void reset(std::size_t n) {
        parent_.resize(n);
        std::iota(parent_.begin(), parent_.end(), std::uint32_t{0});
        size_.assign(n, 1);
    }

    std::uint32_t find(std::uint32_t x) {
        // Path halving: point every other node on the path at its grandparent.
        while (parent_[x] != x) {
            parent_[x] = parent_[parent_[x]];
            x = parent_[x];
        }
        return x;
    }

    // Returns true if a and b were in different sets.
    bool unite(std::uint32_t a, std::uint32_t b) {
        a = find(a);
        b = find(b);
        if (a == b) return false;
        if (size_[a] < size_[b]) std::swap(a, b);  // attach the smaller tree under the larger
        parent_[b] = a;
        size_[a] += size_[b];
        return true;
    }

private:
    std::vector<std::uint32_t> parent_;
    std::vector<std::uint32_t> size_;
};

}  // namespace dmw::detail
