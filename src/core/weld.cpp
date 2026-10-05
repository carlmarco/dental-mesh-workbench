#include "core/weld.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <unordered_map>
#include <vector>

namespace dmw {
namespace {

constexpr std::uint32_t kUnassigned = std::numeric_limits<std::uint32_t>::max();

// boost-style hash_combine over three 64-bit words. The result is narrowed to size_t,
// which is 32 bits on wasm32 and 64 bits natively, hence the explicit cast.
template <typename Word>
struct Hash3 {
    std::size_t operator()(const std::array<Word, 3>& k) const noexcept {
        std::uint64_t h = 0;
        for (Word w : k) {
            std::uint64_t x = 0;
            std::memcpy(&x, &w, sizeof w);
            h ^= x + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        }
        return static_cast<std::size_t>(h);
    }
};

// Exact key: the coordinates' bit patterns, with -0.0 normalized to +0.0 so the two
// zeros (equal as numbers, different sign bit) hash and compare the same.
std::array<std::uint64_t, 3> exact_key(const Vec3& p) {
    std::array<std::uint64_t, 3> key{};
    const double c[3] = {p.x, p.y, p.z};
    for (std::size_t i = 0; i < 3; ++i) {
        const double v = (c[i] == 0.0) ? 0.0 : c[i];  // -0.0 == 0.0 is true, so this maps it to +0.0
        std::memcpy(&key[i], &v, sizeof v);
    }
    return key;
}

// Grid cell of coordinate v for cell size eps. Clamping keeps the double -> int64 cast
// defined for extreme v/eps; it can only put far-apart points in the same cell, which
// costs extra distance checks but never changes the result.
std::int64_t cell_of(double v, double eps) {
    constexpr double kLimit = 4.0e18;  // < 2^63
    return static_cast<std::int64_t>(std::clamp(std::floor(v / eps), -kLimit, kLimit));
}

double dist2(const Vec3& a, const Vec3& b) {
    const double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

}  // namespace

TriMesh weld_vertices(const TriMesh& soup, double epsilon) {
    TriMesh out;
    out.triangles.reserve(soup.triangles.size());
    // remap[old vertex] = new vertex; vertices are visited in triangle-corner order, so
    // output is numbered by first appearance and unreferenced input vertices are dropped.
    std::vector<std::uint32_t> remap(soup.positions.size(), kUnassigned);

    // Returns the output index for position p, creating a new representative if needed.
    std::unordered_map<std::array<std::uint64_t, 3>, std::uint32_t, Hash3<std::uint64_t>> exact;
    std::unordered_map<std::array<std::int64_t, 3>, std::vector<std::uint32_t>, Hash3<std::int64_t>> grid;
    const double eps2 = epsilon * epsilon;

    auto find_or_add = [&](const Vec3& p) -> std::uint32_t {
        const auto next = static_cast<std::uint32_t>(out.positions.size());
        if (epsilon == 0.0) {
            const auto [it, inserted] = exact.try_emplace(exact_key(p), next);
            if (inserted) out.positions.push_back(p);
            return it->second;
        }
        // |a - b| <= eps implies each coordinate differs by <= eps, so the cell indices
        // differ by at most 1 per axis: the 3x3x3 neighborhood contains every candidate.
        const std::array<std::int64_t, 3> c{cell_of(p.x, epsilon), cell_of(p.y, epsilon),
                                            cell_of(p.z, epsilon)};
        std::uint32_t best = kUnassigned;  // smallest index = earliest representative (D20)
        for (std::int64_t dx = -1; dx <= 1; ++dx) {
            for (std::int64_t dy = -1; dy <= 1; ++dy) {
                for (std::int64_t dz = -1; dz <= 1; ++dz) {
                    const auto it = grid.find({c[0] + dx, c[1] + dy, c[2] + dz});
                    if (it == grid.end()) continue;
                    for (std::uint32_t r : it->second) {
                        if (r < best && dist2(p, out.positions[r]) <= eps2) best = r;
                    }
                }
            }
        }
        if (best != kUnassigned) return best;
        grid[c].push_back(next);
        out.positions.push_back(p);
        return next;
    };

    for (const auto& tri : soup.triangles) {
        std::array<std::uint32_t, 3> welded{};
        for (std::size_t k = 0; k < 3; ++k) {
            std::uint32_t& slot = remap[tri[k]];
            if (slot == kUnassigned) slot = find_or_add(soup.positions[tri[k]]);
            welded[k] = slot;
        }
        out.triangles.push_back(welded);  // kept even if collapsed (D21)
    }
    return out;
}

}  // namespace dmw
