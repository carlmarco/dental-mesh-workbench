#pragma once
// Private to dmw_core.

#include <array>
#include <cstdint>
#include <cstring>

namespace dmw::detail {

// boost-style hash_combine over three 32- or 64-bit words (bit patterns, via memcpy).
// The result is narrowed to size_t, which is 32 bits on wasm32 and 64 bits natively,
// hence the explicit cast.
template <typename Word>
struct Hash3 {
    std::size_t operator()(const std::array<Word, 3>& k) const noexcept {
        static_assert(sizeof(Word) <= 8);
        std::uint64_t h = 0;
        for (Word w : k) {
            std::uint64_t x = 0;
            std::memcpy(&x, &w, sizeof w);
            h ^= x + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        }
        return static_cast<std::size_t>(h);
    }
};

// Undirected/directed edge key: u in the high 32 bits, v in the low 32.
inline std::uint64_t pack_pair(std::uint32_t u, std::uint32_t v) {
    return (static_cast<std::uint64_t>(u) << 32) | v;
}

}  // namespace dmw::detail
