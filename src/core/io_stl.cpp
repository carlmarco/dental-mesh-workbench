#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>
#include <string_view>

#include "core/io.h"
#include "detail/parse_number.h"

namespace dmw {
namespace {

using detail::parse_double;

// Binary STL layout (all little-endian):
//   [0, 80)   header (arbitrary bytes; may well start with "solid")
//   [80, 84)  uint32 triangle count N
//   then N records of 50 bytes: float32 normal[3], float32 v0[3], v1[3], v2[3], uint16 attr
constexpr std::size_t kHeaderBytes = 80;
constexpr std::size_t kPreambleBytes = 84;
constexpr std::size_t kRecordBytes = 50;
constexpr std::size_t kNormalBytes = 12;
constexpr std::size_t kVertexBytes = 12;

static_assert(sizeof(float) == 4 && std::numeric_limits<float>::is_iec559,
              "binary STL decoding assumes IEEE-754 float32");

// Assemble from individual bytes with shifts: correct on any host byte order, and no
// unaligned load (record offsets 84 + 50*i are not 4-byte aligned for odd i).
std::uint32_t read_u32_le(const std::uint8_t* p) {
    return static_cast<std::uint32_t>(p[0]) | (static_cast<std::uint32_t>(p[1]) << 8) |
           (static_cast<std::uint32_t>(p[2]) << 16) | (static_cast<std::uint32_t>(p[3]) << 24);
}

// memcpy is the defined way to reinterpret bits; casting a uint8_t* to float* would
// violate strict aliasing (and alignment). Compilers turn this into a single move.
float read_f32_le(const std::uint8_t* p) {
    const std::uint32_t bits = read_u32_le(p);
    float f = 0.0f;
    std::memcpy(&f, &bits, sizeof f);
    return f;
}

// D19: binary iff the size matches the layout exactly. 64-bit arithmetic: with
// N = 2^32-1, 84 + 50*N ~ 2.1e11 overflows 32 bits but not 64.
bool is_binary_stl(std::span<const std::uint8_t> bytes) {
    if (bytes.size() < kPreambleBytes) return false;
    const std::uint64_t n = read_u32_le(bytes.data() + kHeaderBytes);
    return static_cast<std::uint64_t>(bytes.size()) == kPreambleBytes + kRecordBytes * n;
}

LoadResult parse_binary(std::span<const std::uint8_t> bytes) {
    const std::size_t n = read_u32_le(bytes.data() + kHeaderBytes);
    // 3 soup vertices per triangle must fit in uint32 indices (D13).
    if (n > std::numeric_limits<std::uint32_t>::max() / 3) {
        return {TriMesh{}, "binary STL: too many triangles for 32-bit indices"};
    }
    LoadResult result;
    result.mesh.positions.reserve(3 * n);
    result.mesh.triangles.reserve(n);
    for (std::size_t t = 0; t < n; ++t) {
        const std::uint8_t* rec = bytes.data() + kPreambleBytes + kRecordBytes * t;
        for (std::size_t c = 0; c < 3; ++c) {
            const std::uint8_t* v = rec + kNormalBytes + kVertexBytes * c;  // normal skipped (D19)
            const float x = read_f32_le(v), y = read_f32_le(v + 4), z = read_f32_le(v + 8);
            if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
                return {TriMesh{}, "triangle " + std::to_string(t) + ": non-finite coordinate"};
            }
            result.mesh.positions.push_back({x, y, z});  // float -> double is exact
        }
        const auto base = static_cast<std::uint32_t>(3 * t);
        result.mesh.triangles.push_back({base, base + 1, base + 2});
    }
    return result;
}

// Whitespace-delimited tokens that remember their line number for error messages.
class Tokenizer {
public:
    explicit Tokenizer(std::string_view text) : text_(text) {}

    std::string_view next() {
        while (pos_ < text_.size() && is_space(text_[pos_])) {
            if (text_[pos_] == '\n') ++line_;
            ++pos_;
        }
        token_line_ = line_;
        const std::size_t begin = pos_;
        while (pos_ < text_.size() && !is_space(text_[pos_])) ++pos_;
        return text_.substr(begin, pos_ - begin);
    }

    // Discards the rest of the current line (the free-form name after "solid").
    void skip_line() {
        while (pos_ < text_.size() && text_[pos_] != '\n') ++pos_;
    }

    std::size_t line() const { return token_line_; }

private:
    static bool is_space(char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\n'; }

    std::string_view text_;
    std::size_t pos_ = 0;
    std::size_t line_ = 1;
    std::size_t token_line_ = 1;
};

// Grammar:
//   solid [name]
//   { facet normal nx ny nz  outer loop  (vertex x y z){3}  endloop  endfacet }*
//   endsolid [name]
LoadResult parse_ascii(std::string_view text) {
    Tokenizer tok(text);
    LoadResult result;
    auto fail = [&](std::string_view why) {
        return LoadResult{TriMesh{}, "line " + std::to_string(tok.line()) + ": " + std::string(why)};
    };
    auto expect = [&](std::string_view word) { return tok.next() == word; };
    auto read3 = [&](Vec3& p) {
        return parse_double(tok.next(), p.x) && parse_double(tok.next(), p.y) &&
               parse_double(tok.next(), p.z);
    };

    if (!expect("solid")) {
        return {TriMesh{}, "not STL: size doesn't match the binary layout and text doesn't start with 'solid'"};
    }
    tok.skip_line();

    for (;;) {
        const std::string_view word = tok.next();
        if (word == "endsolid") break;
        if (word.empty()) return fail("unexpected end of file, expected 'endsolid'");
        if (word != "facet") return fail("expected 'facet' or 'endsolid'");

        Vec3 normal;  // parsed for validation, then ignored (D19)
        if (!expect("normal") || !read3(normal)) return fail("expected 'normal nx ny nz'");
        if (!expect("outer") || !expect("loop")) return fail("expected 'outer loop'");
        for (int c = 0; c < 3; ++c) {
            Vec3 p;
            if (!expect("vertex")) return fail("expected 'vertex'");
            if (!read3(p)) return fail("vertex needs 3 numeric coordinates");
            result.mesh.positions.push_back(p);
        }
        if (!expect("endloop")) return fail("expected 'endloop'");
        if (!expect("endfacet")) return fail("expected 'endfacet'");

        if (result.mesh.positions.size() > std::numeric_limits<std::uint32_t>::max()) {
            return fail("too many vertices for 32-bit indices");
        }
        const auto base = static_cast<std::uint32_t>(result.mesh.positions.size() - 3);
        result.mesh.triangles.push_back({base, base + 1, base + 2});
    }
    return result;
}

}  // namespace

LoadResult parse_stl(std::span<const std::uint8_t> bytes) {
    if (is_binary_stl(bytes)) return parse_binary(bytes);
    // Viewing bytes as chars is allowed: char may alias any object type.
    return parse_ascii({reinterpret_cast<const char*>(bytes.data()), bytes.size()});
}

}  // namespace dmw
