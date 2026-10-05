#include <charconv>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "core/io.h"
#include "detail/parse_number.h"

namespace dmw {
namespace {

using detail::parse_double;

// '\r' counts as whitespace so CRLF files never leave a stray '\r' on the last token.
constexpr bool is_space(char c) { return c == ' ' || c == '\t' || c == '\r'; }

// Pops the next whitespace-delimited token off the front of `s`; empty when none remain.
std::string_view next_token(std::string_view& s) {
    std::size_t begin = 0;
    while (begin < s.size() && is_space(s[begin])) ++begin;
    std::size_t end = begin;
    while (end < s.size() && !is_space(s[end])) ++end;
    const std::string_view tok = s.substr(begin, end - begin);
    s.remove_prefix(end);
    return tok;
}

// Resolves one face corner ("i", "i/t", "i/t/n", "i//n") to a 0-based vertex index.
// OBJ indices are 1-based; negative indices count back from the last vertex defined
// so far (-1 = most recent). Returns an empty string on success, else the reason.
std::string resolve_corner(std::string_view corner, std::size_t num_vertices, std::uint32_t& out) {
    const std::string_view idx_text = corner.substr(0, corner.find('/'));
    std::int64_t idx = 0;
    const char* first = idx_text.data();
    const char* last = first + idx_text.size();
    const auto [ptr, ec] = std::from_chars(first, last, idx);  // integer from_chars is fine everywhere
    if (ec != std::errc{} || ptr != last) {
        return "bad vertex index '" + std::string(corner) + "'";
    }
    const auto n = static_cast<std::int64_t>(num_vertices);
    const std::int64_t resolved = idx > 0 ? idx - 1 : n + idx;
    if (idx == 0 || resolved < 0 || resolved >= n) {
        return "vertex index " + std::to_string(idx) + " out of range (" + std::to_string(n) +
               " vertices defined so far)";
    }
    out = static_cast<std::uint32_t>(resolved);
    return {};
}

LoadResult error_at(std::size_t line_no, std::string_view why) {
    // Mesh left empty: never hand back a half-parsed mesh alongside an error.
    return {TriMesh{}, "line " + std::to_string(line_no) + ": " + std::string(why)};
}

}  // namespace

LoadResult parse_obj(std::string_view text) {
    LoadResult result;
    std::vector<std::uint32_t> poly;  // reused across faces to avoid reallocating
    std::size_t line_no = 0;
    std::size_t pos = 0;

    while (pos < text.size()) {
        std::size_t nl = text.find('\n', pos);
        if (nl == std::string_view::npos) nl = text.size();
        std::string_view line = text.substr(pos, nl - pos);
        pos = nl + 1;
        ++line_no;

        if (const auto hash = line.find('#'); hash != std::string_view::npos) {
            line = line.substr(0, hash);
        }

        const std::string_view keyword = next_token(line);
        if (keyword == "v") {
            if (result.mesh.positions.size() >= std::numeric_limits<std::uint32_t>::max()) {
                return error_at(line_no, "too many vertices for 32-bit indices");
            }
            Vec3 p;
            if (!parse_double(next_token(line), p.x) || !parse_double(next_token(line), p.y) ||
                !parse_double(next_token(line), p.z)) {
                return error_at(line_no, "vertex needs 3 numeric coordinates");
            }
            // Further values (w, or the common "r g b" vertex-color extension) are ignored.
            result.mesh.positions.push_back(p);
        } else if (keyword == "f") {
            poly.clear();
            for (auto corner = next_token(line); !corner.empty(); corner = next_token(line)) {
                std::uint32_t v = 0;
                if (const auto why = resolve_corner(corner, result.mesh.positions.size(), v);
                    !why.empty()) {
                    return error_at(line_no, why);
                }
                poly.push_back(v);
            }
            if (poly.size() < 3) return error_at(line_no, "face needs at least 3 corners");
            // Fan triangulation (D16): (v0, vi, vi+1) keeps the polygon's cyclic order,
            // hence its orientation. Degenerate faces (repeated indices) are syntactically
            // valid and pass through; flagging them is the diagnostics' job (M3).
            for (std::size_t i = 1; i + 1 < poly.size(); ++i) {
                result.mesh.triangles.push_back({poly[0], poly[i], poly[i + 1]});
            }
        }
        // Everything else (vt, vn, o, g, s, usemtl, mtllib, blank lines) is ignored.
    }
    return result;
}

}  // namespace dmw
