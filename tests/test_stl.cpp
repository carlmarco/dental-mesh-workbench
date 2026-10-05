#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string_view>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "core/io.h"

using dmw::parse_stl;
using Catch::Matchers::ContainsSubstring;
using Tris = std::vector<std::array<std::uint32_t, 3>>;

namespace {

std::vector<std::uint8_t> bytes_of(std::string_view s) { return {s.begin(), s.end()}; }

// Explicit little-endian writes, independent of the host's byte order.
void put_u32(std::vector<std::uint8_t>& out, std::uint32_t v) {
    for (int i = 0; i < 4; ++i) out.push_back(static_cast<std::uint8_t>(v >> (8 * i)));
}
void put_f32(std::vector<std::uint8_t>& out, float f) {
    std::uint32_t bits = 0;
    std::memcpy(&bits, &f, sizeof bits);
    put_u32(out, bits);
}

// Binary STL layout: 80-byte header | uint32 triangle count | per triangle:
// 12 x float32 (normal, v0, v1, v2) + uint16 attribute byte count. Little-endian.
std::vector<std::uint8_t> binary_stl(const std::vector<std::array<float, 9>>& tris,
                                     std::string_view header = "") {
    std::vector<std::uint8_t> out(80, 0);
    std::copy_n(header.begin(), std::min<std::size_t>(header.size(), 80), out.begin());
    put_u32(out, static_cast<std::uint32_t>(tris.size()));
    for (const auto& t : tris) {
        for (int i = 0; i < 3; ++i) put_f32(out, 0.0f);  // normal: ignored by the loader
        for (float c : t) put_f32(out, c);
        out.push_back(0);
        out.push_back(0);
    }
    return out;
}

}  // namespace

TEST_CASE("STL ASCII: single facet -> 3 soup vertices, normal ignored", "[stl]") {
    const auto r = parse_stl(bytes_of(
        "solid t\n"
        "  facet normal 0 0 1\n"
        "    outer loop\n"
        "      vertex 0 0 0\n"
        "      vertex 1 0 0\n"
        "      vertex 0 1 0\n"
        "    endloop\n"
        "  endfacet\n"
        "endsolid t\n"));
    REQUIRE(r.ok());
    REQUIRE(r.mesh.positions.size() == 3);
    CHECK(r.mesh.positions[1].x == 1.0);
    CHECK(r.mesh.triangles == Tris{{0, 1, 2}});
}

TEST_CASE("STL ASCII: facet with too few vertices reports its line", "[stl]") {
    const auto r = parse_stl(bytes_of(
        "solid t\n"
        "facet normal 0 0 1\n"
        "outer loop\n"
        "vertex 0 0 0\n"
        "vertex 1 0 0\n"
        "endloop\n"  // line 6: a third vertex was expected here
        "endfacet\n"
        "endsolid t\n"));
    REQUIRE_FALSE(r.ok());
    CHECK_THAT(r.error, ContainsSubstring("line 6"));
}

TEST_CASE("STL binary: single triangle; float32 converts exactly to double", "[stl]") {
    const auto r = parse_stl(binary_stl({{0.1f, 0, 0, 1, 0, 0, 0, 1, 0}}));
    REQUIRE(r.ok());
    REQUIRE(r.mesh.positions.size() == 3);
    // The double holds the float32 value, which is NOT 0.1 (0.1 isn't representable).
    CHECK(r.mesh.positions[0].x == static_cast<double>(0.1f));
    CHECK(r.mesh.positions[1].x == 1.0);
    CHECK(r.mesh.triangles == Tris{{0, 1, 2}});
}

TEST_CASE("STL binary: header beginning with 'solid' is still detected as binary", "[stl]") {
    // Many CAD exporters do this; prefix-sniffing would misparse it as ASCII.
    const auto r = parse_stl(binary_stl({{0, 0, 0, 1, 0, 0, 0, 1, 0}}, "solid exported by CAD"));
    REQUIRE(r.ok());
    CHECK(r.mesh.positions.size() == 3);
}

TEST_CASE("STL binary: zero triangles is a valid empty mesh", "[stl]") {
    const auto r = parse_stl(binary_stl({}));
    REQUIRE(r.ok());
    CHECK(r.mesh.triangles.empty());
}

TEST_CASE("STL: rejected inputs", "[stl]") {
    SECTION("truncated binary (size doesn't match the triangle count)") {
        auto b = binary_stl({{0, 0, 0, 1, 0, 0, 0, 1, 0}});
        b.pop_back();
        CHECK_FALSE(parse_stl(b).ok());
    }
    SECTION("huge triangle count in a tiny file: size check must not overflow") {
        auto b = binary_stl({});
        std::fill(b.begin() + 80, b.begin() + 84, std::uint8_t{0xFF});  // N = 2^32 - 1
        CHECK_FALSE(parse_stl(b).ok());
    }
    SECTION("NaN coordinate") {
        const float nan = std::numeric_limits<float>::quiet_NaN();
        CHECK_FALSE(parse_stl(binary_stl({{0, 0, 0, nan, 0, 0, 0, 1, 0}})).ok());
    }
    SECTION("too short to be either format") {
        CHECK_FALSE(parse_stl(bytes_of("abc")).ok());
    }
    SECTION("empty input (unlike OBJ, STL requires a header or 'solid')") {
        CHECK_FALSE(parse_stl({}).ok());
    }
}
