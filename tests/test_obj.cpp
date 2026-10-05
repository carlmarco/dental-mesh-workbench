#include <array>
#include <cstdint>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "core/io.h"

using dmw::parse_obj;
using Catch::Matchers::ContainsSubstring;
using Tris = std::vector<std::array<std::uint32_t, 3>>;

// Unit tetrahedron, all faces wound counter-clockwise seen from outside.
// (Checked by hand: each face's (b-a)x(c-a) points away from the interior.)
// Reused by the half-edge tests in 2c.
static constexpr const char* kTetra =
    "# unit tetrahedron\n"
    "v 0 0 0\n"
    "v 1 0 0\n"
    "v 0 1 0\n"
    "v 0 0 1\n"
    "f 1 3 2\n"
    "f 1 2 4\n"
    "f 1 4 3\n"
    "f 2 3 4\n";

TEST_CASE("OBJ: tetrahedron - positions read, indices converted 1-based -> 0-based", "[obj]") {
    const auto r = parse_obj(kTetra);
    REQUIRE(r.ok());
    REQUIRE(r.mesh.positions.size() == 4);
    CHECK(r.mesh.positions[1].x == 1.0);
    CHECK(r.mesh.positions[3].z == 1.0);
    CHECK(r.mesh.triangles == Tris{{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}});
}

TEST_CASE("OBJ: unit cube of quads - fan triangulation keeps winding", "[obj]") {
    const char* cube =
        "v 0 0 0\nv 1 0 0\nv 1 1 0\nv 0 1 0\n"
        "v 0 0 1\nv 1 0 1\nv 1 1 1\nv 0 1 1\n"
        "f 1 4 3 2\n"   // bottom (z=0), normal -z
        "f 5 6 7 8\n"   // top
        "f 1 2 6 5\n"   // front (y=0)
        "f 3 4 8 7\n"   // back
        "f 1 5 8 4\n"   // left (x=0)
        "f 2 3 7 6\n";  // right
    const auto r = parse_obj(cube);
    REQUIRE(r.ok());
    REQUIRE(r.mesh.positions.size() == 8);
    REQUIRE(r.mesh.triangles.size() == 12);
    // Quad (a,b,c,d) -> (a,b,c), (a,c,d): same cyclic order, so same orientation.
    CHECK(r.mesh.triangles[0] == std::array<std::uint32_t, 3>{0, 3, 2});
    CHECK(r.mesh.triangles[1] == std::array<std::uint32_t, 3>{0, 2, 1});
}

TEST_CASE("OBJ: face corner formats v, v/t, v/t/n, v//n all yield the vertex index", "[obj]") {
    const char* text =
        "v 0 0 0\nv 1 0 0\nv 0 1 0\n"
        "vt 0 0\nvt 1 0\nvt 0 1\n"
        "vn 0 0 1\n"
        "f 1 2 3\n"
        "f 1/1 2/2 3/3\n"
        "f 1/1/1 2/2/1 3/3/1\n"
        "f 1//1 2//1 3//1\n";
    const auto r = parse_obj(text);
    REQUIRE(r.ok());
    CHECK(r.mesh.triangles == Tris{{0, 1, 2}, {0, 1, 2}, {0, 1, 2}, {0, 1, 2}});
}

TEST_CASE("OBJ: negative indices are relative to vertices defined so far", "[obj]") {
    const char* text =
        "v 0 0 0\nv 1 0 0\nv 0 1 0\n"
        "f -3 -2 -1\n"   // -> 0 1 2
        "v 0 0 1\n"
        "f -3 -2 -1\n";  // -> 1 2 3
    const auto r = parse_obj(text);
    REQUIRE(r.ok());
    CHECK(r.mesh.triangles == Tris{{0, 1, 2}, {1, 2, 3}});
}

TEST_CASE("OBJ: tolerates real-world noise", "[obj]") {
    // CRLF endings, blank lines, leading whitespace and tabs, trailing comments,
    // optional w coordinate, scientific notation, ignored directives.
    const char* text =
        "mtllib x.mtl\r\n"
        "o thing\r\n"
        "\r\n"
        "  v\t1e-3 -2.5 3  # comment\r\n"
        "v 1 0 0 1.0\r\n"
        "v 0 1 0\r\n"
        "g group\r\ns off\r\nusemtl m\r\n"
        "f 1 2 3 # trailing\r\n";
    const auto r = parse_obj(text);
    REQUIRE(r.ok());
    REQUIRE(r.mesh.positions.size() == 3);
    CHECK(r.mesh.positions[0].x == 1e-3);
    CHECK(r.mesh.positions[0].y == -2.5);
    CHECK(r.mesh.triangles == Tris{{0, 1, 2}});
}

TEST_CASE("OBJ: empty input is a valid, empty mesh", "[obj]") {
    const auto r = parse_obj("");
    REQUIRE(r.ok());
    CHECK(r.mesh.positions.empty());
    CHECK(r.mesh.triangles.empty());
}

TEST_CASE("OBJ: malformed input is rejected with the offending line number", "[obj]") {
    const char* tri3 = "v 0 0 0\nv 1 0 0\nv 0 1 0\n";  // lines 1-3

    SECTION("index 0 is invalid (OBJ is 1-based)") {
        const auto r = parse_obj(std::string(tri3) + "f 1 2 0\n");
        REQUIRE_FALSE(r.ok());
        CHECK_THAT(r.error, ContainsSubstring("line 4"));
    }
    SECTION("positive index out of range") {
        const auto r = parse_obj(std::string(tri3) + "f 1 2 9\n");
        REQUIRE_FALSE(r.ok());
        CHECK_THAT(r.error, ContainsSubstring("line 4"));
    }
    SECTION("negative index out of range") {
        const auto r = parse_obj(std::string(tri3) + "f -4 -2 -1\n");
        REQUIRE_FALSE(r.ok());
        CHECK_THAT(r.error, ContainsSubstring("line 4"));
    }
    SECTION("face with fewer than 3 corners") {
        const auto r = parse_obj(std::string(tri3) + "f 1 2\n");
        REQUIRE_FALSE(r.ok());
        CHECK_THAT(r.error, ContainsSubstring("line 4"));
    }
    SECTION("vertex with missing coordinate") {
        const auto r = parse_obj("v 1 2\n");
        REQUIRE_FALSE(r.ok());
        CHECK_THAT(r.error, ContainsSubstring("line 1"));
    }
    SECTION("vertex with non-numeric coordinate") {
        const auto r = parse_obj("v 0 0 0\nv 1 abc 3\n");
        REQUIRE_FALSE(r.ok());
        CHECK_THAT(r.error, ContainsSubstring("line 2"));
    }
}
