#include <cstdint>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "core/generate.h"
#include "core/io.h"
#include "core/topology.h"
#include "core/weld.h"

using namespace dmw;

TEST_CASE("io: OBJ round-trip is exact (17 significant digits)", "[io]") {
    const auto m = make_torus(9, 5);  // irrational-looking coordinates from sin/cos
    const auto r = parse_obj(write_obj(m));
    REQUIRE(r.ok());
    REQUIRE(r.mesh.positions.size() == m.positions.size());
    for (std::size_t v = 0; v < m.positions.size(); ++v) {
        CHECK(r.mesh.positions[v].x == m.positions[v].x);  // bit-exact
        CHECK(r.mesh.positions[v].y == m.positions[v].y);
        CHECK(r.mesh.positions[v].z == m.positions[v].z);
    }
    CHECK(r.mesh.triangles == m.triangles);
}

TEST_CASE("io: binary STL round-trip welds back to the same topology", "[io]") {
    const auto m = make_plate_with_handle();
    const auto bytes = write_stl_binary(m);
    CHECK(bytes.size() == 84 + 50 * m.triangles.size());
    const auto r = parse_stl(bytes);
    REQUIRE(r.ok());
    CHECK(r.mesh.positions.size() == 3 * m.triangles.size());  // a soup, as STL stores it
    const auto welded = weld_vertices(r.mesh);  // exact weld: float32 corners are bit-identical
    CHECK(welded.positions.size() == m.positions.size());
    const auto topo = analyze_topology(welded);
    REQUIRE(topo.components.size() == 1);
    CHECK(topo.components[0].genus == std::optional<std::uint32_t>{1});
    CHECK(static_cast<float>(welded.positions[0].x) == static_cast<float>(m.positions[0].x));
}
