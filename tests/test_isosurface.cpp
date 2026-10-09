#include <cmath>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/generate.h"
#include "core/halfedge.h"
#include "core/isosurface.h"
#include "core/sdf.h"
#include "core/topology.h"

using namespace dmw;
using Catch::Matchers::WithinAbs;

namespace {

Grid3 analytic(double h, double half, double (*f)(const Vec3&)) {
    Grid3 g;
    g.h = h;
    const auto n = static_cast<std::size_t>(std::ceil(2.0 * half / h)) + 1;
    g.n = {n, n, n};
    g.origin = {-half, -half, -half};
    g.values.resize(n * n * n);
    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t j = 0; j < n; ++j)
            for (std::size_t i = 0; i < n; ++i) g.values[g.index(i, j, k)] = f(g.position(i, j, k));
    return g;
}
double sphere_f(const Vec3& p) { return std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z) - 1.0; }
double torus_f(const Vec3& p) {
    const double q = std::sqrt(p.x * p.x + p.y * p.y) - 1.0;
    return std::sqrt(q * q + p.z * p.z) - 0.4;
}
double signed_volume(const TriMesh& m) {
    double v = 0.0;
    for (const auto& t : m.triangles) {
        const Vec3 &a = m.positions[t[0]], &b = m.positions[t[1]], &c = m.positions[t[2]];
        v += (a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x) + a.z * (b.x * c.y - b.y * c.x)) / 6.0;
    }
    return v;
}

}  // namespace

TEST_CASE("isosurface: sphere field - closed genus-0 manifold, vertices on the radius, outward normals", "[isosurface]") {
    const Grid3 g = analytic(0.05, 1.5, sphere_f);
    for (double iso : {0.0, 0.2}) {
        const TriMesh m = extract_isosurface(g, iso);
        const auto topo = analyze_topology(m);
        REQUIRE(topo.components.size() == 1);
        CHECK(topo.components[0].manifold);
        CHECK(topo.components[0].boundary_loops == 0);
        CHECK(topo.components[0].genus.value_or(99u) == 0);
        CHECK(build_halfedge(m).ok());  // consistently oriented
        double worst = 0.0;
        for (const Vec3& p : m.positions) worst = std::max(worst, std::abs(std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z) - (1.0 + iso)));
        CHECK(worst < 0.01);  // linear interpolation of a smooth field: O(h^2)
        CHECK_THAT(signed_volume(m), WithinAbs(4.0 / 3.0 * std::acos(-1.0) * std::pow(1.0 + iso, 3), 0.02));  // > 0: outward
    }
}

TEST_CASE("isosurface: torus field - genus 1", "[isosurface]") {
    const TriMesh m = extract_isosurface(analytic(0.05, 1.6, torus_f), 0.0);
    const auto topo = analyze_topology(m);
    REQUIRE(topo.components.size() == 1);
    CHECK(topo.components[0].genus.value_or(99u) == 1);
    CHECK(build_halfedge(m).ok());
}

TEST_CASE("isosurface: offsets of a mesh sphere through the signed heat distance", "[isosurface]") {
    SignedHeatParams prm;
    prm.h = 0.04;
    const Grid3 phi = signed_heat_distance(make_icosphere(5), prm);
    for (double d : {0.2, -0.3}) {
        const TriMesh m = extract_isosurface(phi, d);
        double sum = 0.0;
        for (const Vec3& p : m.positions) sum += std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        CHECK_THAT(sum / double(m.positions.size()), WithinAbs(1.0 + d, 0.01));
        const auto topo = analyze_topology(m);
        REQUIRE(topo.components.size() == 1);
        CHECK(topo.components[0].genus.value_or(99u) == 0);
    }
}
