#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/generate.h"
#include "core/geodesic.h"

using namespace dmw;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using Ids = std::vector<std::uint32_t>;

namespace {

HalfEdgeMesh he(const TriMesh& m) {
    auto r = build_halfedge(m);
    REQUIRE(r.ok());
    return std::move(r.mesh);
}

double great_circle(const Vec3& a, const Vec3& b, double r) {
    const double c = (a.x * b.x + a.y * b.y + a.z * b.z) / (r * r);
    return r * std::acos(std::clamp(c, -1.0, 1.0));
}

}  // namespace

// ---- Exact properties -------------------------------------------------------------------

TEST_CASE("dijkstra: exact along grid lines, never below Euclidean on a flat convex domain", "[geodesic]") {
    const std::uint32_t n = 8;
    const auto mesh = he(make_grid(n, n));
    const auto d = dijkstra_distance(mesh, Ids{0});  // corner (0, 0)
    CHECK(d[0] == 0.0);
    CHECK_THAT(d[n], WithinAbs(1.0, 1e-12));             // (1, 0): straight along the bottom row
    CHECK_THAT(d[n * (n + 1)], WithinAbs(1.0, 1e-12));   // (0, 1)
    CHECK_THAT(d[n * (n + 1) + n], WithinAbs(std::sqrt(2.0), 1e-12));  // (1, 1): along the diagonals
    for (std::size_t v = 0; v < d.size(); ++v) {
        const Vec3& p = mesh.positions[v];
        CHECK(d[v] >= std::hypot(p.x, p.y) - 1e-12);
    }
}

TEST_CASE("heat: zero at the source, finite and non-negative elsewhere", "[geodesic]") {
    const auto mesh = he(make_icosphere(3));
    const auto r = HeatGeodesics(mesh).distance(Ids{0});
    REQUIRE(r.ok());
    CHECK(r.distance[0] == 0.0);
    for (double d : r.distance) {
        CHECK(std::isfinite(d));
        CHECK(d > -1e-3);  // tiny negative values near the source are discretization noise
    }
}

TEST_CASE("heat and dijkstra: unreachable components are NaN", "[geodesic]") {
    auto m = make_icosphere(2);
    const auto first_far = static_cast<std::uint32_t>(m.positions.size());
    auto far = make_icosphere(1);
    for (auto& p : far.positions) p.x += 5.0;
    append(m, far);
    const auto mesh = he(m);
    const auto heat = HeatGeodesics(mesh).distance(Ids{0});
    REQUIRE(heat.ok());
    const auto dij = dijkstra_distance(mesh, Ids{0});
    CHECK(std::isfinite(heat.distance[1]));
    CHECK(std::isnan(heat.distance[first_far]));
    CHECK(std::isnan(dij[first_far]));
}

TEST_CASE("heat: multiple sources give the distance to the nearest one", "[geodesic]") {
    const auto mesh = he(make_icosphere(3));
    const Ids both{0, 3};  // icosahedron vertices 0 = (-1, phi, 0) and 3 = (1, -phi, 0): antipodal
    const auto r = HeatGeodesics(mesh).distance(both);
    REQUIRE(r.ok());
    CHECK_THAT(r.distance[0], WithinAbs(0.0, 1e-9));
    CHECK_THAT(r.distance[3], WithinAbs(0.0, 0.05));  // sources are shifted to mean zero
    double far = 0.0;
    for (double d : r.distance) far = std::max(far, d);
    CHECK(far < 0.6 * std::numbers::pi);  // single-source max would be ~pi
}

TEST_CASE("heat: invalid sources are rejected", "[geodesic]") {
    const auto mesh = he(make_icosphere(1));
    const HeatGeodesics solver(mesh);
    CHECK_FALSE(solver.distance(Ids{}).ok());
    CHECK_FALSE(solver.distance(Ids{100000}).ok());
}

// ---- Accuracy (thresholds from measured values, DECISIONS.md D50, with margin) ------------

namespace {
struct Errors {
    double mean = 0.0, max = 0.0;
};
template <typename Ref>
Errors errors(const std::vector<double>& d, Ref ref) {
    Errors e;
    for (std::size_t v = 0; v < d.size(); ++v) {
        const double x = std::abs(d[v] - ref(v));
        e.mean += x;
        e.max = std::max(e.max, x);
    }
    e.mean /= static_cast<double>(d.size());
    return e;
}
}  // namespace

TEST_CASE("heat: sphere distance converges to great-circle distance; Dijkstra does not", "[geodesic][convergence]") {
    // Source vertex 5 on purpose: not the vertex pinned in the Poisson solve.
    const Ids src{5};
    std::vector<double> heat_mean, dij_mean;
    for (std::uint32_t s = 2; s <= 4; ++s) {
        const auto mesh = he(make_icosphere(s, 1.0));
        auto ref = [&](std::size_t v) { return great_circle(mesh.positions[5], mesh.positions[v], 1.0); };
        const auto r = HeatGeodesics(mesh).distance(src);
        REQUIRE(r.ok());
        heat_mean.push_back(errors(r.distance, ref).mean);
        dij_mean.push_back(errors(dijkstra_distance(mesh, src), ref).mean);
        if (s == 4) {
            CHECK(errors(r.distance, ref).mean < 2e-2);  // measured 1.32e-2
            CHECK(errors(r.distance, ref).max < 4e-2);   // measured 2.92e-2
        }
    }
    CHECK(heat_mean[1] < heat_mean[0]);
    CHECK(heat_mean[2] < heat_mean[1]);
    for (std::size_t i = 0; i < heat_mean.size(); ++i) CHECK(heat_mean[i] < 0.5 * dij_mean[i]);
    CHECK(dij_mean[2] > 0.1);  // edge paths zigzag: ~0.12 at every resolution (measured)
}

TEST_CASE("heat: flat grid distance converges to Euclidean distance", "[geodesic][convergence]") {
    std::vector<double> mean;
    for (std::uint32_t n : {8u, 16u, 32u}) {
        const auto mesh = he(make_grid(n, n));
        const Ids center{(n / 2) * (n + 1) + n / 2};
        const auto r = HeatGeodesics(mesh).distance(center);
        REQUIRE(r.ok());
        mean.push_back(errors(r.distance, [&](std::size_t v) {
                           return std::hypot(mesh.positions[v].x - 0.5, mesh.positions[v].y - 0.5);
                       }).mean);
    }
    CHECK(mean[1] < mean[0]);
    CHECK(mean[2] < mean[1]);
    CHECK(mean[2] < 1e-2);  // measured 7.78e-3
}

TEST_CASE("heat: irregular (jittered) sphere stays accurate with the direct solver", "[geodesic]") {
    // Regression for D47: with CG the error grew under refinement; with LDL^T it does not.
    auto m = make_icosphere(4, 1.0);
    const double h = 1.1 / 16.0;
    for (std::size_t v = 1; v < m.positions.size(); ++v) {
        auto& p = m.positions[v];
        const double k = static_cast<double>(v);
        auto j = [](double x) { const double s = std::sin(x * 12.9898) * 43758.5453; return 2 * (s - std::floor(s)) - 1; };
        p.x += 0.25 * h * j(3 * k), p.y += 0.25 * h * j(3 * k + 1), p.z += 0.25 * h * j(3 * k + 2);
        const double n = std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        p = {p.x / n, p.y / n, p.z / n};
    }
    const auto mesh = he(m);
    const auto r = HeatGeodesics(mesh).distance(Ids{5});
    REQUIRE(r.ok());
    CHECK(errors(r.distance, [&](std::size_t v) { return great_circle(mesh.positions[5], mesh.positions[v], 1.0); }).mean < 2e-2);
}
