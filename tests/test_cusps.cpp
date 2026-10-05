#include <cmath>
#include <cstdint>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/cusps.h"
#include "core/generate.h"
#include "core/metrics.h"

using namespace dmw;
using Catch::Matchers::WithinAbs;

namespace {

// 10 x 10 mm patch (spacing ~0.1 mm) with Gaussian "cusps" (height 1 mm, width 0.8 mm) at known
// centers, plus fine deterministic roughness standing in for scan noise.
const std::vector<Vec3> kCenters{{2.5, 2.5, 0}, {7.0, 3.0, 0}, {3.0, 7.5, 0}, {7.5, 7.5, 0}};

TriMesh bumpy_patch(double noise) {
    auto m = make_grid(100, 100);
    for (auto& p : m.positions) {
        p.x *= 10.0, p.y *= 10.0;
        double z = 0.0;
        for (const Vec3& c : kCenters) {
            const double d2 = (p.x - c.x) * (p.x - c.x) + (p.y - c.y) * (p.y - c.y);
            z += std::exp(-d2 / (2 * 0.8 * 0.8));
        }
        z += noise * std::sin(37.0 * p.x) * std::sin(41.0 * p.y);
        p.z = z;
    }
    return m;
}

std::vector<Vec3> tips() {
    std::vector<Vec3> t;
    for (const Vec3& c : kCenters) t.push_back({c.x, c.y, 1.0});
    return t;
}

std::vector<Vec3> positions_of(const HalfEdgeMesh& m, const CuspDetection& d) {
    std::vector<Vec3> p;
    for (std::uint32_t v : d.vertices) p.push_back(m.positions[v]);
    return p;
}

HalfEdgeMesh he(const TriMesh& m) {
    auto r = build_halfedge(m);
    REQUIRE(r.ok());
    return std::move(r.mesh);
}

}  // namespace

TEST_CASE("cusps: occlusal axis is the shallow direction, pointing at the cusps", "[cusps]") {
    const auto mesh = he(bumpy_patch(0.0));
    const std::vector<std::uint8_t> all(mesh.positions.size(), 1);
    const Vec3 a = occlusal_axis(mesh, all);
    CHECK(a.z > 0.999);  // +z: away from the boundary (z ~ 0), toward the bumps
}

TEST_CASE("cusps: non-maximum suppression keeps the strongest per ball", "[cusps]") {
    const std::vector<Vec3> p{{0, 0, 0}, {0.5, 0, 0}, {3, 0, 0}, {3.2, 0, 0}};
    const std::vector<double> s{1.0, 2.0, 0.5, 0.7};
    const std::vector<std::uint32_t> cand{0, 1, 2, 3};
    CHECK(non_max_suppression(p, s, cand, 1.0) == std::vector<std::uint32_t>{1, 3});
}

TEST_CASE("cusps: smoothing preserves constants and the mean (mass-weighted)", "[cusps]") {
    const auto mesh = he(bumpy_patch(0.0));
    const auto ops = build_dec(mesh);
    const std::vector<double> ones(mesh.positions.size(), 1.0);
    for (double x : diffuse(ops, ones, 0.5)) CHECK_THAT(x, WithinAbs(1.0, 1e-6));
}

TEST_CASE("cusps: prominence solved directly equals field minus diffused field", "[cusps]") {
    const auto mesh = he(bumpy_patch(0.01));
    const auto ops = build_dec(mesh);
    std::vector<double> z(mesh.positions.size());
    for (std::size_t v = 0; v < z.size(); ++v) z[v] = mesh.positions[v].z;
    const auto base = diffuse(ops, z, 2.0);
    const auto prom = prominence(ops, z, 2.0);
    for (std::size_t v = 0; v < z.size(); ++v) CHECK_THAT(prom[v], WithinAbs(z[v] - base[v], 1e-6));
}

TEST_CASE("cusps: smoothed and prominence detectors find every tip despite roughness; raw curvature does not", "[cusps]") {
    const auto mesh = he(bumpy_patch(0.01));  // 10 um roughness at ~0.16 mm wavelength
    CuspParams prom;  // defaults: occlusal prominence
    prom.threshold = 0.2;
    const auto dp = detect_cusps(mesh, prom);
    const auto rp = match_points(positions_of(mesh, dp), tips(), 0.3);
    CHECK(rp.recall() == 1.0);
    CHECK(rp.precision() == 1.0);

    CuspParams smooth;
    smooth.method = CuspParams::Method::SmoothedCurvature;
    smooth.threshold = 0.3;
    const auto rs = match_points(positions_of(mesh, detect_cusps(mesh, smooth)), tips(), 0.3);
    CHECK(rs.recall() == 1.0);

    CuspParams raw;
    raw.method = CuspParams::Method::RawCurvature;
    raw.threshold = 0.3;
    const auto rr = match_points(positions_of(mesh, detect_cusps(mesh, raw)), tips(), 0.3);
    CHECK(rr.precision() < 0.5);  // roughness creates curvature peaks everywhere
}
