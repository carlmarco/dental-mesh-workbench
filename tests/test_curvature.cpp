#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <numeric>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/curvature.h"
#include "core/generate.h"
#include "core/halfedge.h"

using namespace dmw;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using std::numbers::pi;

namespace {

HalfEdgeMesh he(const TriMesh& m) {
    auto r = build_halfedge(m);
    REQUIRE(r.ok());
    return std::move(r.mesh);
}

double sum(const std::vector<double>& v) { return std::accumulate(v.begin(), v.end(), 0.0); }

double total_area(const TriMesh& m) {
    double a = 0.0;
    for (const auto& t : m.triangles) {
        const Vec3 &p = m.positions[t[0]], &q = m.positions[t[1]], &r = m.positions[t[2]];
        const double ux = q.x - p.x, uy = q.y - p.y, uz = q.z - p.z;
        const double vx = r.x - p.x, vy = r.y - p.y, vz = r.z - p.z;
        const double cx = uy * vz - uz * vy, cy = uz * vx - ux * vz, cz = ux * vy - uy * vx;
        a += 0.5 * std::sqrt(cx * cx + cy * cy + cz * cz);
    }
    return a;
}

// Largest |H - expected| over interior vertices.
double max_mean_error(const CurvatureField& c, double expected) {
    double e = 0.0;
    for (std::size_t v = 0; v < c.mean.size(); ++v) {
        if (!c.boundary[v]) e = std::max(e, std::abs(c.mean[v] - expected));
    }
    return e;
}

}  // namespace

// ---- Exact identities -------------------------------------------------------------------

TEST_CASE("curvature: Gauss-Bonnet on closed meshes, sum of defects = 2 pi chi", "[curvature]") {
    // Exact for ANY triangulation: sum(2pi - sum theta) = 2pi V - pi F, and 3F = 2E when
    // closed, so it equals 2pi (V - E + F). Only floating-point roundoff remains.
    for (std::uint32_t s = 0; s <= 3; ++s) {
        CHECK_THAT(sum(compute_curvature(he(make_icosphere(s))).angle_defect), WithinAbs(4 * pi, 1e-10));
    }
    CHECK_THAT(sum(compute_curvature(he(make_torus(16, 8))).angle_defect), WithinAbs(0.0, 1e-10));
}

TEST_CASE("curvature: Gauss-Bonnet with boundary uses pi - sum theta at boundary vertices", "[curvature]") {
    // Disk: chi = 1. The four grid corners contribute pi - pi/2 each; edges contribute 0.
    CHECK_THAT(sum(compute_curvature(he(make_grid(5, 4))).angle_defect), WithinAbs(2 * pi * 1, 1e-10));
    // Two square holes: each hole corner has sum theta = 3pi/2, so -pi/2 x 4 = -2pi per hole.
    CHECK_THAT(sum(compute_curvature(he(make_grid_with_holes(6, 6, {{1, 1}, {4, 3}}))).angle_defect),
               WithinAbs(2 * pi * -1, 1e-10));
    CHECK_THAT(sum(compute_curvature(he(make_cylinder(10, 3))).angle_defect), WithinAbs(0.0, 1e-10));
}

TEST_CASE("curvature: mixed areas tile the surface exactly", "[curvature]") {
    for (const TriMesh& m : {make_icosphere(2), make_torus(12, 6), make_grid(4, 4), make_cylinder(9, 2)}) {
        CHECK_THAT(sum(compute_curvature(he(m)).mixed_area), WithinRel(total_area(m), 1e-12));
    }
}

TEST_CASE("curvature: mixed area on a non-obtuse triangle is the Voronoi split", "[curvature]") {
    // Equilateral: by symmetry each vertex's Voronoi region is one third of the area.
    const double h = std::sqrt(3.0) / 2;
    const auto c = compute_curvature(he({{{0, 0, 0}, {1, 0, 0}, {0.5, h, 0}}, {{0, 1, 2}}}));
    const double third = (0.5 * h) / 3;
    for (double a : c.mixed_area) CHECK_THAT(a, WithinRel(third, 1e-12));
}

TEST_CASE("curvature: mixed area on an obtuse triangle is area/2 at the obtuse vertex, area/4 else", "[curvature]") {
    // Vertex 2 at (2, 0.5): angle there is obtuse. Area = 1/2 * 4 * 0.5 = 1. The Voronoi
    // region would extend outside the triangle (circumcenter outside), hence the fallback.
    const auto c = compute_curvature(he({{{0, 0, 0}, {4, 0, 0}, {2, 0.5, 0}}, {{0, 1, 2}}}));
    CHECK_THAT(c.mixed_area[0], WithinRel(0.25, 1e-12));
    CHECK_THAT(c.mixed_area[1], WithinRel(0.25, 1e-12));
    CHECK_THAT(c.mixed_area[2], WithinRel(0.5, 1e-12));
}

TEST_CASE("curvature: flat interior vertices have zero curvature, even on irregular meshes", "[curvature]") {
    // Jitter interior vertices within the plane. The cotan formula reproduces linear
    // functions exactly, so the mean-curvature normal of any planar interior vertex is 0.
    auto m = make_grid(6, 6);
    for (std::size_t v = 0; v < m.positions.size(); ++v) {
        auto& p = m.positions[v];
        if (p.x > 0 && p.x < 1 && p.y > 0 && p.y < 1) {
            p.x += 0.03 * std::sin(17.0 * static_cast<double>(v));
            p.y += 0.03 * std::cos(11.0 * static_cast<double>(v));
        }
    }
    const auto c = compute_curvature(he(m));
    for (std::size_t v = 0; v < m.positions.size(); ++v) {
        if (c.boundary[v]) continue;
        CHECK_THAT(c.gaussian[v], WithinAbs(0.0, 1e-9));
        CHECK_THAT(c.mean[v], WithinAbs(0.0, 1e-9));
        CHECK_THAT(c.mean_normal[v].x, WithinAbs(0.0, 1e-9));
        CHECK_THAT(c.mean_normal[v].y, WithinAbs(0.0, 1e-9));
    }
}

// ---- Boundary handling and degeneracy --------------------------------------------------

TEST_CASE("curvature: boundary vertices are flagged and have no pointwise curvature", "[curvature]") {
    const auto c = compute_curvature(he(make_grid(3, 3)));
    // 4x4 vertices; the 2x2 interior block is 5, 6, 9, 10.
    for (std::uint32_t v = 0; v < 16; ++v) {
        const bool interior = (v == 5 || v == 6 || v == 9 || v == 10);
        CHECK(c.boundary[v] == (interior ? 0 : 1));
        CHECK(std::isnan(c.gaussian[v]) == !interior);
        CHECK(std::isnan(c.mean[v]) == !interior);
    }
}

TEST_CASE("curvature: zero-area faces are listed and don't poison the sums", "[curvature]") {
    // Face 2 is collinear (0,0)-(2,0)-(1,0). Vertex 1 is interior: 90 + 90 + 180 degrees.
    const TriMesh m{{{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {1, 1, 0}}, {{0, 1, 3}, {1, 2, 3}, {0, 2, 1}}};
    const auto c = compute_curvature(he(m));
    CHECK(c.degenerate_faces == std::vector<std::uint32_t>{2});
    REQUIRE_FALSE(c.boundary[1]);
    CHECK_THAT(c.angle_defect[1], WithinAbs(0.0, 1e-12));
    CHECK(std::isfinite(c.mean[1]));
    CHECK(std::isfinite(c.mean_normal[1].x));
}

// ---- Signs -------------------------------------------------------------------------------

TEST_CASE("curvature: torus has K > 0 on the outer equator and K < 0 on the inner", "[curvature]") {
    // Vertex (i=0, j) sits at tube angle v = 2 pi j / nv: j = 0 outer equator, j = nv/2 inner.
    const std::uint32_t nu = 24, nv = 12;
    const auto c = compute_curvature(he(make_torus(nu, nv)));
    CHECK(c.gaussian[0] > 0.0);
    CHECK(c.gaussian[(nv / 2) * nu] < 0.0);
}

TEST_CASE("curvature: sign of H follows the orientation", "[curvature]") {
    auto m = make_icosphere(2);
    const auto outward = compute_curvature(he(m));
    for (auto& t : m.triangles) std::swap(t[1], t[2]);  // wind inward
    const auto inward = compute_curvature(he(m));
    for (std::size_t v = 0; v < m.positions.size(); ++v) {
        CHECK(outward.mean[v] > 0.0);
        CHECK_THAT(inward.mean[v], WithinRel(-outward.mean[v], 1e-12));
        CHECK_THAT(inward.gaussian[v], WithinRel(outward.gaussian[v], 1e-12));  // K is intrinsic
    }
}

// ---- Convergence (measured values recorded in DECISIONS.md D39) -------------------------
//
// Observed order of accuracy p from errors e1, e2 at spacings h, h/2: p = log2(e1 / e2).

namespace {
double order(double coarse, double fine) { return std::log2(coarse / fine); }
}  // namespace

TEST_CASE("curvature: sphere - H exact at s=1 by symmetry, then K and H second order", "[curvature][convergence]") {
    // Identity: if all neighbors lie on a sphere (center c) through x_i, then
    // |x_i - x_j|^2 = 2 (x_i - c).(x_i - x_j), which makes K(x_i).n = 2/r EXACTLY whenever
    // the mixed area is the pure Voronoi area and n is radial. At s = 1 every one-ring is
    // symmetric, so the area-weighted normal is radial; from s = 2 on, the remaining H error
    // is entirely the normal's deviation from radial (D37).
    const double r = 2.0;
    std::vector<double> eh, ek;
    for (std::uint32_t s = 1; s <= 5; ++s) {
        const auto c = compute_curvature(he(make_icosphere(s, r)));
        double k_err = 0.0;
        for (double k : c.gaussian) k_err = std::max(k_err, std::abs(k - 1 / (r * r)));
        eh.push_back(max_mean_error(c, 1 / r));
        ek.push_back(k_err);
    }
    CHECK(eh[0] < 1e-12);
    for (std::size_t i = 1; i + 1 < eh.size(); ++i) CHECK(order(eh[i], eh[i + 1]) > 1.8);
    for (std::size_t i = 0; i + 1 < ek.size(); ++i) CHECK(order(ek[i], ek[i + 1]) > 1.8);
}

TEST_CASE("curvature: cylinder - K = 0 and H = 1/(2r) exactly at every resolution", "[curvature]") {
    // Not convergence but an identity: the faces are flat strips, the cotan weights reproduce
    // linear functions on the unrolled strip, and that balances the horizontal and vertical
    // terms so that H = 1/(2r) up to roundoff. (Errors measured: 1e-16 to 2e-13.)
    const double r = 0.5;
    for (std::uint32_t nu : {12u, 24u, 48u, 96u}) {
        const auto nv = static_cast<std::uint32_t>(std::lround(nu / (2 * pi * r)));
        const auto c = compute_curvature(he(make_cylinder(nu, nv, r, 1.0)));
        for (std::size_t v = 0; v < c.gaussian.size(); ++v) {
            if (c.boundary[v]) continue;
            CHECK_THAT(c.gaussian[v], WithinAbs(0.0, 1e-9));
            CHECK_THAT(c.mean[v], WithinAbs(1 / (2 * r), 1e-10));
            CHECK_THAT(c.k1[v], WithinAbs(1 / r, 1e-9));  // principal curvatures 1/r and 0
            CHECK_THAT(c.k2[v], WithinAbs(0.0, 1e-9));
        }
    }
}

TEST_CASE("curvature: torus - H and K converge at second order to the closed forms", "[curvature][convergence]") {
    // p(u,v) = ((R + r cos v) cos u, (R + r cos v) sin u, r sin v):
    //   K = cos v / (r (R + r cos v)),   H = (R + 2 r cos v) / (2 r (R + r cos v)).
    const double R = 1.0, r = 0.3;
    std::vector<double> eh, ek;
    for (std::uint32_t nv : {8u, 16u, 32u, 64u}) {
        const std::uint32_t nu = 3 * nv;  // roughly square cells for R/r ~ 3
        const auto c = compute_curvature(he(make_torus(nu, nv, R, r)));
        double h_err = 0.0, k_err = 0.0;
        for (std::uint32_t j = 0; j < nv; ++j) {
            const double v = 2 * pi * j / nv, ring = R + r * std::cos(v);
            for (std::uint32_t i = 0; i < nu; ++i) {
                h_err = std::max(h_err, std::abs(c.mean[j * nu + i] - (R + 2 * r * std::cos(v)) / (2 * r * ring)));
                k_err = std::max(k_err, std::abs(c.gaussian[j * nu + i] - std::cos(v) / (r * ring)));
            }
        }
        eh.push_back(h_err);
        ek.push_back(k_err);
    }
    for (std::size_t i = 0; i + 1 < eh.size(); ++i) {
        CHECK(order(eh[i], eh[i + 1]) > 1.8);
        CHECK(order(ek[i], ek[i + 1]) > 1.8);
    }
}

TEST_CASE("curvature: KNOWN LIMITATION - no pointwise convergence on irregular meshes", "[curvature][limitation]") {
    // Characterization test, not a goal: perturb icosphere vertices along the sphere by
    // ~0.25 h. The max H error then stops shrinking under refinement (measured: ~0.28 at
    // s = 4, 5, 6). Kept so M6b (intrinsic Delaunay) can show a measured improvement.
    // Deterministic hash jitter: <random> distributions differ between standard libraries.
    auto jitter = [](double k) { const double s = std::sin(k * 12.9898) * 43758.5453; return 2 * (s - std::floor(s)) - 1; };
    auto m = make_icosphere(5, 1.0);
    const double h = 1.1 / 32.0;  // ~ edge length at s = 5
    for (std::size_t v = 0; v < m.positions.size(); ++v) {
        auto& p = m.positions[v];
        const double k = static_cast<double>(v);
        p.x += 0.25 * h * jitter(3 * k);
        p.y += 0.25 * h * jitter(3 * k + 1);
        p.z += 0.25 * h * jitter(3 * k + 2);
        const double n = std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        p = {p.x / n, p.y / n, p.z / n};
    }
    CHECK(max_mean_error(compute_curvature(he(m)), 1.0) > 0.05);
}
