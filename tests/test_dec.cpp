#include <cmath>
#include <cstdint>
#include <numeric>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/curvature.h"
#include "core/dec.h"
#include "core/generate.h"

using namespace dmw;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

HalfEdgeMesh he(const TriMesh& m) {
    auto r = build_halfedge(m);
    REQUIRE(r.ok());
    return std::move(r.mesh);
}

// Deterministic pseudo-random values in [-1, 1] (portable, unlike <random> distributions).
double hash01(double k) {
    const double s = std::sin(k * 12.9898) * 43758.5453;
    return 2 * (s - std::floor(s)) - 1;
}
std::vector<double> pseudo_random(std::size_t n, double seed) {
    std::vector<double> v(n);
    for (std::size_t i = 0; i < n; ++i) v[i] = hash01(seed + static_cast<double>(i));
    return v;
}

// Flat 8x8 grid with interior vertices jittered within the plane (irregular, still planar).
TriMesh jittered_grid() {
    auto m = make_grid(8, 8);
    for (std::size_t v = 0; v < m.positions.size(); ++v) {
        auto& p = m.positions[v];
        if (p.x > 0 && p.x < 1 && p.y > 0 && p.y < 1) {
            p.x += 0.03 * hash01(2.0 * static_cast<double>(v));
            p.y += 0.03 * hash01(2.0 * static_cast<double>(v) + 1);
        }
    }
    return m;
}

}  // namespace

TEST_CASE("dec: d1 d0 = 0 exactly (curl of a gradient vanishes)", "[dec]") {
    for (const TriMesh& m : {make_icosphere(2), make_torus(10, 6), make_grid_with_holes(6, 6, {{2, 2}})}) {
        const auto ops = build_dec(he(m));
        for (double seed : {1.0, 2.0, 3.0}) {
            const auto f = pseudo_random(m.positions.size(), seed * 1000);
            const auto curl_grad = ops.d1.multiply(ops.d0.multiply(f));
            // Entries are +-1, so each face sums the same three differences with opposite
            // signs: the cancellation is exact even in floating point.
            for (double c : curl_grad) CHECK(c == 0.0);
        }
    }
}

TEST_CASE("dec: d0 kills constants, so L kills constants (H^0 = constants)", "[dec]") {
    const auto m = make_torus(10, 6);
    const auto ops = build_dec(he(m));
    const std::vector<double> ones(m.positions.size(), 1.0);
    for (double e : ops.d0.multiply(ones)) CHECK(e == 0.0);
    for (double l : ops.laplacian.multiply(ones)) CHECK_THAT(l, WithinAbs(0.0, 1e-12));
}

TEST_CASE("dec: L is symmetric negative semidefinite", "[dec]") {
    const auto ops = build_dec(he(make_icosphere(2)));
    const auto& L = ops.laplacian;
    for (std::uint32_t r = 0; r < L.rows; ++r) {
        for (std::uint32_t k = L.row_start[r]; k < L.row_start[r + 1]; ++k) CHECK(L.at(L.col[k], r) == L.value[k]);
    }
    for (double seed : {1.0, 7.0, 42.0}) {
        const auto u = pseudo_random(L.rows, seed * 100);
        const auto lu = L.multiply(u);
        CHECK(std::inner_product(u.begin(), u.end(), lu.begin(), 0.0) <= 0.0);  // = -sum *1 (d0 u)^2
    }
}

TEST_CASE("dec: *0 (lumped areas) sums to the surface area; counts match V, E, F", "[dec]") {
    const auto m = make_torus(12, 6);
    const auto ops = build_dec(he(m));
    CHECK(ops.d0.rows == ops.edges.size());
    CHECK(ops.d0.cols == m.positions.size());
    CHECK(ops.d1.rows == m.triangles.size());
    CHECK(ops.d1.cols == ops.edges.size());
    const auto curv = compute_curvature(he(m));  // mixed areas also tile the surface
    const double area = std::accumulate(curv.mixed_area.begin(), curv.mixed_area.end(), 0.0);
    CHECK_THAT(std::accumulate(ops.star0.begin(), ops.star0.end(), 0.0), WithinRel(area, 1e-12));
}

TEST_CASE("dec: obtuse angle gives a negative *1 weight (the cotan caveat)", "[dec]") {
    const auto ops = build_dec(he({{{0, 0, 0}, {4, 0, 0}, {2, 0.5, 0}}, {{0, 1, 2}}}));
    for (std::size_t e = 0; e < ops.edges.size(); ++e) {
        if (ops.edges[e].v0 == 0 && ops.edges[e].v1 == 1) CHECK(ops.star1[e] < 0.0);  // opposite the obtuse corner
    }
}

TEST_CASE("dec: L reproduces linear functions on flat irregular meshes (linear precision)", "[dec]") {
    const auto m = jittered_grid();
    const auto mesh = he(m);
    const auto ops = build_dec(mesh);
    std::vector<double> x(m.positions.size()), y(m.positions.size());
    for (std::size_t v = 0; v < x.size(); ++v) x[v] = m.positions[v].x, y[v] = m.positions[v].y;
    const auto lx = ops.laplacian.multiply(x), ly = ops.laplacian.multiply(y);
    const auto curv = compute_curvature(mesh);
    for (std::size_t v = 0; v < x.size(); ++v) {
        if (curv.boundary[v]) continue;
        CHECK_THAT(lx[v], WithinAbs(0.0, 1e-12));
        CHECK_THAT(ly[v], WithinAbs(0.0, 1e-12));
    }
}

TEST_CASE("dec: L x = -A K(x) ties the Laplacian to M4's mean-curvature normal", "[dec]") {
    // (L x)_i = 1/2 sum (cot a + cot b)(x_j - x_i) = -A_i K(x_i)  (Meyer et al.'s formula).
    const auto mesh = he(make_torus(16, 8));
    const auto ops = build_dec(mesh);
    const auto curv = compute_curvature(mesh);
    const std::size_t n = mesh.positions.size();
    std::vector<double> coord(n);
    for (int axis = 0; axis < 3; ++axis) {
        for (std::size_t v = 0; v < n; ++v) {
            const Vec3& p = mesh.positions[v];
            coord[v] = axis == 0 ? p.x : axis == 1 ? p.y : p.z;
        }
        const auto l = ops.laplacian.multiply(coord);
        for (std::size_t v = 0; v < n; ++v) {
            const Vec3& k = curv.mean_normal[v];
            const double kc = axis == 0 ? k.x : axis == 1 ? k.y : k.z;
            CHECK_THAT(l[v], WithinAbs(-curv.mixed_area[v] * kc, 1e-12));
        }
    }
}

TEST_CASE("dec: gradient of a linear function is exact on every face", "[dec]") {
    const auto m = jittered_grid();
    std::vector<double> f(m.positions.size());
    for (std::size_t v = 0; v < f.size(); ++v) f[v] = 2.0 * m.positions[v].x - 3.0 * m.positions[v].y + 5.0;
    for (const Vec3& g : face_gradient(he(m), f)) {
        CHECK_THAT(g.x, WithinAbs(2.0, 1e-12));
        CHECK_THAT(g.y, WithinAbs(-3.0, 1e-12));
        CHECK_THAT(g.z, WithinAbs(0.0, 1e-12));
    }
}

TEST_CASE("dec: div(grad u) = L u exactly, for any u, including boundary vertices", "[dec]") {
    for (const TriMesh& m : {make_icosphere(2), make_grid_with_holes(6, 6, {{2, 2}}), jittered_grid()}) {
        const auto mesh = he(m);
        const auto u = pseudo_random(m.positions.size(), 3.0);
        const auto lu = build_dec(mesh).laplacian.multiply(u);
        const auto dg = vertex_divergence(mesh, face_gradient(mesh, u));
        for (std::size_t v = 0; v < u.size(); ++v) CHECK_THAT(dg[v], WithinAbs(lu[v], 1e-10));
    }
}
