#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/cholesky.h"
#include "core/curvature.h"
#include "core/dec.h"
#include "core/geodesic.h"
#include "core/generate.h"
#include "core/intrinsic.h"

using namespace dmw;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace {

HalfEdgeMesh he(const TriMesh& m) {
    auto r = build_halfedge(m);
    REQUIRE(r.ok());
    return std::move(r.mesh);
}

double hash01(double k) {
    const double s = std::sin(k * 12.9898) * 43758.5453;
    return 2 * (s - std::floor(s)) - 1;
}

// Icosphere with vertices pushed ~amp*h along the sphere: irregular, many obtuse angles.
TriMesh jittered_sphere(std::uint32_t s, double amp) {
    auto m = make_icosphere(s, 1.0);
    const double h = 1.1 / static_cast<double>(1u << s);
    for (std::size_t v = 1; v < m.positions.size(); ++v) {
        auto& p = m.positions[v];
        const double k = static_cast<double>(v);
        p.x += amp * h * hash01(3 * k), p.y += amp * h * hash01(3 * k + 1), p.z += amp * h * hash01(3 * k + 2);
        const double n = std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        p = {p.x / n, p.y / n, p.z / n};
    }
    return m;
}

// Flat grid whose interior vertices are sheared, producing long obtuse triangles.
TriMesh sheared_grid() {
    auto m = make_grid(10, 10);
    for (auto& p : m.positions) {
        if (p.x > 0 && p.x < 1 && p.y > 0 && p.y < 1) p.x += 0.6 * (p.y - 0.5) * 0.09 * std::sin(40 * p.y);
    }
    return m;
}

double min_weight(const DecOperators& ops) { return *std::min_element(ops.star1.begin(), ops.star1.end()); }

std::size_t count_edges(const IntrinsicTriangulation& t) {
    std::size_t n = 0;
    for (std::uint32_t h = 0; h < t.connectivity.twin.size(); ++h) {
        if (t.connectivity.twin[h] == kInvalid || h < t.connectivity.twin[h]) ++n;
    }
    return n;
}

}  // namespace

TEST_CASE("intrinsic: Euclidean lengths reproduce the extrinsic DEC operators", "[intrinsic]") {
    // Same Laplacian and mass matrix from lengths alone as from 3D positions: the cotan
    // Laplacian is intrinsic.
    for (const TriMesh& m : {make_torus(10, 6), make_grid_with_holes(6, 6, {{2, 2}}), jittered_sphere(2, 0.25)}) {
        const auto mesh = he(m);
        const auto ext = build_dec(mesh);
        const auto in = build_dec(intrinsic_from_mesh(mesh));
        REQUIRE(in.laplacian.nonzeros() == ext.laplacian.nonzeros());
        for (std::size_t k = 0; k < ext.laplacian.value.size(); ++k) {
            CHECK(in.laplacian.col[k] == ext.laplacian.col[k]);
            CHECK_THAT(in.laplacian.value[k], WithinAbs(ext.laplacian.value[k], 1e-12));
        }
        for (std::size_t v = 0; v < ext.star0.size(); ++v) CHECK_THAT(in.star0[v], WithinAbs(ext.star0[v], 1e-14));
    }
}

TEST_CASE("intrinsic: div(grad u) = L u from lengths alone", "[intrinsic]") {
    const auto t = intrinsic_delaunay(he(jittered_sphere(2, 0.25)));
    std::vector<double> u(t.connectivity.positions.size());
    for (std::size_t v = 0; v < u.size(); ++v) u[v] = hash01(static_cast<double>(v) + 0.5);
    const auto lu = build_dec(t).laplacian.multiply(u);
    const auto dg = intrinsic_divergence(t, intrinsic_gradient(t, u));
    for (std::size_t v = 0; v < u.size(); ++v) CHECK_THAT(dg[v], WithinAbs(lu[v], 1e-10));
}

TEST_CASE("intrinsic: regular icosphere is already Delaunay (no flips)", "[intrinsic]") {
    const auto t = intrinsic_delaunay(he(make_icosphere(3)));
    CHECK(t.flips == 0);
    CHECK(min_weight(build_dec(t)) >= 0.0);
}

TEST_CASE("intrinsic: flipping makes every cotan weight non-negative", "[intrinsic]") {
    for (const TriMesh& m : {jittered_sphere(3, 0.25), sheared_grid()}) {
        const auto mesh = he(m);
        const auto before = build_dec(intrinsic_from_mesh(mesh));
        const auto t = intrinsic_delaunay(mesh);
        CHECK(min_weight(before) < 0.0);  // the input really had obtuse-pair edges
        CHECK(t.flips > 0);
        CHECK(t.skipped == 0);
        CHECK(min_weight(build_dec(t)) >= -1e-12);
        for (std::uint32_t h = 0; h < t.length.size(); ++h) CHECK(is_delaunay(t, h));
    }
}

TEST_CASE("intrinsic: flips preserve the surface (area, cone angles, edge count)", "[intrinsic]") {
    const auto mesh = he(jittered_sphere(3, 0.25));
    const auto before = intrinsic_from_mesh(mesh);
    const auto after = intrinsic_delaunay(mesh);
    const auto a0 = build_dec(before).star0, a1 = build_dec(after).star0;
    CHECK_THAT(std::accumulate(a1.begin(), a1.end(), 0.0), WithinRel(std::accumulate(a0.begin(), a0.end(), 0.0), 1e-12));
    const auto s0 = angle_sums(before), s1 = angle_sums(after);
    for (std::size_t v = 0; v < s0.size(); ++v) CHECK_THAT(s1[v], WithinAbs(s0[v], 1e-10));  // K is intrinsic
    CHECK(count_edges(after) == count_edges(before));  // a flip replaces one edge by one edge
}

TEST_CASE("intrinsic: the flipped connectivity is a valid half-edge mesh", "[intrinsic]") {
    // Rebuild from the flipped faces: build_halfedge re-checks every invariant from scratch.
    const auto t = intrinsic_delaunay(he(sheared_grid()));
    TriMesh faces{t.connectivity.positions, {}};
    for (std::size_t f = 0; f < t.connectivity.origin.size() / 3; ++f) {
        faces.triangles.push_back({t.connectivity.origin[3 * f], t.connectivity.origin[3 * f + 1], t.connectivity.origin[3 * f + 2]});
    }
    const auto rebuilt = build_halfedge(faces);
    REQUIRE(rebuilt.ok());
    CHECK(rebuilt.mesh.twin == t.connectivity.twin);
    for (std::uint32_t v = 0; v < faces.positions.size(); ++v) {
        CHECK((t.connectivity.vertex_halfedge[v] == kInvalid) == (rebuilt.mesh.vertex_halfedge[v] == kInvalid));
        if (t.connectivity.vertex_halfedge[v] != kInvalid) {
            CHECK(t.connectivity.origin[t.connectivity.vertex_halfedge[v]] == v);
        }
    }
}

TEST_CASE("intrinsic: on a flat mesh, flipping keeps linear precision", "[intrinsic]") {
    // Flat: intrinsic edges are straight segments in the plane, so L x = 0 at interior vertices.
    const auto mesh = he(sheared_grid());
    const auto t = intrinsic_delaunay(mesh);
    const auto L = build_dec(t).laplacian;
    std::vector<double> x(mesh.positions.size());
    for (std::size_t v = 0; v < x.size(); ++v) x[v] = mesh.positions[v].x;
    const auto lx = L.multiply(x);
    for (std::size_t v = 0; v < x.size(); ++v) {
        const auto& p = mesh.positions[v];
        if (p.x > 0 && p.x < 1 && p.y > 0 && p.y < 1) CHECK_THAT(lx[v], WithinAbs(0.0, 1e-10));
    }
}

// ---- Measured effect (DECISIONS.md D55, thresholds with margin) ---------------------------

namespace {

std::vector<double> heat_from(const IntrinsicTriangulation& t, std::uint32_t src) {
    const auto ops = build_dec(t);
    const double tt = ops.mean_edge_length * ops.mean_edge_length;
    const EnvelopeLdlt f(scaled_plus_diagonal(ops.laplacian, -tt, ops.star0));
    REQUIRE(f.ok());
    std::vector<double> d(ops.star0.size(), 0.0);
    d[src] = 1.0;
    return f.solve(d);
}

}  // namespace

TEST_CASE("intrinsic: restores the maximum principle (no negative heat)", "[intrinsic]") {
    // Heat from a point source can't go negative. With negative cotan weights it does.
    const std::uint32_t n = 20, src = (n / 2) * (n + 1) + n / 2;
    const auto mesh = he(make_brick_grid(n, 0.45));
    const auto cot = heat_from(intrinsic_from_mesh(mesh), src);
    const auto idt = heat_from(intrinsic_delaunay(mesh), src);
    CHECK(*std::min_element(cot.begin(), cot.end()) < 0.0);   // measured: 54 negative values
    CHECK(*std::min_element(idt.begin(), idt.end()) > 0.0);   // measured: none
}

TEST_CASE("intrinsic: heat geodesics converge on an obtuse mesh where cotan fails", "[intrinsic][convergence]") {
    std::vector<double> cot_err, idt_err;
    for (std::uint32_t n : {20u, 40u}) {
        const std::uint32_t src = (n / 2) * (n + 1) + n / 2;
        const auto mesh = he(make_brick_grid(n, 0.45));
        for (bool idt : {false, true}) {
            const auto r = HeatGeodesics(mesh, 1.0, idt).distance(std::vector<std::uint32_t>{src});
            REQUIRE(r.ok());
            double mean = 0.0;
            for (std::size_t v = 0; v < r.distance.size(); ++v) {
                const auto &a = mesh.positions[src], &b = mesh.positions[v];
                mean += std::abs(r.distance[v] - std::hypot(b.x - a.x, b.y - a.y));
            }
            (idt ? idt_err : cot_err).push_back(mean / static_cast<double>(r.distance.size()));
        }
    }
    CHECK(cot_err[1] > cot_err[0]);  // cotan gets WORSE under refinement (measured 4.3e-2 -> 1.07e-1)
    CHECK(idt_err[1] < idt_err[0]);  // intrinsic Delaunay converges (measured 9.9e-3 -> 5.9e-3)
    CHECK(idt_err[1] < 1e-2);
    CHECK(idt_err[1] < 0.2 * cot_err[1]);
}

TEST_CASE("intrinsic: mean curvature blow-ups on a badly jittered sphere are reduced", "[intrinsic]") {
    // Not a convergence claim: max |H - 1| still grows with refinement (D55). It drops ~8x.
    const auto mesh = he(jittered_sphere(4, 0.45));
    auto max_err = [](const CurvatureField& c) {
        double e = 0.0;
        for (double h : c.mean) e = std::max(e, std::abs(h - 1.0));
        return e;
    };
    const auto cot = compute_curvature(mesh), idt = compute_curvature(mesh, true);
    CHECK(idt.intrinsic_flips > 0);
    CHECK(max_err(idt) < 0.25 * max_err(cot));
    for (std::size_t v = 0; v < cot.angle_defect.size(); ++v) {
        CHECK_THAT(idt.angle_defect[v], WithinAbs(cot.angle_defect[v], 1e-10));  // K untouched
    }
}

TEST_CASE("intrinsic: a flip that would duplicate an edge is skipped, not performed", "[intrinsic]") {
    // Flattened tetrahedron: O and Z sit near the midpoint of edge X-Y, so the angles opposite
    // X-Y are ~169 degrees each (non-Delaunay). Flipping X-Y would create O-Z, which already
    // exists (a tetrahedron's graph is K4): not representable, so it must be skipped (D53).
    const TriMesh m{{{0, 0.1, 0}, {-1, 0, 0}, {1, 0, 0}, {0, 0, 0.1}},  // O, X, Y, Z
                    {{0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 2, 3}}};
    const auto mesh = he(m);
    std::uint32_t xy = kInvalid;
    for (std::uint32_t h = 0; h < mesh.origin.size(); ++h) {
        if (mesh.origin[h] == 1 && dest(mesh, h) == 2) xy = h;
    }
    REQUIRE(xy != kInvalid);
    REQUIRE_FALSE(is_delaunay(intrinsic_from_mesh(mesh), xy));
    const auto t = intrinsic_delaunay(mesh);
    CHECK(t.flips == 0);
    CHECK(t.skipped >= 1);
    CHECK(t.connectivity.twin == mesh.twin);  // untouched
}
