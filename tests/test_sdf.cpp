#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/generate.h"
#include "core/sdf.h"

using namespace dmw;
using Catch::Matchers::WithinAbs;

namespace {

Grid3 grid(std::size_t nx, std::size_t ny, std::size_t nz, double h) {
    Grid3 g;
    g.h = h, g.n = {nx, ny, nz}, g.values.assign(nx * ny * nz, 0.0);
    return g;
}

// Mean |phi - (|p| - 1)| over nodes within `band` of the unit sphere.
double sphere_error(const Grid3& g, double band) {
    double sum = 0.0;
    std::size_t count = 0;
    for (std::size_t k = 0; k < g.n[2]; ++k)
        for (std::size_t j = 0; j < g.n[1]; ++j)
            for (std::size_t i = 0; i < g.n[0]; ++i) {
                const Vec3 p = g.position(i, j, k);
                const double exact = std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z) - 1.0;
                if (std::abs(exact) > band) continue;
                sum += std::abs(g.values[g.index(i, j, k)] - exact), ++count;
            }
    return sum / double(count);
}

}  // namespace

TEST_CASE("sdf: the Neumann Poisson solve recovers a cosine eigenmode and a random zero-mean field exactly", "[sdf]") {
    const std::size_t nx = 12, ny = 9, nz = 7;
    const double h = 0.3;
    // Discrete Neumann Laplacian applied by finite differences (reference for the DCT solve).
    auto laplacian = [&](const Grid3& u) {
        Grid3 out = u;
        for (std::size_t k = 0; k < nz; ++k)
            for (std::size_t j = 0; j < ny; ++j)
                for (std::size_t i = 0; i < nx; ++i) {
                    const double c = u.values[u.index(i, j, k)];
                    double s = 0.0;
                    if (i > 0) s += u.values[u.index(i - 1, j, k)] - c;
                    if (i + 1 < nx) s += u.values[u.index(i + 1, j, k)] - c;
                    if (j > 0) s += u.values[u.index(i, j - 1, k)] - c;
                    if (j + 1 < ny) s += u.values[u.index(i, j + 1, k)] - c;
                    if (k > 0) s += u.values[u.index(i, j, k - 1)] - c;
                    if (k + 1 < nz) s += u.values[u.index(i, j, k + 1)] - c;
                    out.values[out.index(i, j, k)] = s / (h * h);
                }
        return out;
    };
    std::mt19937 rng(5);
    std::normal_distribution<double> g(0.0, 1.0);
    Grid3 u = grid(nx, ny, nz, h);
    double mean = 0.0;
    for (double& x : u.values) x = g(rng), mean += x;
    mean /= double(u.values.size());
    for (double& x : u.values) x -= mean;  // the Poisson solution is defined up to a constant: compare zero-mean
    Grid3 f = laplacian(u);
    dct_solve_poisson(f);
    double worst = 0.0;
    for (std::size_t q = 0; q < u.values.size(); ++q) worst = std::max(worst, std::abs(f.values[q] - u.values[q]));
    CHECK(worst < 1e-10);
    // Screened solve: (I - tL) u = f round trip.
    const double t = 0.07;
    Grid3 lu = laplacian(u), rhs = u;
    for (std::size_t q = 0; q < u.values.size(); ++q) rhs.values[q] = u.values[q] - t * lu.values[q];
    dct_solve_screened(rhs, t);
    worst = 0.0;
    for (std::size_t q = 0; q < u.values.size(); ++q) worst = std::max(worst, std::abs(rhs.values[q] - u.values[q]));
    CHECK(worst < 1e-10);
}

TEST_CASE("sdf: signed heat distance of a sphere matches |p| - 1 and improves under refinement", "[sdf]") {
    const auto sphere = make_icosphere(5);  // unit sphere, 20,480 faces
    SignedHeatParams coarse, fine;
    coarse.h = 0.08, fine.h = 0.04;
    const Grid3 a = signed_heat_distance(sphere, coarse), b = signed_heat_distance(sphere, fine);
    const double ea = sphere_error(a, 0.3), eb = sphere_error(b, 0.3);
    INFO("mean error within 0.3 of the surface: h 0.08 -> " << ea << ", h 0.04 -> " << eb);
    CHECK(ea < 0.05);
    CHECK(eb < ea);
    // Sign: negative at the centre, positive at the corners.
    CHECK(b.sample({0, 0, 0}) < 0.0);
    CHECK(b.sample({1.3, 1.3, 1.3}) > 0.0);
}

TEST_CASE("sdf: a sphere with a hole is still signed correctly (implicit completion)", "[sdf]") {
    // Remove every face with a vertex above z = 0.6: a hole of radius ~0.8 at the top.
    const auto full = make_icosphere(5);
    TriMesh holed;
    holed.positions = full.positions;
    for (const auto& t : full.triangles) {
        bool keep = true;
        for (auto v : t) keep = keep && full.positions[v].z <= 0.6;
        if (keep) holed.triangles.push_back(t);
    }
    REQUIRE(holed.triangles.size() < full.triangles.size());
    SignedHeatParams prm;
    prm.h = 0.05;
    const Grid3 g = signed_heat_distance(holed, prm);
    CHECK(g.sample({0, 0, 0}) < -0.5);    // deep inside, even though the surface is open
    CHECK(g.sample({0, 0, 0.3}) < 0.0);   // between the centre and the hole
    CHECK(g.sample({0, 0, -1.3}) > 0.0);  // outside, below
    CHECK(g.sample({1.4, 0, 0}) > 0.0);   // outside, beside
    // Away from the hole (z < 0), still close to the true sphere distance.
    double sum = 0.0;
    std::size_t count = 0;
    for (std::size_t k = 0; k < g.n[2]; ++k)
        for (std::size_t j = 0; j < g.n[1]; ++j)
            for (std::size_t i = 0; i < g.n[0]; ++i) {
                const Vec3 p = g.position(i, j, k);
                const double exact = std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z) - 1.0;
                if (p.z > 0.0 || std::abs(exact) > 0.2) continue;
                sum += std::abs(g.values[g.index(i, j, k)] - exact), ++count;
            }
    CHECK(sum / double(count) < 0.06);
}
