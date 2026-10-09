#include "core/sdf.h"

#include <algorithm>
#include <cmath>
#include <functional>
#ifndef __EMSCRIPTEN__
#include <thread>
#endif

#include "detail/vec.h"

namespace dmw {
namespace {

using namespace detail;

// Run body(begin, end) over [0, count) on hardware threads (native) or serially (WebAssembly).
void parallel_for(std::size_t count, const std::function<void(std::size_t, std::size_t)>& body) {
#ifdef __EMSCRIPTEN__
    body(0, count);
#else
    const std::size_t t = count < 64 ? 1 : std::max<std::size_t>(1, std::thread::hardware_concurrency());
    if (t == 1) return body(0, count);
    std::vector<std::thread> pool;
    const std::size_t chunk = (count + t - 1) / t;
    for (std::size_t b = 0; b < count; b += chunk) pool.emplace_back(body, b, std::min(count, b + chunk));
    for (auto& th : pool) th.join();
#endif
}

// Orthonormal DCT-II matrix: C[k][i] = s_k cos(pi k (i + 1/2) / n); its inverse is the transpose.
std::vector<double> dct_matrix(std::size_t n) {
    std::vector<double> c(n * n);
    const double pi = std::acos(-1.0);
    for (std::size_t k = 0; k < n; ++k)
        for (std::size_t i = 0; i < n; ++i)
            c[k * n + i] = (k == 0 ? std::sqrt(1.0 / double(n)) : std::sqrt(2.0 / double(n))) *
                           std::cos(pi * double(k) * (double(i) + 0.5) / double(n));
    return c;
}

// Apply the DCT (forward) or its inverse along one axis of the grid, line by line.
void dct_axis(Grid3& g, int axis, bool forward) {
    const std::size_t n = g.n[static_cast<std::size_t>(axis)];
    const std::vector<double> c = dct_matrix(n);
    const std::size_t stride = axis == 0 ? 1 : axis == 1 ? g.n[0] : g.n[0] * g.n[1];
    const std::size_t lines = g.values.size() / n;
    parallel_for(lines, [&](std::size_t begin, std::size_t end) {
        std::vector<double> in(n), out(n);
        for (std::size_t line = begin; line < end; ++line) {
            // First element of this line: decompose `line` into the two other coordinates.
            std::size_t base;
            if (axis == 0) base = line * g.n[0];
            else if (axis == 1) base = (line / g.n[0]) * g.n[0] * g.n[1] + line % g.n[0];
            else base = line;
            for (std::size_t i = 0; i < n; ++i) in[i] = g.values[base + i * stride];
            for (std::size_t k = 0; k < n; ++k) {
                double s = 0.0;
                if (forward)
                    for (std::size_t i = 0; i < n; ++i) s += c[k * n + i] * in[i];
                else
                    for (std::size_t i = 0; i < n; ++i) s += c[i * n + k] * in[i];
                out[k] = s;
            }
            for (std::size_t i = 0; i < n; ++i) g.values[base + i * stride] = out[i];
        }
    });
}

// Multiply each DCT coefficient by scale(sum of the three axis eigenvalues of -L).
void dct_diagonal_solve(Grid3& g, const std::function<double(double)>& scale) {
    for (int a = 0; a < 3; ++a) dct_axis(g, a, true);
    std::array<std::vector<double>, 3> lam;
    for (std::size_t a = 0; a < 3; ++a) {
        lam[a].resize(g.n[a]);
        for (std::size_t k = 0; k < g.n[a]; ++k)
            lam[a][k] = (2.0 - 2.0 * std::cos(std::acos(-1.0) * double(k) / double(g.n[a]))) / (g.h * g.h);
    }
    for (std::size_t k = 0; k < g.n[2]; ++k)
        for (std::size_t j = 0; j < g.n[1]; ++j)
            for (std::size_t i = 0; i < g.n[0]; ++i) g.values[g.index(i, j, k)] *= scale(lam[0][i] + lam[1][j] + lam[2][k]);
    for (int a = 0; a < 3; ++a) dct_axis(g, a, false);
}

}  // namespace

Vec3 Grid3::position(std::size_t i, std::size_t j, std::size_t k) const {
    return {origin.x + h * double(i), origin.y + h * double(j), origin.z + h * double(k)};
}

double Grid3::sample(const Vec3& p) const {
    const double fx = std::clamp((p.x - origin.x) / h, 0.0, double(n[0] - 1));
    const double fy = std::clamp((p.y - origin.y) / h, 0.0, double(n[1] - 1));
    const double fz = std::clamp((p.z - origin.z) / h, 0.0, double(n[2] - 1));
    const std::size_t i = std::min(static_cast<std::size_t>(fx), n[0] - 2), j = std::min(static_cast<std::size_t>(fy), n[1] - 2),
                      k = std::min(static_cast<std::size_t>(fz), n[2] - 2);
    const double u = fx - double(i), v = fy - double(j), w = fz - double(k);
    double s = 0.0;
    for (int c = 0; c < 8; ++c) {
        const std::size_t di = c & 1, dj = (c >> 1) & 1, dk = (c >> 2) & 1;
        s += (di ? u : 1 - u) * (dj ? v : 1 - v) * (dk ? w : 1 - w) * values[index(i + di, j + dj, k + dk)];
    }
    return s;
}

void dct_solve_screened(Grid3& field, double t) {
    dct_diagonal_solve(field, [t](double lam) { return 1.0 / (1.0 + t * lam); });
}

void dct_solve_poisson(Grid3& field) {
    dct_diagonal_solve(field, [](double lam) { return lam > 0.0 ? -1.0 / lam : 0.0; });  // L u = f, L = -(sum lam)
}

Grid3 signed_heat_distance(const TriMesh& m, const SignedHeatParams& prm) {
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (const auto& t : m.triangles)
        for (auto v : t) {
            const Vec3& p = m.positions[v];
            lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
            hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
        }
    const double h = prm.h > 0.0 ? prm.h : norm(hi - lo) / 64.0;
    const double pad = prm.padding > 0.0 ? prm.padding : 8.0 * h;
    Grid3 g;
    g.h = h;
    g.origin = lo - Vec3{pad, pad, pad};
    const Vec3 span = hi - lo;
    g.n = {static_cast<std::size_t>(std::ceil((span.x + 2 * pad) / h)) + 1, static_cast<std::size_t>(std::ceil((span.y + 2 * pad) / h)) + 1,
           static_cast<std::size_t>(std::ceil((span.z + 2 * pad) / h)) + 1};
    const std::size_t total = g.n[0] * g.n[1] * g.n[2];
    g.values.assign(total, 0.0);

    // --- Source: area-weighted unit normals at sample points on each triangle (spacing <= h / 2), splatted
    // trilinearly as a density (divide by the cell volume h^3). Samples are kept for the final shift.
    std::array<Grid3, 3> y{g, g, g};
    struct Sample {
        Vec3 p;
        double area;
    };
    std::vector<Sample> samples;
    for (const auto& t : m.triangles) {
        const Vec3 &a = m.positions[t[0]], &b = m.positions[t[1]], &c = m.positions[t[2]];
        const Vec3 n2 = cross(b - a, c - a);
        const double area = 0.5 * norm(n2);
        if (area <= 0.0) continue;
        const Vec3 nrm = (1.0 / (2.0 * area)) * n2;
        const double longest = std::max({norm(b - a), norm(c - b), norm(a - c)});
        const auto sub = static_cast<std::size_t>(std::max(1.0, std::ceil(longest / (0.5 * h))));
        const double sub_area = area / double(sub * sub);
        // Sub-triangles of a regular subdivision: centroids of the "up" and "down" triangles.
        for (std::size_t r = 0; r < sub; ++r)
            for (std::size_t q = 0; q + r < sub; ++q) {
                for (int down = 0; down < 2; ++down) {
                    if (down && q + r + 1 >= sub) continue;
                    const double u0 = down ? (double(q) + 2.0 / 3.0) / double(sub) : (double(q) + 1.0 / 3.0) / double(sub);
                    const double v0 = down ? (double(r) + 2.0 / 3.0) / double(sub) : (double(r) + 1.0 / 3.0) / double(sub);
                    const Vec3 p = a + u0 * (b - a) + v0 * (c - a);
                    samples.push_back({p, sub_area});
                    const double fx = (p.x - g.origin.x) / h, fy = (p.y - g.origin.y) / h, fz = (p.z - g.origin.z) / h;
                    const auto i = static_cast<std::size_t>(fx), j = static_cast<std::size_t>(fy), k = static_cast<std::size_t>(fz);
                    const double u = fx - double(i), v = fy - double(j), w = fz - double(k);
                    for (int cc = 0; cc < 8; ++cc) {
                        const std::size_t di = cc & 1, dj = (cc >> 1) & 1, dk = (cc >> 2) & 1;
                        const double wt = (di ? u : 1 - u) * (dj ? v : 1 - v) * (dk ? w : 1 - w) * sub_area / (h * h * h);
                        const std::size_t id = g.index(i + di, j + dj, k + dk);
                        y[0].values[id] += wt * nrm.x, y[1].values[id] += wt * nrm.y, y[2].values[id] += wt * nrm.z;
                    }
                }
            }
    }
    // --- 1. Diffuse: (I - t L) Y = N, componentwise (Euclidean domain).
    const double t = prm.t_coef * h * h;
    for (auto& comp : y) dct_solve_screened(comp, t);
    // --- 2. Normalize.
    for (std::size_t q = 0; q < total; ++q) {
        const double len = std::sqrt(y[0].values[q] * y[0].values[q] + y[1].values[q] * y[1].values[q] + y[2].values[q] * y[2].values[q]);
        const double s = len > 1e-300 ? 1.0 / len : 0.0;
        y[0].values[q] *= s, y[1].values[q] *= s, y[2].values[q] *= s;
    }
    // --- 3. Integrate: L phi = -D^T X_edge, X_edge = mean of the two node values along the edge's axis.
    g.values.assign(total, 0.0);
    for (std::size_t a = 0; a < 3; ++a) {
        const std::size_t stride = a == 0 ? 1 : a == 1 ? g.n[0] : g.n[0] * g.n[1];
        for (std::size_t k = 0; k < g.n[2]; ++k)
            for (std::size_t j = 0; j < g.n[1]; ++j)
                for (std::size_t i = 0; i < g.n[0]; ++i) {
                    const std::size_t c = g.index(i, j, k), coord = a == 0 ? i : a == 1 ? j : k;
                    if (coord + 1 >= g.n[a]) continue;
                    const double xe = 0.5 * (y[a].values[c] + y[a].values[c + stride]);  // edge c -> c + stride
                    // D^T: the edge contributes -xe/h to its start node and +xe/h to its end node; f = -D^T X.
                    g.values[c] += xe / h;
                    g.values[c + stride] -= xe / h;
                }
    }
    dct_solve_poisson(g);
    // --- 4. Shift: area-weighted mean over the input surface = 0.
    double num = 0.0, den = 0.0;
    for (const Sample& s : samples) num += s.area * g.sample(s.p), den += s.area;
    const double shift = den > 0.0 ? num / den : 0.0;
    for (double& v : g.values) v -= shift;
    return g;
}

}  // namespace dmw
