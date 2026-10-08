#include "core/undercut.h"

#include <algorithm>
#include <cmath>

#include "detail/vec.h"

namespace dmw {
namespace {

using namespace detail;

Vec3 unit(const Vec3& v) {
    const double n = norm(v);
    return n > 0.0 ? (1.0 / n) * v : Vec3{0, 0, 1};
}

// Rotate `v` about unit axis `k` by `angle` (Rodrigues).
Vec3 rotate(const Vec3& v, const Vec3& k, double angle) {
    const double c = std::cos(angle), s = std::sin(angle);
    return c * v + s * cross(k, v) + (dot(k, v) * (1.0 - c)) * k;
}

// Two unit vectors completing `a` to an orthonormal frame.
void frame(const Vec3& a, Vec3& u, Vec3& w) {
    const Vec3 t = std::abs(a.x) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
    u = unit(cross(a, t));
    w = cross(a, u);
}

double angle_between(const Vec3& a, const Vec3& b) { return std::acos(std::clamp(dot(a, b), -1.0, 1.0)); }

}  // namespace

UndercutResult undercut_map(const TriMesh& m, const Bvh& bvh, const Vec3& dir, std::span<const std::uint8_t> region,
                            double tolerance) {
    const Vec3 d = unit(dir);
    const std::size_t nf = m.triangles.size();
    UndercutResult r;
    r.undercut.assign(nf, 0);
    // Scale for the self-intersection offsets: rays start a hair off the surface and skip their own face.
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (const Vec3& p : m.positions)
        lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)}, hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
    const double eps = 1e-9 * std::max(norm(hi - lo), 1e-12);
    for (std::size_t f = 0; f < nf; ++f) {
        if (!region.empty() && !region[f]) continue;
        const Vec3 &a = m.positions[m.triangles[f][0]], &b = m.positions[m.triangles[f][1]], &c = m.positions[m.triangles[f][2]];
        const Vec3 n2 = cross(b - a, c - a);  // |n2| = 2 * area
        const double area = 0.5 * norm(n2);
        if (area <= 0.0) continue;
        r.region_area += area;
        bool under = dot(n2, d) / (2.0 * area) < -tolerance;  // back-facing
        if (!under) {
            // Start a hair off the surface along the face normal: a face parallel to d (n . d = 0, e.g. a vertical wall)
            // would otherwise cast its ray inside its own plane and graze its neighbours. Real overhangs still hit.
            const Vec3 origin = (1.0 / 3.0) * (a + b + c) + (100.0 * eps / (2.0 * area)) * n2;
            under = bvh.occluded(origin, d, eps, std::numeric_limits<double>::infinity(), static_cast<std::uint32_t>(f));
        }
        if (under) r.undercut[f] = 1, r.undercut_area += area;
    }
    return r;
}

InsertionAxis best_insertion_axis(const TriMesh& m, const Bvh& bvh, std::span<const std::uint8_t> region, const Vec3& hint_in,
                                  double max_tilt_degrees, int samples, double min_step_degrees, double tolerance) {
    const double deg = std::acos(-1.0) / 180.0, max_tilt = max_tilt_degrees * deg;
    const Vec3 hint = unit(hint_in);
    Vec3 u, w;
    frame(hint, u, w);
    InsertionAxis best;
    double best_tilt = 0.0;
    auto consider = [&](const Vec3& dir) {
        const Vec3 d = unit(dir);
        const double tilt = angle_between(d, hint);
        if (tilt > max_tilt + 1e-12) return false;
        UndercutResult r = undercut_map(m, bvh, d, region, tolerance);
        ++best.evaluations;
        // Lexicographic: less undercut area, then closer to the hint (a surveyor prefers the familiar path).
        const double scale = std::max(r.region_area, 1e-300) * 1e-12;
        const bool better = best.evaluations == 1 || r.undercut_area < best.result.undercut_area - scale ||
                            (std::abs(r.undercut_area - best.result.undercut_area) <= scale && tilt < best_tilt);
        if (better) best.axis = d, best.result = std::move(r), best_tilt = tilt;
        return better;
    };
    // Global: Fibonacci spiral over the spherical cap of half-angle max_tilt (area-uniform), plus the hint itself.
    consider(hint);
    const double golden = std::acos(-1.0) * (3.0 - std::sqrt(5.0));
    for (int i = 0; i < samples; ++i) {
        const double cz = 1.0 - (1.0 - std::cos(max_tilt)) * (i + 0.5) / samples, sz = std::sqrt(std::max(0.0, 1.0 - cz * cz));
        const double phi = golden * i;
        consider(cz * hint + (sz * std::cos(phi)) * u + (sz * std::sin(phi)) * w);
    }
    // Local: pattern search around the best direction, rotating about two axes perpendicular to it.
    double step = std::max(2.0 * max_tilt / std::sqrt(std::max(samples, 1)), min_step_degrees * deg);
    while (step >= min_step_degrees * deg) {
        Vec3 a, b;
        frame(best.axis, a, b);
        const Vec3 c = best.axis;
        bool moved = false;
        for (const Vec3& k : {a, b})
            for (double sgn : {1.0, -1.0}) moved = consider(rotate(c, k, sgn * step)) || moved;
        if (!moved) step *= 0.5;
    }
    return best;
}

}  // namespace dmw
