#include "core/undercut.h"

#include <algorithm>
#include <cmath>
#include <limits>
#ifndef __EMSCRIPTEN__
#include <thread>
#endif

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


// Per region face, everything an evaluation needs (computed once per query, reused for every direction).
struct FaceSample {
    Vec3 origin;  // centroid, lifted a hair along the normal (see undercut_map)
    Vec3 normal;  // unit
    double area;
    std::uint32_t face;
};

std::vector<FaceSample> collect_faces(const TriMesh& m, std::span<const std::uint8_t> region) {
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (const Vec3& p : m.positions)
        lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)}, hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
    const double lift = 1e-7 * std::max(norm(hi - lo), 1e-12);
    std::vector<FaceSample> out;
    for (std::size_t f = 0; f < m.triangles.size(); ++f) {
        if (!region.empty() && !region[f]) continue;
        const Vec3 &a = m.positions[m.triangles[f][0]], &b = m.positions[m.triangles[f][1]], &c = m.positions[m.triangles[f][2]];
        const Vec3 n2 = cross(b - a, c - a);
        const double area = 0.5 * norm(n2);
        if (area <= 0.0) continue;
        const Vec3 n = (1.0 / (2.0 * area)) * n2;
        // Start a hair off the surface along the face normal: a face parallel to d (n . d = 0, e.g. a vertical wall)
        // would otherwise cast its ray inside its own plane and graze its neighbours. Real overhangs still hit.
        out.push_back({(1.0 / 3.0) * (a + b + c) + lift * n, n, area, static_cast<std::uint32_t>(f)});
    }
    return out;
}

int thread_count(int threads, std::size_t work) {
#ifdef __EMSCRIPTEN__
    (void)threads, (void)work;
    return 1;
#else
    if (work < 4096) return 1;  // not worth the thread start-up
    const int hw = static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
    return threads <= 0 ? hw : std::min(threads, hw);
#endif
}

// Undercut area of `faces` along unit d; optionally flags undercut faces (by face id) in `flags`.
double evaluate(const std::vector<FaceSample>& faces, const Bvh& bvh, const Vec3& d, double tolerance, int threads,
                std::vector<std::uint8_t>* flags) {
    const int t = thread_count(threads, faces.size());
    auto run = [&](std::size_t begin, std::size_t end) {
        double area = 0.0;
        for (std::size_t i = begin; i < end; ++i) {
            const FaceSample& s = faces[i];
            const bool under = dot(s.normal, d) < -tolerance ||
                               bvh.occluded(s.origin, d, 1e-12, std::numeric_limits<double>::infinity(), s.face);
            if (under) {
                area += s.area;
                if (flags) (*flags)[s.face] = 1;  // distinct faces per thread: no data race
            }
        }
        return area;
    };
    if (t == 1) return run(0, faces.size());
#ifndef __EMSCRIPTEN__
    std::vector<double> part(static_cast<std::size_t>(t), 0.0);
    std::vector<std::thread> pool;
    const std::size_t chunk = (faces.size() + static_cast<std::size_t>(t) - 1) / static_cast<std::size_t>(t);
    for (int k = 0; k < t; ++k) {
        const std::size_t b = static_cast<std::size_t>(k) * chunk, e = std::min(faces.size(), b + chunk);
        pool.emplace_back([&, k, b, e] { part[static_cast<std::size_t>(k)] = b < e ? run(b, e) : 0.0; });
    }
    for (auto& th : pool) th.join();
    double sum = 0.0;
    for (double x : part) sum += x;  // fixed order: deterministic
    return sum;
#else
    return run(0, faces.size());
#endif
}

double total_area(const std::vector<FaceSample>& faces) {
    double a = 0.0;
    for (const FaceSample& s : faces) a += s.area;
    return a;
}

}  // namespace

UndercutResult undercut_map(const TriMesh& m, const Bvh& bvh, const Vec3& dir, std::span<const std::uint8_t> region,
                            double tolerance, int threads) {
    const std::vector<FaceSample> faces = collect_faces(m, region);
    UndercutResult r;
    r.undercut.assign(m.triangles.size(), 0);
    r.region_area = total_area(faces);
    r.undercut_area = evaluate(faces, bvh, unit(dir), tolerance, threads, &r.undercut);
    return r;
}

InsertionAxis best_insertion_axis(const TriMesh& m, const Bvh& bvh, std::span<const std::uint8_t> region, const Vec3& hint_in,
                                  double max_tilt_degrees, int samples, double min_step_degrees, double tolerance,
                                  std::size_t coarse_faces, int threads) {
    const double deg = std::acos(-1.0) / 180.0, max_tilt = max_tilt_degrees * deg;
    const Vec3 hint = unit(hint_in);
    const std::vector<FaceSample> all = collect_faces(m, region);
    std::vector<FaceSample> coarse;
    if (coarse_faces > 0 && all.size() > coarse_faces) {
        const std::size_t stride = (all.size() + coarse_faces - 1) / coarse_faces;
        for (std::size_t i = 0; i < all.size(); i += stride) coarse.push_back(all[i]);
    }
    InsertionAxis best;
    best.axis = hint;
    double best_area = 0.0, best_tilt = 0.0;
    const std::vector<FaceSample>* current = coarse.empty() ? &all : &coarse;
    bool any = false;
    auto consider = [&](const Vec3& dir) {
        const Vec3 d = unit(dir);
        const double tilt = angle_between(d, hint);
        if (tilt > max_tilt + 1e-12) return false;
        const double area = evaluate(*current, bvh, d, tolerance, threads, nullptr);
        ++best.evaluations;
        // Lexicographic: less undercut area, then closer to the hint (a surveyor prefers the familiar path).
        const double scale = std::max(total_area(*current), 1e-300) * 1e-12;
        const bool better = !any || area < best_area - scale || (std::abs(area - best_area) <= scale && tilt < best_tilt);
        if (better) best.axis = d, best_area = area, best_tilt = tilt, any = true;
        return better;
    };
    auto pattern_search = [&](double step, double min_step) {
        while (step >= min_step) {
            Vec3 a, b;
            frame(best.axis, a, b);
            const Vec3 c = best.axis;
            bool moved = false;
            for (const Vec3& k : {a, b})
                for (double sgn : {1.0, -1.0}) moved = consider(rotate(c, k, sgn * step)) || moved;
            if (!moved) step *= 0.5;
        }
    };
    // Global: Fibonacci spiral over the spherical cap of half-angle max_tilt (area-uniform), plus the hint itself.
    Vec3 u, w;
    frame(hint, u, w);
    consider(hint);
    const double golden = std::acos(-1.0) * (3.0 - std::sqrt(5.0));
    for (int i = 0; i < samples; ++i) {
        const double cz = 1.0 - (1.0 - std::cos(max_tilt)) * (i + 0.5) / samples, sz = std::sqrt(std::max(0.0, 1.0 - cz * cz));
        const double phi = golden * i;
        consider(cz * hint + (sz * std::cos(phi)) * u + (sz * std::sin(phi)) * w);
    }
    // Local: pattern search around the best direction, rotating about two axes perpendicular to it.
    pattern_search(std::max(2.0 * max_tilt / std::sqrt(std::max(samples, 1)), min_step_degrees * deg), min_step_degrees * deg);
    if (!coarse.empty()) {
        // Fine: re-score the incumbent on every face, then a short full-resolution search around it.
        current = &all;
        any = false;
        consider(best.axis);
        pattern_search(std::max(1.0 * deg, min_step_degrees * deg), min_step_degrees * deg);
    }
    best.result = undercut_map(m, bvh, best.axis, region, tolerance, threads);
    return best;
}

}  // namespace dmw
