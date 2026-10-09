#include "core/thickness.h"

#include <algorithm>
#include <cmath>
#ifndef __EMSCRIPTEN__
#include <thread>
#endif

#include "detail/vec.h"

namespace dmw {
namespace {

using namespace detail;

Vec3 unit(const Vec3& v) {
    const double n = norm(v);
    return n > 0.0 ? (1.0 / n) * v : Vec3{0, 0, 0};
}

}  // namespace

ThicknessField wall_thickness(const TriMesh& m, const Bvh& bvh, const ThicknessParams& prm) {
    const std::size_t nv = m.positions.size();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    ThicknessField out;
    out.along_normal.assign(nv, nan), out.cone_min.assign(nv, nan), out.cone_median.assign(nv, nan);
    // Area-weighted vertex normals and unit face normals (outward orientation assumed).
    std::vector<Vec3> vn(nv, Vec3{}), fn(m.triangles.size());
    for (std::size_t f = 0; f < m.triangles.size(); ++f) {
        const auto& t = m.triangles[f];
        const Vec3 n2 = cross(m.positions[t[1]] - m.positions[t[0]], m.positions[t[2]] - m.positions[t[0]]);
        fn[f] = unit(n2);
        for (auto v : t) vn[v] += n2;
    }
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (const Vec3& p : m.positions)
        lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)}, hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
    const double t_min = 1e-6 * std::max(norm(hi - lo), 1e-12);  // skips the faces around the ray's own vertex
    // Unit offsets of the cone rays around the local z axis (Fibonacci spiral over the cap, area-uniform).
    const double cone = prm.cone_degrees * std::acos(-1.0) / 180.0, golden = std::acos(-1.0) * (3.0 - std::sqrt(5.0));
    std::vector<Vec3> local;
    for (int i = 0; i < prm.rays; ++i) {
        const double cz = 1.0 - (1.0 - std::cos(cone)) * (i + 0.5) / prm.rays, sz = std::sqrt(std::max(0.0, 1.0 - cz * cz));
        local.push_back({sz * std::cos(golden * i), sz * std::sin(golden * i), cz});
    }
    // Length of a ray into the solid if it leaves through the opposite wall (hit face pointing along the ray).
    auto shoot = [&](const Vec3& o, const Vec3& d) {
        const RayHit h = bvh.intersect(o, d, t_min, prm.max_distance);
        return h.hit && dot(fn[h.face], d) > 0.0 ? h.t : nan;
    };
    auto run = [&](std::size_t begin, std::size_t end) {
        std::vector<double> accepted;
        for (std::size_t v = begin; v < end; ++v) {
            const Vec3 d0 = -1.0 * unit(vn[v]);
            if (norm(d0) == 0.0) continue;
            const Vec3& o = m.positions[v];
            out.along_normal[v] = shoot(o, d0);
            // Frame around d0 for the cone.
            const Vec3 t = std::abs(d0.x) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
            const Vec3 u = unit(cross(d0, t)), w = cross(d0, u);
            accepted.clear();
            if (std::isfinite(out.along_normal[v])) accepted.push_back(out.along_normal[v]);
            for (const Vec3& l : local) {
                const double len = shoot(o, l.x * u + l.y * w + l.z * d0);
                if (std::isfinite(len)) accepted.push_back(len);
            }
            if (accepted.empty()) continue;
            std::sort(accepted.begin(), accepted.end());
            out.cone_min[v] = accepted.front();
            out.cone_median[v] = accepted[accepted.size() / 2];
        }
    };
#ifdef __EMSCRIPTEN__
    run(0, nv);
#else
    const int hw = static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
    const int threads = nv < 2048 ? 1 : (prm.threads <= 0 ? hw : std::min(prm.threads, hw));
    if (threads == 1) {
        run(0, nv);
    } else {
        std::vector<std::thread> pool;
        const std::size_t chunk = (nv + static_cast<std::size_t>(threads) - 1) / static_cast<std::size_t>(threads);
        for (int k = 0; k < threads; ++k) {
            const std::size_t b = static_cast<std::size_t>(k) * chunk, e = std::min(nv, b + chunk);
            if (b < e) pool.emplace_back(run, b, e);  // distinct vertices per thread: no data race
        }
        for (auto& th : pool) th.join();
    }
#endif
    return out;
}

}  // namespace dmw
