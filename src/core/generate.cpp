#include "core/generate.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <set>
#include <unordered_map>

namespace dmw {
namespace {

constexpr double kTwoPi = 2.0 * std::numbers::pi;

Vec3 scaled_to(const Vec3& p, double radius) {
    const double s = radius / std::sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
    return {p.x * s, p.y * s, p.z * s};
}

}  // namespace

TriMesh make_grid(std::uint32_t nx, std::uint32_t ny) {
    return make_grid_with_holes(nx, ny, {});
}

TriMesh make_grid_with_holes(std::uint32_t nx, std::uint32_t ny,
                             const std::vector<std::array<std::uint32_t, 2>>& hole_cells) {
    if (nx == 0 || ny == 0) return {};
    std::set<std::array<std::uint32_t, 2>> holes;
    for (const auto& c : hole_cells) {
        if (c[0] >= nx || c[1] >= ny) return {};
        holes.insert(c);
    }
    TriMesh m;
    const std::uint32_t row = nx + 1;  // vertex (i, j) has index j * row + i
    for (std::uint32_t j = 0; j <= ny; ++j) {
        for (std::uint32_t i = 0; i <= nx; ++i) {
            m.positions.push_back({double(i) / nx, double(j) / ny, 0.0});
        }
    }
    for (std::uint32_t j = 0; j < ny; ++j) {
        for (std::uint32_t i = 0; i < nx; ++i) {
            if (holes.count({i, j}) != 0) continue;
            // Cell corners counter-clockwise from +z: a (i,j), b (i+1,j), c (i+1,j+1), d (i,j+1).
            const std::uint32_t a = j * row + i, b = a + 1, c = a + row + 1, d = a + row;
            m.triangles.push_back({a, b, c});
            m.triangles.push_back({a, c, d});
        }
    }
    return m;
}

TriMesh make_torus(std::uint32_t nu, std::uint32_t nv, double major_radius, double minor_radius) {
    if (nu < 3 || nv < 3) return {};
    TriMesh m;
    // p(u, v) = ((R + r cos v) cos u, (R + r cos v) sin u, r sin v). The partial derivatives
    // satisfy p_u x p_v = outward normal, so a cell traversed +u then +v is outward-wound.
    for (std::uint32_t j = 0; j < nv; ++j) {
        const double v = kTwoPi * j / nv;
        for (std::uint32_t i = 0; i < nu; ++i) {
            const double u = kTwoPi * i / nu;
            const double ring = major_radius + minor_radius * std::cos(v);
            m.positions.push_back({ring * std::cos(u), ring * std::sin(u), minor_radius * std::sin(v)});
        }
    }
    auto idx = [nu, nv](std::uint32_t i, std::uint32_t j) { return (j % nv) * nu + (i % nu); };
    for (std::uint32_t j = 0; j < nv; ++j) {
        for (std::uint32_t i = 0; i < nu; ++i) {
            const std::uint32_t a = idx(i, j), b = idx(i + 1, j), c = idx(i + 1, j + 1), d = idx(i, j + 1);
            m.triangles.push_back({a, b, c});
            m.triangles.push_back({a, c, d});
        }
    }
    return m;
}

TriMesh make_cylinder(std::uint32_t nu, std::uint32_t nv, double radius, double height) {
    if (nu < 3 || nv < 1) return {};
    TriMesh m;
    // p(u, z) = (r cos u, r sin u, z); p_u x p_z = (r cos u, r sin u, 0) points outward.
    for (std::uint32_t j = 0; j <= nv; ++j) {
        for (std::uint32_t i = 0; i < nu; ++i) {
            const double u = kTwoPi * i / nu;
            m.positions.push_back({radius * std::cos(u), radius * std::sin(u), height * j / nv});
        }
    }
    auto idx = [nu](std::uint32_t i, std::uint32_t j) { return j * nu + (i % nu); };
    for (std::uint32_t j = 0; j < nv; ++j) {
        for (std::uint32_t i = 0; i < nu; ++i) {
            const std::uint32_t a = idx(i, j), b = idx(i + 1, j), c = idx(i + 1, j + 1), d = idx(i, j + 1);
            m.triangles.push_back({a, b, c});
            m.triangles.push_back({a, c, d});
        }
    }
    return m;
}

TriMesh make_icosphere(std::uint32_t subdivisions, double radius) {
    if (subdivisions > 10) return {};  // 10 * 4^10 + 2 ~ 10.5M vertices; beyond is a mistake
    // Regular icosahedron: vertices are cyclic permutations of (0, +-1, +-phi).
    const double t = std::numbers::phi;
    TriMesh m;
    for (const Vec3& p : {Vec3{-1, t, 0}, Vec3{1, t, 0}, Vec3{-1, -t, 0}, Vec3{1, -t, 0},
                          Vec3{0, -1, t}, Vec3{0, 1, t}, Vec3{0, -1, -t}, Vec3{0, 1, -t},
                          Vec3{t, 0, -1}, Vec3{t, 0, 1}, Vec3{-t, 0, -1}, Vec3{-t, 0, 1}}) {
        m.positions.push_back(scaled_to(p, radius));
    }
    m.triangles = {{0, 11, 5}, {0, 5, 1},  {0, 1, 7},   {0, 7, 10}, {0, 10, 11},
                   {1, 5, 9},  {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
                   {3, 9, 4},  {3, 4, 2},  {3, 2, 6},   {3, 6, 8},  {3, 8, 9},
                   {4, 9, 5},  {2, 4, 11}, {6, 2, 10},  {8, 6, 7},  {9, 8, 1}};

    for (std::uint32_t s = 0; s < subdivisions; ++s) {
        // One new vertex per undirected edge, shared by both adjacent triangles. Keyed on
        // (min, max) so the two triangles, which see the edge in opposite directions, agree.
        std::unordered_map<std::uint64_t, std::uint32_t> midpoint;
        auto mid = [&](std::uint32_t a, std::uint32_t b) {
            const std::uint64_t key = (std::uint64_t{std::min(a, b)} << 32) | std::max(a, b);
            const auto [it, inserted] =
                midpoint.try_emplace(key, static_cast<std::uint32_t>(m.positions.size()));
            if (inserted) {
                const Vec3 &p = m.positions[a], &q = m.positions[b];
                m.positions.push_back(scaled_to({p.x + q.x, p.y + q.y, p.z + q.z}, radius));
            }
            return it->second;
        };
        std::vector<std::array<std::uint32_t, 3>> next;
        next.reserve(4 * m.triangles.size());
        for (const auto& [a, b, c] : m.triangles) {
            const std::uint32_t ab = mid(a, b), bc = mid(b, c), ca = mid(c, a);
            // Corner triangles keep the parent's winding; the center one (ab, bc, ca) too.
            next.push_back({a, ab, ca});
            next.push_back({b, bc, ab});
            next.push_back({c, ca, bc});
            next.push_back({ab, bc, ca});
        }
        m.triangles = std::move(next);
    }
    return m;
}

TriMesh make_mobius(std::uint32_t segments, double radius, double half_width) {
    if (segments < 3) return {};
    TriMesh m;
    // Center circle c(theta); the cross-section direction d rotates by theta/2, so after a
    // full turn d(2 pi) = -d(0): the strip comes back flipped (the half twist).
    for (std::uint32_t i = 0; i < segments; ++i) {
        const double th = kTwoPi * i / segments;
        const double cx = radius * std::cos(th), cy = radius * std::sin(th);
        const double dr = std::cos(th / 2), dz = std::sin(th / 2);  // d = dr * radial + dz * z
        const double ox = half_width * dr * std::cos(th), oy = half_width * dr * std::sin(th),
                     oz = half_width * dz;
        m.positions.push_back({cx + ox, cy + oy, oz});  // top_i    = 2i
        m.positions.push_back({cx - ox, cy - oy, -oz});  // bottom_i = 2i + 1
    }
    for (std::uint32_t i = 0; i < segments; ++i) {
        const std::uint32_t t = 2 * i, b = 2 * i + 1;
        // At the seam, top_n coincides with bottom_0 and bottom_n with top_0.
        const bool seam = (i + 1 == segments);
        const std::uint32_t tn = seam ? 1 : 2 * (i + 1), bn = seam ? 0 : 2 * (i + 1) + 1;
        m.triangles.push_back({t, b, bn});
        m.triangles.push_back({t, bn, tn});
    }
    return m;
}

void append(TriMesh& dst, const TriMesh& src) {
    const auto offset = static_cast<std::uint32_t>(dst.positions.size());
    dst.positions.insert(dst.positions.end(), src.positions.begin(), src.positions.end());
    for (const auto& [a, b, c] : src.triangles) dst.triangles.push_back({a + offset, b + offset, c + offset});
}

}  // namespace dmw
