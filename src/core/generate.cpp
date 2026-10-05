#include "core/generate.h"

#include <algorithm>
#include <array>
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

TriMesh make_plate_with_handle(std::uint32_t n, std::uint32_t tube_segments, double arch_height) {
    if (n < 8 || tube_segments < 2) return {};
    const std::uint32_t j = n / 2, ia = n / 4, ib = 3 * n / 4 - 1;
    TriMesh m = make_grid_with_holes(n, n, {{ia, j}, {ib, j}});
    const std::uint32_t row = n + 1;
    auto corner = [row, j](std::uint32_t i, int dx, int dy) {
        return (j + static_cast<std::uint32_t>(dy)) * row + i + static_cast<std::uint32_t>(dx);
    };
    // Hole corners counter-clockwise seen from +z: (-,-), (+,-), (+,+), (-,+).
    const std::array<std::uint32_t, 4> a{corner(ia, 0, 0), corner(ia, 1, 0), corner(ia, 1, 1), corner(ia, 0, 1)};
    const std::array<std::uint32_t, 4> b{corner(ib, 0, 0), corner(ib, 1, 0), corner(ib, 1, 1), corner(ib, 0, 1)};

    // The plate traverses each hole rim clockwise (seen from +z), so the tube must traverse
    // both rims counter-clockwise. Ring t sits on an arch from hole A's center to hole B's;
    // its corner offsets rotate by pi * t about the y axis, so at hole B the ring arrives
    // upside down: corner k of A lands on B's corner order [b1, b0, b3, b2]. That mirrored
    // order is exactly what keeps the winding consistent (checked by the topology tests).
    const double h = 1.0 / n;  // cell size
    const double cax = (ia + 0.5) * h, cbx = (ib + 0.5) * h, cy = (j + 0.5) * h;
    const double ox[4] = {-h / 2, h / 2, h / 2, -h / 2}, oy[4] = {-h / 2, -h / 2, h / 2, h / 2};
    std::vector<std::array<std::uint32_t, 4>> rings{a};
    for (std::uint32_t t = 1; t < tube_segments; ++t) {
        const double s = static_cast<double>(t) / tube_segments, phi = std::numbers::pi * s;
        const double cx = cax + (cbx - cax) * (1 - std::cos(phi)) / 2, cz = arch_height * std::sin(phi);
        std::array<std::uint32_t, 4> ring{};
        for (std::size_t k = 0; k < 4; ++k) {
            ring[k] = static_cast<std::uint32_t>(m.positions.size());
            m.positions.push_back({cx + ox[k] * std::cos(phi), cy + oy[k], cz - ox[k] * std::sin(phi)});
        }
        rings.push_back(ring);
    }
    rings.push_back({b[1], b[0], b[3], b[2]});
    for (std::size_t t = 0; t + 1 < rings.size(); ++t) {
        for (std::size_t k = 0; k < 4; ++k) {
            const std::uint32_t p0 = rings[t][k], p1 = rings[t][(k + 1) % 4];
            const std::uint32_t q0 = rings[t + 1][k], q1 = rings[t + 1][(k + 1) % 4];
            m.triangles.push_back({p0, p1, q1});  // contains p0 -> p1: counter-clockwise on rim A
            m.triangles.push_back({p0, q1, q0});
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

TriMesh make_brick_grid(std::uint32_t n, double shift) {
    TriMesh m = make_grid(n, n);
    for (auto& p : m.positions) {
        const long j = std::lround(p.y * n);
        if (p.x > 0 && p.x < 1) p.x += (j % 2 != 0 ? 1.0 : -1.0) * shift / n;
    }
    return m;
}

TriMesh make_defect_showcase() {
    TriMesh m = make_grid(4, 4);                  // vertices 0..24, faces 0..31
    m.triangles[2 * (2 * 4 + 2)] = {12, 18, 13};  // cell (2,2), first triangle, flipped
    const auto apex = static_cast<std::uint32_t>(m.positions.size());
    m.positions.push_back({0.375, 0.25, 0.3});
    m.triangles.push_back({6, 7, apex});  // fin: edge 6-7 now has 3 faces
    const auto p = static_cast<std::uint32_t>(m.positions.size());
    m.positions.push_back({1.25, 1.0, 0.0});
    m.positions.push_back({1.25, 1.25, 0.0});
    m.triangles.push_back({24, p, p + 1});         // touches the grid only at vertex 24
    m.triangles.push_back({16, 21, 15});           // same vertex set as face 24 (15, 16, 21)
    m.triangles.push_back({0, 0, 1});              // repeated index
    m.positions.push_back({-0.5, -0.5, 0.0});      // referenced by nothing
    TriMesh strip = make_mobius(12, 0.35, 0.12);
    for (auto& q : strip.positions) q.x += 2.0, q.y += 0.5;
    append(m, strip);
    return m;
}

void append(TriMesh& dst, const TriMesh& src) {
    const auto offset = static_cast<std::uint32_t>(dst.positions.size());
    dst.positions.insert(dst.positions.end(), src.positions.begin(), src.positions.end());
    for (const auto& [a, b, c] : src.triangles) dst.triangles.push_back({a + offset, b + offset, c + offset});
}

}  // namespace dmw
