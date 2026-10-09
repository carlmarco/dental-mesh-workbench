#include "core/isosurface.h"

#include <array>
#include <cmath>
#include <unordered_map>

#include "detail/vec.h"

namespace dmw {
namespace {

using namespace detail;

// The 6 tetrahedra of the Kuhn subdivision of a cube along its diagonal from corner 0 to corner 7 (corner bits:
// 1 = +x, 2 = +y, 4 = +z): one per ordering of the three axes, walking 0 -> e_a -> e_a + e_b -> 7.
constexpr std::array<std::array<int, 4>, 6> kTets{{{0, 1, 3, 7}, {0, 1, 5, 7}, {0, 2, 3, 7}, {0, 2, 6, 7}, {0, 4, 5, 7}, {0, 4, 6, 7}}};

}  // namespace

TriMesh extract_isosurface(const Grid3& g, double iso) {
    TriMesh out;
    if (g.n[0] < 2 || g.n[1] < 2 || g.n[2] < 2) return out;
    std::unordered_map<std::uint64_t, std::uint32_t> edge_vertex;  // key: (min node, max node)
    auto vertex_on = [&](std::size_t a, std::size_t b, const Vec3& pa, const Vec3& pb, double sa, double sb) {
        // A crossing exactly at a node (the field equals iso there; nodes with s == 0 count as positive) is one shared
        // vertex keyed by the node, not one per incident edge: otherwise coincident vertices make zero-area triangles
        // whose orientation is undefined (found on a sphere field whose radius passes through grid nodes).
        if (sb == 0.0 || sa == 0.0) {
            const std::size_t at = sb == 0.0 ? b : a;
            const std::uint64_t key = static_cast<std::uint64_t>(g.values.size()) * static_cast<std::uint64_t>(g.values.size()) + at;
            const auto it = edge_vertex.find(key);
            if (it != edge_vertex.end()) return it->second;
            const auto id = static_cast<std::uint32_t>(out.positions.size());
            out.positions.push_back(sb == 0.0 ? pb : pa);
            edge_vertex.emplace(key, id);
            return id;
        }
        const std::uint64_t lo = std::min(a, b), hi = std::max(a, b);
        const std::uint64_t key = lo * static_cast<std::uint64_t>(g.values.size()) + hi;
        const auto it = edge_vertex.find(key);
        if (it != edge_vertex.end()) return it->second;
        // Interpolate from the lower node id so the shared vertex is bit-identical from either tetrahedron.
        const bool a_first = a < b;
        const Vec3 &p0 = a_first ? pa : pb, &p1 = a_first ? pb : pa;
        const double s0 = a_first ? sa : sb, s1 = a_first ? sb : sa;
        const double t = s0 / (s0 - s1);
        const auto id = static_cast<std::uint32_t>(out.positions.size());
        out.positions.push_back(p0 + t * (p1 - p0));
        edge_vertex.emplace(key, id);
        return id;
    };
    for (std::size_t k = 0; k + 1 < g.n[2]; ++k)
        for (std::size_t j = 0; j + 1 < g.n[1]; ++j)
            for (std::size_t i = 0; i + 1 < g.n[0]; ++i) {
                std::array<std::size_t, 8> node;
                std::array<Vec3, 8> pos;
                std::array<double, 8> s;
                bool any_neg = false, any_pos = false;
                for (int c = 0; c < 8; ++c) {
                    const std::size_t ci = i + (c & 1), cj = j + ((c >> 1) & 1), ck = k + ((c >> 2) & 1);
                    node[static_cast<std::size_t>(c)] = g.index(ci, cj, ck);
                    pos[static_cast<std::size_t>(c)] = g.position(ci, cj, ck);
                    // Values within 1e-10 h of the iso value snap to exactly 0, so every near-node crossing takes the shared
                    // node vertex below: otherwise a crossing that rounds onto a node (t = 1 in floating point) duplicates
                    // it, and whether that happens depends on the platform's rounding (x86 vs ARM fused multiply-add).
                    double sv = g.values[node[static_cast<std::size_t>(c)]] - iso;
                    if (std::abs(sv) < 1e-10 * g.h) sv = 0.0;
                    s[static_cast<std::size_t>(c)] = sv;
                    (s[static_cast<std::size_t>(c)] < 0.0 ? any_neg : any_pos) = true;
                }
                if (!any_neg || !any_pos) continue;
                for (const auto& tet : kTets) {
                    std::array<std::size_t, 4> neg{}, posv{};  // local corner ids (0..3) by sign
                    std::size_t nn = 0, np = 0;
                    for (std::size_t q = 0; q < 4; ++q) (s[static_cast<std::size_t>(tet[q])] < 0.0 ? neg[nn++] : posv[np++]) = q;
                    if (nn == 0 || np == 0) continue;
                    auto corner = [&](std::size_t q) { return static_cast<std::size_t>(tet[q]); };
                    auto cut = [&](std::size_t qa, std::size_t qb) {
                        const std::size_t a = corner(qa), b = corner(qb);
                        return vertex_on(node[a], node[b], pos[a], pos[b], s[a], s[b]);
                    };
                    // Orientation reference: the gradient of the linear interpolant over this tetrahedron (exact, unlike a
                    // centroid difference, which can point the wrong way in a slanted tetrahedron). With e_q = p_q - p_0
                    // and d_q = s_q - s_0, grad = (d1 (e2 x e3) + d2 (e3 x e1) + d3 (e1 x e2)) / (e1 . (e2 x e3)); only
                    // its direction matters, so the sign of the denominator is applied and its size dropped.
                    const Vec3 e1 = pos[corner(1)] - pos[corner(0)], e2 = pos[corner(2)] - pos[corner(0)], e3 = pos[corner(3)] - pos[corner(0)];
                    const double d1 = s[corner(1)] - s[corner(0)], d2 = s[corner(2)] - s[corner(0)], d3 = s[corner(3)] - s[corner(0)];
                    const double vol = dot(e1, cross(e2, e3));
                    const Vec3 up = (vol < 0.0 ? -1.0 : 1.0) * (d1 * cross(e2, e3) + d2 * cross(e3, e1) + d3 * cross(e1, e2));
                    auto emit = [&](std::uint32_t a, std::uint32_t b, std::uint32_t c) {
                        if (a == b || b == c || a == c) return;
                        const Vec3 n = cross(out.positions[b] - out.positions[a], out.positions[c] - out.positions[a]);
                        if (dot(n, up) < 0.0) std::swap(b, c);
                        out.triangles.push_back({a, b, c});
                    };
                    if (nn == 1 || np == 1) {  // one triangle around the lone vertex
                        const bool lone_neg = nn == 1;
                        const std::size_t lone = lone_neg ? neg[0] : posv[0];
                        const auto& others = lone_neg ? posv : neg;
                        emit(cut(lone, others[0]), cut(lone, others[1]), cut(lone, others[2]));
                    } else {  // two and two: a quad, split into two triangles
                        const std::uint32_t v00 = cut(neg[0], posv[0]), v01 = cut(neg[0], posv[1]), v11 = cut(neg[1], posv[1]),
                                            v10 = cut(neg[1], posv[0]);
                        emit(v00, v01, v11);
                        emit(v00, v11, v10);
                    }
                }
            }
    return out;
}

}  // namespace dmw
