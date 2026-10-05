#include "core/intrinsic.h"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

#include "detail/vec.h"

namespace dmw {
namespace {

// Side lengths of face f: l[k] = |v_k v_{k+1}| (the length of half-edge 3f + k).
std::array<double, 3> sides(const IntrinsicTriangulation& t, std::size_t f) {
    return {t.length[3 * f], t.length[3 * f + 1], t.length[3 * f + 2]};
}

// Triangle area from side lengths, Kahan's numerically stable form of Heron's formula
// (plain Heron loses all precision on slivers). Roundoff-violated triangle inequality -> 0.
double area_from_sides(double a, double b, double c) {
    if (a < b) std::swap(a, b);
    if (b < c) std::swap(b, c);
    if (a < b) std::swap(a, b);  // now a >= b >= c
    const double p = (a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c));
    return p > 0.0 ? 0.25 * std::sqrt(p) : 0.0;
}

// Corner k of a face with sides l: adjacent sides l[k] (to k+1) and l[k+2] (from k+2), opposite
// side l[k+1]. Law of cosines: cos = (a^2 + b^2 - c^2) / 2ab, sin = 2A / ab, so
//   cot = (a^2 + b^2 - c^2) / 4A,  angle = atan2(4A, a^2 + b^2 - c^2).
double corner_numerator(const std::array<double, 3>& l, std::size_t k) {
    const double a = l[k], b = l[(k + 2) % 3], c = l[(k + 1) % 3];
    return a * a + b * b - c * c;
}

// 2D layout of face f: corner 0 at the origin, corner 1 on +x, corner 2 above.
std::array<Vec2, 3> layout(const std::array<double, 3>& l) {
    const double x = (l[0] * l[0] + l[2] * l[2] - l[1] * l[1]) / (2.0 * l[0]);
    return {Vec2{0.0, 0.0}, Vec2{l[0], 0.0}, Vec2{x, std::sqrt(std::max(l[2] * l[2] - x * x, 0.0))}};
}

Vec2 sub(const Vec2& a, const Vec2& b) { return {a[0] - b[0], a[1] - b[1]}; }
double dot2(const Vec2& a, const Vec2& b) { return a[0] * b[0] + a[1] * b[1]; }
double cross2(const Vec2& a, const Vec2& b) { return a[0] * b[1] - a[1] * b[0]; }
double cot2(const Vec2& p, const Vec2& q, const Vec2& r) {  // angle at p
    const Vec2 u = sub(q, p), w = sub(r, p);
    const double s = std::abs(cross2(u, w));
    return s > 0.0 ? dot2(u, w) / s : 0.0;
}

// Cotangent of the angle opposite half-edge h, inside h's face.
double cot_opposite(const IntrinsicTriangulation& t, std::uint32_t h) {
    const auto l = sides(t, face(h));
    const double area = area_from_sides(l[0], l[1], l[2]);
    const std::size_t k = (h + 2) % 3;  // corner opposite side h % 3 is (h % 3) + 2
    return area > 0.0 ? corner_numerator(l, k) / (4.0 * area) : 0.0;
}

bool edge_exists(const HalfEdgeMesh& m, std::uint32_t c, std::uint32_t d) {
    const auto ring = one_ring(m, c);
    return std::find(ring.begin(), ring.end(), d) != ring.end();
}

// Flip the interior edge of half-edge h. Before: f = (a, b, c) with h: a->b, g = (b, a, d)
// with twin: b->a. The quad's boundary, counter-clockwise, is a->d->b->c. After:
//   f = (c, a, d): slots 3f+0 c->a, 3f+1 a->d, 3f+2 d->c
//   g = (d, b, c): slots 3g+0 d->b, 3g+1 b->c, 3g+2 c->d   (new diagonal: 3f+2 / 3g+2)
void flip(IntrinsicTriangulation& t, std::uint32_t h) {
    HalfEdgeMesh& m = t.connectivity;
    std::vector<double>& len = t.length;
    const std::uint32_t tw = m.twin[h];
    const std::uint32_t hn = next(h), hp = prev(h), tn = next(tw), tp = prev(tw);
    const std::uint32_t a = m.origin[h], b = m.origin[tw], c = m.origin[hp], d = m.origin[tp];
    const double l_ab = len[h], l_bc = len[hn], l_ca = len[hp], l_ad = len[tn], l_db = len[tp];

    // New diagonal length: unfold both triangles into the plane, a = (0,0), b = (l_ab, 0);
    // c lies above the x axis (f is counter-clockwise), d below (g is too, traversed b -> a).
    const double xc = (l_ab * l_ab + l_ca * l_ca - l_bc * l_bc) / (2.0 * l_ab);
    const double yc = std::sqrt(std::max(l_ca * l_ca - xc * xc, 0.0));
    const double xd = (l_ab * l_ab + l_ad * l_ad - l_db * l_db) / (2.0 * l_ab);
    const double yd = -std::sqrt(std::max(l_ad * l_ad - xd * xd, 0.0));
    const double l_cd = std::hypot(xc - xd, yc - yd);

    const std::uint32_t outer_twin[4] = {m.twin[hp], m.twin[tn], m.twin[tp], m.twin[hn]};
    const std::uint32_t f0 = 3 * face(h), g0 = 3 * face(tw);
    // Remap stored outgoing half-edges of a, b, c, d whose slot moved (captured before rewrite).
    auto remap = [&](std::uint32_t old) -> std::uint32_t {
        if (old == hp) return f0;                  // c->a
        if (old == tn || old == h) return f0 + 1;  // a->d (a's old a->b is gone)
        if (old == tp) return g0;                  // d->b
        if (old == hn || old == tw) return g0 + 1;  // b->c (b's old b->a is gone)
        return old;
    };
    for (std::uint32_t v : {a, b, c, d}) m.vertex_halfedge[v] = remap(m.vertex_halfedge[v]);

    auto set = [&](std::uint32_t slot, std::uint32_t origin, double l, std::uint32_t twin) {
        m.origin[slot] = origin;
        len[slot] = l;
        m.twin[slot] = twin;
        if (twin != kInvalid) m.twin[twin] = slot;  // redirect the neighbor (or the new diagonal)
    };
    set(f0, c, l_ca, outer_twin[0]);
    set(f0 + 1, a, l_ad, outer_twin[1]);
    set(g0, d, l_db, outer_twin[2]);
    set(g0 + 1, b, l_bc, outer_twin[3]);
    set(f0 + 2, d, l_cd, g0 + 2);
    set(g0 + 2, c, l_cd, f0 + 2);
}

}  // namespace

IntrinsicTriangulation intrinsic_from_mesh(const HalfEdgeMesh& mesh) {
    IntrinsicTriangulation t;
    t.connectivity = mesh;
    t.length.resize(mesh.origin.size());
    for (std::uint32_t h = 0; h < mesh.origin.size(); ++h) {
        t.length[h] = detail::norm(detail::operator-(mesh.positions[dest(mesh, h)], mesh.positions[mesh.origin[h]]));
    }
    return t;
}

bool is_delaunay(const IntrinsicTriangulation& t, std::uint32_t h) {
    const std::uint32_t tw = t.connectivity.twin[h];
    if (tw == kInvalid) return true;
    return cot_opposite(t, h) + cot_opposite(t, tw) >= -1e-12;  // tolerance: no flip-flopping on ties
}

void flip_to_delaunay(IntrinsicTriangulation& t) {
    const HalfEdgeMesh& m = t.connectivity;
    // Work list of half-edge slots. Slots get reused by flips, so an entry may later name a
    // different edge; re-testing it is harmless. Every edge whose status a flip can change
    // (the quad's four sides) is pushed again, so nothing is missed.
    std::vector<std::uint32_t> work;
    for (std::uint32_t h = 0; h < m.origin.size(); ++h) {
        if (m.twin[h] != kInvalid && h < m.twin[h]) work.push_back(h);
    }
    while (!work.empty()) {
        const std::uint32_t h = work.back();
        work.pop_back();
        if (is_delaunay(t, h)) continue;
        const std::uint32_t c = m.origin[prev(h)], d = m.origin[prev(m.twin[h])];
        // D53: the flipped edge c-d must be new and not a loop, or the result would need
        // a multi-edge / self-loop, which vertex-triple connectivity cannot represent.
        if (c == d || edge_exists(m, c, d)) continue;
        const std::uint32_t f0 = 3 * face(h), g0 = 3 * face(m.twin[h]);
        flip(t, h);
        ++t.flips;
        for (std::uint32_t s : {f0, f0 + 1, g0, g0 + 1}) work.push_back(s);
    }
    t.skipped = 0;
    for (std::uint32_t h = 0; h < m.origin.size(); ++h) {
        if (m.twin[h] != kInvalid && h < m.twin[h] && !is_delaunay(t, h)) ++t.skipped;
    }
}

IntrinsicTriangulation intrinsic_delaunay(const HalfEdgeMesh& mesh) {
    IntrinsicTriangulation t = intrinsic_from_mesh(mesh);
    flip_to_delaunay(t);
    return t;
}

DecOperators build_dec(const IntrinsicTriangulation& t) {
    const HalfEdgeMesh& m = t.connectivity;
    const auto nv = static_cast<std::uint32_t>(m.positions.size());
    const auto nh = static_cast<std::uint32_t>(m.origin.size());
    const std::uint32_t nf = nh / 3;
    DecOperators ops;

    // Edges from half-edges: one per twin pair (at the smaller id) or boundary half-edge,
    // oriented low -> high vertex index.
    std::vector<std::uint32_t> edge_of(nh, kInvalid);
    double total_length = 0.0;
    for (std::uint32_t h = 0; h < nh; ++h) {
        if (edge_of[h] != kInvalid) continue;
        const std::uint32_t a = m.origin[h], b = dest(m, h);
        const auto e = static_cast<std::uint32_t>(ops.edges.size());
        ops.edges.push_back({std::min(a, b), std::max(a, b)});
        edge_of[h] = e;
        if (m.twin[h] != kInvalid) edge_of[m.twin[h]] = e;
        total_length += t.length[h];
    }
    const auto ne = static_cast<std::uint32_t>(ops.edges.size());
    ops.mean_edge_length = ne > 0 ? total_length / ne : 0.0;

    std::vector<Triplet> d0, d1;
    for (std::uint32_t e = 0; e < ne; ++e) {
        d0.push_back({e, ops.edges[e].v0, -1.0});
        d0.push_back({e, ops.edges[e].v1, +1.0});
    }
    ops.star0.assign(nv, 0.0);
    ops.star1.assign(ne, 0.0);
    for (std::uint32_t f = 0; f < nf; ++f) {
        const auto l = sides(t, f);
        const double area = area_from_sides(l[0], l[1], l[2]);
        for (std::uint32_t k = 0; k < 3; ++k) {
            const std::uint32_t h = 3 * f + k, e = edge_of[h];
            d1.push_back({f, e, m.origin[h] == ops.edges[e].v0 ? 1.0 : -1.0});
            // The angle opposite side k is at corner k+2; a zero-area face contributes nothing.
            if (area > 0.0) ops.star1[e] += 0.5 * corner_numerator(l, (k + 2) % 3) / (4.0 * area);
            ops.star0[m.origin[h]] += area / 3.0;
        }
    }
    ops.d0 = SparseMatrix::from_triplets(ne, nv, std::move(d0));
    ops.d1 = SparseMatrix::from_triplets(nf, ne, std::move(d1));
    ops.laplacian = weighted_gram(ops.d0, ops.star1);  // L = -d0^T *1 d0
    for (double& v : ops.laplacian.value) v = -v;
    return ops;
}

std::vector<double> mixed_area(const IntrinsicTriangulation& t) {
    const HalfEdgeMesh& m = t.connectivity;
    std::vector<double> area(m.positions.size(), 0.0);
    for (std::size_t f = 0; f < m.origin.size() / 3; ++f) {
        const auto l = sides(t, f);
        const double a = area_from_sides(l[0], l[1], l[2]);
        if (a <= 0.0) continue;
        std::size_t obtuse = 3;  // 3 = none (a triangle has at most one obtuse corner)
        for (std::size_t k = 0; k < 3; ++k) {
            if (corner_numerator(l, k) < 0.0) obtuse = k;
        }
        for (std::size_t k = 0; k < 3; ++k) {
            const std::uint32_t v = m.origin[3 * f + k];
            if (obtuse < 3) {
                area[v] += (k == obtuse) ? a / 2 : a / 4;
            } else {
                // Voronoi part of corner P = k in PQR: (|PR|^2 cot Q + |PQ|^2 cot R) / 8.
                const double pq = l[k], pr = l[(k + 2) % 3];
                const double cot_q = corner_numerator(l, (k + 1) % 3) / (4 * a);
                const double cot_r = corner_numerator(l, (k + 2) % 3) / (4 * a);
                area[v] += (pr * pr * cot_q + pq * pq * cot_r) / 8.0;
            }
        }
    }
    return area;
}

std::vector<double> angle_sums(const IntrinsicTriangulation& t) {
    const HalfEdgeMesh& m = t.connectivity;
    std::vector<double> sum(m.positions.size(), 0.0);
    for (std::size_t f = 0; f < m.origin.size() / 3; ++f) {
        const auto l = sides(t, f);
        const double a = area_from_sides(l[0], l[1], l[2]);
        for (std::size_t k = 0; k < 3; ++k) {
            sum[m.origin[3 * f + k]] += std::atan2(4.0 * a, corner_numerator(l, k));
        }
    }
    return sum;
}

std::vector<Vec2> intrinsic_gradient(const IntrinsicTriangulation& t, std::span<const double> u) {
    const HalfEdgeMesh& m = t.connectivity;
    const std::size_t nf = m.origin.size() / 3;
    std::vector<Vec2> grad(nf, Vec2{0.0, 0.0});
    for (std::size_t f = 0; f < nf; ++f) {
        const auto p = layout(sides(t, f));
        const double twice_area = cross2(sub(p[1], p[0]), sub(p[2], p[0]));
        if (twice_area <= 0.0) continue;
        Vec2 g{0.0, 0.0};
        for (std::size_t i = 0; i < 3; ++i) {
            const Vec2 e = sub(p[(i + 2) % 3], p[(i + 1) % 3]);  // opposite corner i, counter-clockwise
            const double ui = u[m.origin[3 * f + i]];
            g[0] += ui * -e[1];  // N x e with N = +z rotates e by +90 degrees: (-e_y, e_x)
            g[1] += ui * e[0];
        }
        grad[f] = {g[0] / twice_area, g[1] / twice_area};
    }
    return grad;
}

std::vector<double> intrinsic_divergence(const IntrinsicTriangulation& t, std::span<const Vec2> field) {
    const HalfEdgeMesh& m = t.connectivity;
    std::vector<double> div(m.positions.size(), 0.0);
    for (std::size_t f = 0; f < m.origin.size() / 3; ++f) {
        const auto p = layout(sides(t, f));
        const Vec2& x = field[f];
        for (std::size_t i = 0; i < 3; ++i) {
            const std::size_t j = (i + 1) % 3, k = (i + 2) % 3;
            const Vec2 e1 = sub(p[j], p[i]), e2 = sub(p[k], p[i]);
            div[m.origin[3 * f + i]] +=
                0.5 * (cot2(p[k], p[i], p[j]) * dot2(e1, x) + cot2(p[j], p[k], p[i]) * dot2(e2, x));
        }
    }
    return div;
}

}  // namespace dmw
