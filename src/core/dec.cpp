#include "core/dec.h"

#include <algorithm>
#include <cstdint>
#include <vector>

#include "detail/vec.h"

namespace dmw {
namespace {

using namespace detail;

// Cotangent of the angle at corner p between edges to q and r: (u.w) / |u x w|.
// Zero-area corners return 0 (their infinite cotangent is dropped, as in M4's D38).
double cot_at(const Vec3& p, const Vec3& q, const Vec3& r) {
    const Vec3 u = q - p, w = r - p;
    const double s = norm(cross(u, w));
    return s > 0.0 ? dot(u, w) / s : 0.0;
}

}  // namespace

DecOperators build_dec(const HalfEdgeMesh& m) {
    const auto nv = static_cast<std::uint32_t>(m.positions.size());
    const auto nh = static_cast<std::uint32_t>(m.origin.size());
    const std::uint32_t nf = nh / 3;
    DecOperators ops;

    // Edges from half-edges: one per twin pair (taken at the smaller id) or per boundary
    // half-edge. Each edge is oriented low -> high vertex index.
    std::vector<std::uint32_t> edge_of(nh, kInvalid);
    double total_length = 0.0;
    for (std::uint32_t h = 0; h < nh; ++h) {
        if (edge_of[h] != kInvalid) continue;
        const std::uint32_t a = m.origin[h], b = dest(m, h);
        const auto e = static_cast<std::uint32_t>(ops.edges.size());
        ops.edges.push_back({std::min(a, b), std::max(a, b)});
        edge_of[h] = e;
        if (m.twin[h] != kInvalid) edge_of[m.twin[h]] = e;
        total_length += norm(m.positions[a] - m.positions[b]);
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
        const Vec3 p[3] = {m.positions[m.origin[3 * f]], m.positions[m.origin[3 * f + 1]],
                           m.positions[m.origin[3 * f + 2]]};
        const double area = 0.5 * norm(cross(p[1] - p[0], p[2] - p[0]));
        for (std::uint32_t k = 0; k < 3; ++k) {
            const std::uint32_t h = 3 * f + k, e = edge_of[h];
            // Half-edge h runs corner k -> k+1, agreeing with edge e iff it starts at e.v0.
            d1.push_back({f, e, m.origin[h] == ops.edges[e].v0 ? 1.0 : -1.0});
            // The angle opposite h is at corner k+2.
            ops.star1[e] += 0.5 * cot_at(p[(k + 2) % 3], p[k], p[(k + 1) % 3]);
            ops.star0[m.origin[h]] += area / 3.0;
        }
    }
    ops.d0 = SparseMatrix::from_triplets(ne, nv, std::move(d0));
    ops.d1 = SparseMatrix::from_triplets(nf, ne, std::move(d1));

    // L = -d0^T *1 d0. Expanding one edge (i, j) with weight w: L_ii -= w, L_jj -= w,
    // L_ij += w, L_ji += w. Hence rows sum to zero (constants are in the kernel).
    ops.laplacian = weighted_gram(ops.d0, ops.star1);
    for (double& v : ops.laplacian.value) v = -v;
    return ops;
}

std::vector<Vec3> face_gradient(const HalfEdgeMesh& m, std::span<const double> u) {
    const std::size_t nf = m.origin.size() / 3;
    std::vector<Vec3> grad(nf, Vec3{});
    for (std::size_t f = 0; f < nf; ++f) {
        const std::uint32_t v[3] = {m.origin[3 * f], m.origin[3 * f + 1], m.origin[3 * f + 2]};
        const Vec3 p[3] = {m.positions[v[0]], m.positions[v[1]], m.positions[v[2]]};
        const Vec3 n = cross(p[1] - p[0], p[2] - p[0]);  // = 2A * unit normal
        const double twice_area = norm(n);
        if (twice_area == 0.0) continue;
        const Vec3 unit = (1.0 / twice_area) * n;
        // N x e_i is the in-plane edge vector rotated 90 degrees inward, scaled by |e_i|;
        // u varies linearly, so summing u_i (N x e_i) / 2A recovers its constant gradient.
        Vec3 g{};
        for (int i = 0; i < 3; ++i) {
            const Vec3 e = p[(i + 2) % 3] - p[(i + 1) % 3];  // opposite corner i, counter-clockwise
            g += u[v[i]] * cross(unit, e);
        }
        grad[f] = (1.0 / twice_area) * g;
    }
    return grad;
}

std::vector<double> vertex_divergence(const HalfEdgeMesh& m, std::span<const Vec3> field) {
    std::vector<double> div(m.positions.size(), 0.0);
    const std::size_t nf = m.origin.size() / 3;
    for (std::size_t f = 0; f < nf; ++f) {
        const std::uint32_t v[3] = {m.origin[3 * f], m.origin[3 * f + 1], m.origin[3 * f + 2]};
        const Vec3 p[3] = {m.positions[v[0]], m.positions[v[1]], m.positions[v[2]]};
        const Vec3& x = field[f];
        for (int i = 0; i < 3; ++i) {
            const int j = (i + 1) % 3, k = (i + 2) % 3;
            const Vec3 e1 = p[j] - p[i], e2 = p[k] - p[i];  // edges leaving corner i
            // Angle opposite e1 is at corner k; opposite e2 is at corner j.
            div[v[i]] += 0.5 * (cot_at(p[k], p[i], p[j]) * dot(e1, x) + cot_at(p[j], p[k], p[i]) * dot(e2, x));
        }
    }
    return div;
}

}  // namespace dmw
