#include "core/dec.h"
#include "core/intrinsic.h"

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

// The cotan Laplacian is intrinsic: build everything from the mesh's Euclidean edge lengths,
// with the one implementation in intrinsic.cpp (D52).
DecOperators build_dec(const HalfEdgeMesh& m) { return build_dec(intrinsic_from_mesh(m)); }

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
