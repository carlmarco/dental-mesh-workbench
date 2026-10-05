#include "core/curvature.h"
#include "detail/vec.h"

#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>

namespace dmw {
namespace {

using namespace detail;

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kPi = std::numbers::pi;


}  // namespace

CurvatureField compute_curvature(const HalfEdgeMesh& m) {
    const std::size_t nv = m.positions.size();
    const std::size_t nf = m.origin.size() / 3;

    CurvatureField c;
    c.boundary.assign(nv, 0);
    c.mixed_area.assign(nv, 0.0);
    c.angle_defect.assign(nv, 0.0);
    c.gaussian.assign(nv, kNaN);
    c.mean_normal.assign(nv, Vec3{});
    c.mean.assign(nv, kNaN);
    c.k1.assign(nv, kNaN);
    c.k2.assign(nv, kNaN);

    for (std::uint32_t h = 0; h < m.origin.size(); ++h) {
        if (m.twin[h] == kInvalid) c.boundary[m.origin[h]] = c.boundary[dest(m, h)] = 1;
    }

    std::vector<double> angle_sum(nv, 0.0);
    std::vector<Vec3> cot_sum(nv, Vec3{});     // sum_j (cot a + cot b)(x_i - x_j)
    std::vector<Vec3> normal_sum(nv, Vec3{});  // sum of face normals scaled by 2 * area

    // Per-face scatter: each triangle adds its share to its three corners. Each edge (i, j)
    // lies in two triangles, each contributing the cotangent of the angle opposite the edge,
    // so after all faces cot_sum[i] holds exactly Meyer's one-ring sum.
    for (std::size_t f = 0; f < nf; ++f) {
        const std::uint32_t v[3] = {m.origin[3 * f], m.origin[3 * f + 1], m.origin[3 * f + 2]};
        const Vec3 p[3] = {m.positions[v[0]], m.positions[v[1]], m.positions[v[2]]};
        const Vec3 n = cross(p[1] - p[0], p[2] - p[0]);
        const double twice_area = norm(n);  // = |u x v| at every corner of this triangle

        // Corner k with edges u = p[k+1] - p[k], w = p[k+2] - p[k]:
        //   cos ~ u.w, sin ~ |u x w| = twice_area, so
        //   theta = atan2(twice_area, u.w)  (robust near 0 and pi, unlike acos(u.w / |u||w|))
        //   cot   = u.w / twice_area        (no trig at all)
        double uw[3], cot[3];
        for (int k = 0; k < 3; ++k) {
            const Vec3 u = p[(k + 1) % 3] - p[k], w = p[(k + 2) % 3] - p[k];
            uw[k] = dot(u, w);
            angle_sum[v[k]] += std::atan2(twice_area, uw[k]);
            normal_sum[v[k]] += n;
        }
        if (twice_area == 0.0) {  // D38: cotangents are infinite; angles above are still valid
            c.degenerate_faces.push_back(static_cast<std::uint32_t>(f));
            continue;
        }
        for (int k = 0; k < 3; ++k) cot[k] = uw[k] / twice_area;

        // Cotan Laplacian: the edge opposite corner k joins corners k+1 and k+2.
        for (int k = 0; k < 3; ++k) {
            const int i = (k + 1) % 3, j = (k + 2) % 3;
            const Vec3 d = cot[k] * (p[i] - p[j]);
            cot_sum[v[i]] += d;
            cot_sum[v[j]] += -1.0 * d;
        }

        // Mixed area (Meyer et al., Fig. 4). The angle at corner k is obtuse iff u.w < 0.
        const double area = 0.5 * twice_area;
        int obtuse = -1;
        for (int k = 0; k < 3; ++k) {
            if (uw[k] < 0.0) obtuse = k;  // a triangle has at most one obtuse angle
        }
        if (obtuse < 0) {
            // Voronoi region of corner P in triangle PQR: (|PR|^2 cot Q + |PQ|^2 cot R) / 8.
            for (int k = 0; k < 3; ++k) {
                const int q = (k + 1) % 3, r = (k + 2) % 3;
                const Vec3 pq = p[q] - p[k], pr = p[r] - p[k];
                c.mixed_area[v[k]] += (dot(pr, pr) * cot[q] + dot(pq, pq) * cot[r]) / 8.0;
            }
        } else {
            // Circumcenter lies outside the triangle, so the Voronoi split would go negative.
            for (int k = 0; k < 3; ++k) c.mixed_area[v[k]] += (k == obtuse) ? area / 2 : area / 4;
        }
    }

    for (std::uint32_t i = 0; i < nv; ++i) {
        if (m.vertex_halfedge[i] == kInvalid) continue;  // isolated: no faces, all NaN, defect 0
        c.angle_defect[i] = (c.boundary[i] ? kPi : 2 * kPi) - angle_sum[i];
        const double area = c.mixed_area[i];
        if (c.boundary[i] || area <= 0.0) continue;  // D36: pointwise values interior-only

        const double gauss = c.angle_defect[i] / area;
        const Vec3 kn = (1.0 / (2.0 * area)) * cot_sum[i];  // Meyer's K(x_i) = 2 H n
        const double nn = norm(normal_sum[i]);
        // D37: the sign of H needs a normal independent of K(x_i) (which vanishes on flat
        // regions). The area-weighted average of face normals follows the winding.
        const double mean = nn > 0.0 ? 0.5 * dot(kn, (1.0 / nn) * normal_sum[i]) : 0.0;
        const double disc = std::sqrt(std::max(mean * mean - gauss, 0.0));  // clamp roundoff

        c.gaussian[i] = gauss;
        c.mean_normal[i] = kn;
        c.mean[i] = mean;
        c.k1[i] = mean + disc;
        c.k2[i] = mean - disc;
    }
    return c;
}

}  // namespace dmw
