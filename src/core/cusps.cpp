#include "core/cusps.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>

#include "core/curvature.h"
#include "detail/hash.h"
#include "detail/vec.h"

namespace dmw {
namespace {

using namespace detail;
constexpr double kNegInf = -std::numeric_limits<double>::infinity();

// Eigen-decomposition of a symmetric 3x3 matrix by cyclic Jacobi rotations: each rotation zeroes
// one off-diagonal entry; a few sweeps converge to machine precision for 3x3.
void symmetric_eigen3(double a[3][3], double vec[3][3]) {
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) vec[i][j] = (i == j) ? 1.0 : 0.0;
    for (int sweep = 0; sweep < 50; ++sweep) {
        const double off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
        if (off < 1e-30) break;
        for (int p = 0; p < 2; ++p) {
            for (int q = p + 1; q < 3; ++q) {
                if (std::abs(a[p][q]) < 1e-300) continue;
                const double theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                const double t = (theta >= 0 ? 1.0 : -1.0) / (std::abs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
                for (int k = 0; k < 3; ++k) {  // A <- J^T A J
                    const double akp = a[k][p], akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for (int k = 0; k < 3; ++k) {
                    const double apk = a[p][k], aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for (int k = 0; k < 3; ++k) {  // V <- V J
                    const double vkp = vec[k][p], vkq = vec[k][q];
                    vec[k][p] = c * vkp - s * vkq;
                    vec[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }
}

}  // namespace

std::vector<std::uint8_t> largest_component_mask(const HalfEdgeMesh& m) {
    const std::size_t nv = m.positions.size();
    std::vector<std::uint32_t> comp(nv, kInvalid), size;
    std::vector<std::uint32_t> queue;
    for (std::uint32_t v = 0; v < nv; ++v) {
        if (comp[v] != kInvalid || m.vertex_halfedge[v] == kInvalid) continue;
        const auto id = static_cast<std::uint32_t>(size.size());
        comp[v] = id;
        queue.assign(1, v);
        for (std::size_t q = 0; q < queue.size(); ++q) {
            for (std::uint32_t w : one_ring(m, queue[q])) {
                if (comp[w] == kInvalid) comp[w] = id, queue.push_back(w);
            }
        }
        size.push_back(static_cast<std::uint32_t>(queue.size()));
    }
    std::vector<std::uint8_t> mask(nv, 0);
    if (size.empty()) return mask;
    const auto best = static_cast<std::uint32_t>(std::max_element(size.begin(), size.end()) - size.begin());
    for (std::size_t v = 0; v < nv; ++v) mask[v] = (comp[v] == best);
    return mask;
}

std::vector<std::uint32_t> local_maxima(const HalfEdgeMesh& m, std::span<const double> score, double threshold) {
    std::vector<std::uint32_t> out;
    for (std::uint32_t v = 0; v < m.positions.size(); ++v) {
        if (!(score[v] > threshold) || m.vertex_halfedge[v] == kInvalid) continue;
        bool is_max = true;
        for (std::uint32_t w : one_ring(m, v)) {
            if (!(score[v] > score[w])) {
                is_max = false;
                break;
            }
        }
        if (is_max) out.push_back(v);
    }
    return out;
}

Vec3 occlusal_axis(const HalfEdgeMesh& m, std::span<const std::uint8_t> mask) {
    Vec3 c{};
    double n = 0.0;
    for (std::size_t v = 0; v < m.positions.size(); ++v) {
        if (mask[v]) c += m.positions[v], n += 1.0;
    }
    if (n == 0.0) return {0, 0, 1};
    c = (1.0 / n) * c;
    double cov[3][3] = {};
    for (std::size_t v = 0; v < m.positions.size(); ++v) {
        if (!mask[v]) continue;
        const Vec3 d = m.positions[v] - c;
        const double e[3] = {d.x, d.y, d.z};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) cov[i][j] += e[i] * e[j];
    }
    double vec[3][3];
    symmetric_eigen3(cov, vec);
    int smallest = 0;
    for (int i = 1; i < 3; ++i) {
        if (cov[i][i] < cov[smallest][smallest]) smallest = i;
    }
    Vec3 axis{vec[0][smallest], vec[1][smallest], vec[2][smallest]};
    // Sign: away from the boundary (the scan's cut through the gingiva) toward the crowns.
    Vec3 cb{};
    double nb = 0.0;
    for (std::uint32_t h = 0; h < m.origin.size(); ++h) {
        if (m.twin[h] == kInvalid && mask[m.origin[h]]) cb += m.positions[m.origin[h]], nb += 1.0;
    }
    if (nb > 0.0 && dot(axis, c - (1.0 / nb) * cb) < 0.0) axis = -1.0 * axis;
    return (1.0 / norm(axis)) * axis;
}

std::vector<double> diffuse(const DecOperators& ops, std::span<const double> field, double sigma, double tolerance) {
    const double t = 0.5 * sigma * sigma;  // heat kernel standard deviation sqrt(2t) = sigma
    std::vector<double> mass = ops.star0;
    for (double& a : mass) a = a > 0.0 ? a : 1.0;  // vertices without faces: decoupled
    const SparseMatrix system = scaled_plus_diagonal(ops.laplacian, -t, mass);
    std::vector<double> rhs(field.size());
    for (std::size_t v = 0; v < field.size(); ++v) rhs[v] = mass[v] * field[v];
    return conjugate_gradient(system, rhs, tolerance).x;
}

std::vector<double> prominence(const DecOperators& ops, std::span<const double> field, double sigma, double tolerance) {
    const double t = 0.5 * sigma * sigma;
    std::vector<double> mass = ops.star0;
    for (double& a : mass) a = a > 0.0 ? a : 1.0;
    const SparseMatrix system = scaled_plus_diagonal(ops.laplacian, -t, mass);
    std::vector<double> rhs = ops.laplacian.multiply(field);
    for (double& r : rhs) r *= -t;
    return conjugate_gradient(system, rhs, tolerance).x;
}

std::vector<std::uint32_t> non_max_suppression(std::span<const Vec3> positions, std::span<const double> score,
                                               std::span<const std::uint32_t> candidates, double radius) {
    std::vector<std::uint32_t> order(candidates.begin(), candidates.end());
    std::sort(order.begin(), order.end(), [&](std::uint32_t a, std::uint32_t b) { return score[a] > score[b]; });
    // Uniform grid with cell size = radius: a kept point within `radius` lies in the 27 cells around.
    auto cell = [radius](const Vec3& p) {
        return std::array<std::int64_t, 3>{static_cast<std::int64_t>(std::floor(p.x / radius)),
                                           static_cast<std::int64_t>(std::floor(p.y / radius)),
                                           static_cast<std::int64_t>(std::floor(p.z / radius))};
    };
    std::unordered_map<std::array<std::int64_t, 3>, std::vector<std::uint32_t>, Hash3<std::int64_t>> grid;
    std::vector<std::uint32_t> kept;
    const double r2 = radius * radius;
    for (std::uint32_t v : order) {
        const Vec3& p = positions[v];
        const auto c = cell(p);
        bool suppressed = false;
        for (std::int64_t dx = -1; dx <= 1 && !suppressed; ++dx)
            for (std::int64_t dy = -1; dy <= 1 && !suppressed; ++dy)
                for (std::int64_t dz = -1; dz <= 1 && !suppressed; ++dz) {
                    const auto it = grid.find({c[0] + dx, c[1] + dy, c[2] + dz});
                    if (it == grid.end()) continue;
                    for (std::uint32_t k : it->second) {
                        const Vec3 d = positions[k] - p;
                        if (dot(d, d) < r2) {
                            suppressed = true;
                            break;
                        }
                    }
                }
        if (suppressed) continue;
        kept.push_back(v);
        grid[c].push_back(v);
    }
    return kept;
}

// CG tolerance for the detector's two smoothing solves (D70): measured on a 93.6k-vertex scan, 1e-4
// gives a max prominence error of 5e-4 mm (threshold 1.3 mm; scanner accuracy 10-90 um) at 4x less
// time than 1e-8. Valid because prominence is solved directly (small right-hand side).
constexpr double kDetectorTolerance = 1e-4;

CuspDetection detect_cusps(const HalfEdgeMesh& m, const CuspParams& p) {
    const std::size_t nv = m.positions.size();
    const std::vector<std::uint8_t> arch = largest_component_mask(m);
    CuspDetection out;
    out.occlusal_axis = occlusal_axis(m, arch);

    const CurvatureField curv = compute_curvature(m);
    std::vector<double> mean(nv);
    for (std::size_t v = 0; v < nv; ++v) mean[v] = std::isfinite(curv.mean[v]) ? curv.mean[v] : 0.0;

    std::vector<double> score(nv, kNegInf);
    if (p.method == CuspParams::Method::RawCurvature) {
        for (std::size_t v = 0; v < nv; ++v) {
            if (arch[v] && std::isfinite(curv.mean[v])) score[v] = curv.mean[v];
        }
    } else {
        const DecOperators ops = build_dec(m);
        const std::vector<double> smooth_h = diffuse(ops, mean, p.curvature_scale, kDetectorTolerance);
        if (p.method == CuspParams::Method::SmoothedCurvature) {
            for (std::size_t v = 0; v < nv; ++v) {
                if (arch[v]) score[v] = smooth_h[v];
            }
        } else {
            // Height along the occlusal axis, and its local baseline at scale R.
            std::vector<double> h(nv, 0.0);
            std::vector<double> arch_heights;
            for (std::size_t v = 0; v < nv; ++v) {
                h[v] = dot(m.positions[v], out.occlusal_axis);
                if (arch[v]) arch_heights.push_back(h[v]);
            }
            std::sort(arch_heights.begin(), arch_heights.end());
            const double gate = arch_heights.empty()
                                    ? 0.0
                                    : arch_heights[static_cast<std::size_t>(p.height_quantile * double(arch_heights.size() - 1))];
            const std::vector<double> prom = prominence(ops, h, p.prominence_scale, kDetectorTolerance);
            for (std::size_t v = 0; v < nv; ++v) {
                if (arch[v] && smooth_h[v] > 0.0 && h[v] >= gate) score[v] = prom[v];
            }
        }
    }
    const auto candidates = local_maxima(m, score, p.threshold);
    out.vertices = non_max_suppression(m.positions, score, candidates, p.nms_radius);
    for (std::uint32_t v : out.vertices) out.scores.push_back(score[v]);
    return out;
}

}  // namespace dmw
