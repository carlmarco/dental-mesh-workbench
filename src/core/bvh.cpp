#include "core/bvh.h"

#include <algorithm>
#include <cmath>

#include "detail/vec.h"

namespace dmw {
namespace {

using namespace detail;

// Moller-Trumbore: hit of origin + t dir with the triangle (v0, v0 + e1, v0 + e2); fills t, u, v on a hit.
inline bool hit_triangle(const Vec3& o, const Vec3& d, const Vec3& v0, const Vec3& e1, const Vec3& e2, double t_min, double t_max,
                         double& t, double& u, double& v) {
    const Vec3 p = cross(d, e2);
    const double det = dot(e1, p);
    if (std::abs(det) < 1e-18) return false;  // ray parallel to the triangle's plane (or a degenerate triangle)
    const double inv = 1.0 / det;
    const Vec3 s = o - v0;
    // Barycentric bounds inflated by kEdge (dimensionless): a ray through a shared edge or vertex hits at least one of
    // the triangles instead of slipping between them by rounding (D92: offset shells align vertices along normals).
    constexpr double kEdge = 1e-10;
    u = dot(s, p) * inv;
    if (u < -kEdge || u > 1.0 + kEdge) return false;
    const Vec3 q = cross(s, e1);
    v = dot(d, q) * inv;
    if (v < -kEdge || u + v > 1.0 + kEdge) return false;
    t = dot(e2, q) * inv;
    return t > t_min && t < t_max;
}

struct Box {
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    void grow(const Vec3& p) {
        lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
        hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
    }
    void grow(const Box& b) { grow(b.lo), grow(b.hi); }
    double area() const {
        if (hi.x < lo.x) return 0.0;
        const Vec3 e = hi - lo;
        return 2.0 * (e.x * e.y + e.y * e.z + e.z * e.x);
    }
};
inline double axis(const Vec3& p, int a) { return a == 0 ? p.x : a == 1 ? p.y : p.z; }

// Entry distance of the ray into the box, or +inf if it misses within [t_min, t_max]. Axes with a zero direction
// component are handled explicitly: the ray is inside that slab for all t or for none (computing (lo - o) * inf
// would give NaN = 0 * inf for an origin exactly on the slab plane; D92 found this with axis-aligned rays).
inline double enter(const Vec3& lo, const Vec3& hi, const Vec3& o, const Vec3& inv, double t_min, double t_max) {
    double t0 = t_min, t1 = t_max;
    const double ol[3] = {lo.x, lo.y, lo.z}, oh[3] = {hi.x, hi.y, hi.z}, oo[3] = {o.x, o.y, o.z}, iv[3] = {inv.x, inv.y, inv.z};
    for (int k = 0; k < 3; ++k) {
        if (std::isinf(iv[k])) {  // direction component 0
            if (oo[k] < ol[k] || oo[k] > oh[k]) return std::numeric_limits<double>::infinity();
            continue;
        }
        double a = (ol[k] - oo[k]) * iv[k], b = (oh[k] - oo[k]) * iv[k];
        if (a > b) std::swap(a, b);
        t0 = std::max(t0, a), t1 = std::min(t1, b);
        if (t0 > t1) return std::numeric_limits<double>::infinity();
    }
    return t0;
}

}  // namespace

Bvh::Bvh(std::span<const Vec3> positions, std::span<const std::array<std::uint32_t, 3>> triangles, std::uint32_t leaf_size) {
    const auto n = static_cast<std::uint32_t>(triangles.size());
    v0_.resize(n), e1_.resize(n), e2_.resize(n), order_.resize(n);
    std::vector<Box> box(n);
    std::vector<Vec3> centroid(n);
    for (std::uint32_t f = 0; f < n; ++f) {
        const Vec3 &a = positions[triangles[f][0]], &b = positions[triangles[f][1]], &c = positions[triangles[f][2]];
        v0_[f] = a, e1_[f] = b - a, e2_[f] = c - a, order_[f] = f;
        box[f].grow(a), box[f].grow(b), box[f].grow(c);
        centroid[f] = (1.0 / 3.0) * (a + b + c);
    }
    if (n == 0) return;
    nodes_.reserve(2 * static_cast<std::size_t>(n));
    nodes_.push_back({{}, {}, 0, n});
    std::vector<std::uint32_t> todo{0};
    constexpr int kBins = 12;
    while (!todo.empty()) {
        const std::uint32_t ni = todo.back();
        todo.pop_back();
        const std::uint32_t first = nodes_[ni].first, count = nodes_[ni].count;
        Box bounds, cbounds;
        for (std::uint32_t i = first; i < first + count; ++i) bounds.grow(box[order_[i]]), cbounds.grow(centroid[order_[i]]);
        nodes_[ni].lo = bounds.lo, nodes_[ni].hi = bounds.hi;
        if (count <= leaf_size) continue;  // leaf
        // Binned SAH over all three axes: cost(split) = N_left A_left + N_right A_right.
        int best_axis = -1, best_split = 0;
        double best_cost = std::numeric_limits<double>::infinity();
        for (int a = 0; a < 3; ++a) {
            const double lo = axis(cbounds.lo, a), extent = axis(cbounds.hi, a) - lo;
            if (extent <= 0.0) continue;
            Box bin_box[kBins];
            std::uint32_t bin_count[kBins] = {};
            for (std::uint32_t i = first; i < first + count; ++i) {
                const std::uint32_t f = order_[i];
                const int b = std::min(kBins - 1, static_cast<int>(kBins * (axis(centroid[f], a) - lo) / extent));
                bin_box[b].grow(box[f]), ++bin_count[b];
            }
            double right_area[kBins];
            std::uint32_t right_count[kBins];
            Box acc;
            std::uint32_t cnt = 0;
            for (int b = kBins - 1; b > 0; --b) acc.grow(bin_box[b]), cnt += bin_count[b], right_area[b] = acc.area(), right_count[b] = cnt;
            acc = Box{}, cnt = 0;
            for (int b = 0; b < kBins - 1; ++b) {
                acc.grow(bin_box[b]), cnt += bin_count[b];
                const double cost = double(cnt) * acc.area() + double(right_count[b + 1]) * right_area[b + 1];
                if (cnt > 0 && right_count[b + 1] > 0 && cost < best_cost) best_cost = cost, best_axis = a, best_split = b;
            }
        }
        std::uint32_t mid;
        if (best_axis < 0) {
            mid = first + count / 2;  // all centroids coincide: split the list in half
        } else {
            const double lo = axis(cbounds.lo, best_axis), extent = axis(cbounds.hi, best_axis) - lo;
            const auto begin = order_.begin() + static_cast<std::ptrdiff_t>(first);
            const auto it = std::partition(begin, begin + static_cast<std::ptrdiff_t>(count), [&](std::uint32_t f) {
                return std::min(kBins - 1, static_cast<int>(kBins * (axis(centroid[f], best_axis) - lo) / extent)) <= best_split;
            });
            mid = static_cast<std::uint32_t>(it - order_.begin());
        }
        const auto left = static_cast<std::uint32_t>(nodes_.size());
        nodes_.push_back({{}, {}, first, mid - first});
        nodes_.push_back({{}, {}, mid, first + count - mid});
        nodes_[ni].first = left, nodes_[ni].count = 0;
        todo.push_back(left), todo.push_back(left + 1);
    }
}

template <bool kAny>
RayHit Bvh::trace(const Vec3& o, const Vec3& d, double t_min, double t_max, std::uint32_t ignore) const {
    RayHit best;
    if (nodes_.empty()) return best;
    const Vec3 inv{1.0 / d.x, 1.0 / d.y, 1.0 / d.z};
    std::uint32_t stack[64];
    int top = 0;
    if (enter(nodes_[0].lo, nodes_[0].hi, o, inv, t_min, t_max) == std::numeric_limits<double>::infinity()) return best;
    stack[top++] = 0;
    while (top > 0) {
        const Node& node = nodes_[stack[--top]];
        if (node.count > 0) {
            for (std::uint32_t i = node.first; i < node.first + node.count; ++i) {
                const std::uint32_t f = order_[i];
                double t, u, v;
                if (f == ignore || !hit_triangle(o, d, v0_[f], e1_[f], e2_[f], t_min, t_max, t, u, v)) continue;
                best = {true, t, f, u, v};
                if constexpr (kAny) return best;
                t_max = t;  // shrink: only nearer hits matter now
            }
            continue;
        }
        const double tl = enter(nodes_[node.first].lo, nodes_[node.first].hi, o, inv, t_min, t_max);
        const double tr = enter(nodes_[node.first + 1].lo, nodes_[node.first + 1].hi, o, inv, t_min, t_max);
        const bool hl = tl != std::numeric_limits<double>::infinity(), hr = tr != std::numeric_limits<double>::infinity();
        // Push the farther child first so the nearer one is popped (and shrinks t_max) first.
        if (hl && hr) {
            stack[top++] = tl <= tr ? node.first + 1 : node.first;
            stack[top++] = tl <= tr ? node.first : node.first + 1;
        } else if (hl) {
            stack[top++] = node.first;
        } else if (hr) {
            stack[top++] = node.first + 1;
        }
    }
    return best;
}

Vec3 closest_point_on_triangle(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c) {
    const Vec3 ab = b - a, ac = c - a, ap = p - a;
    const double d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) return a;  // vertex region a
    const Vec3 bp = p - b;
    const double d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) return b;  // vertex region b
    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) return a + (d1 / (d1 - d3)) * ab;  // edge ab
    const Vec3 cp = p - c;
    const double d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) return c;  // vertex region c
    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) return a + (d2 / (d2 - d6)) * ac;  // edge ac
    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) return b + ((d4 - d3) / ((d4 - d3) + (d5 - d6))) * (c - b);  // edge bc
    const double denom = 1.0 / (va + vb + vc);
    return a + (vb * denom) * ab + (vc * denom) * ac;  // face interior
}

namespace {
double box_distance2(const Vec3& lo, const Vec3& hi, const Vec3& p) {
    const double dx = std::max({lo.x - p.x, 0.0, p.x - hi.x}), dy = std::max({lo.y - p.y, 0.0, p.y - hi.y}),
                 dz = std::max({lo.z - p.z, 0.0, p.z - hi.z});
    return dx * dx + dy * dy + dz * dz;
}
}  // namespace

ClosestPoint Bvh::closest(const Vec3& p) const {
    ClosestPoint best;
    if (nodes_.empty()) return best;
    double best2 = std::numeric_limits<double>::infinity();
    std::uint32_t stack[64];
    int top = 0;
    stack[top++] = 0;
    while (top > 0) {
        const Node& node = nodes_[stack[--top]];
        if (box_distance2(node.lo, node.hi, p) >= best2) continue;
        if (node.count > 0) {
            for (std::uint32_t i = node.first; i < node.first + node.count; ++i) {
                const std::uint32_t f = order_[i];
                const Vec3 q = closest_point_on_triangle(p, v0_[f], v0_[f] + e1_[f], v0_[f] + e2_[f]);
                const Vec3 d = q - p;
                const double d2 = dot(d, d);
                if (d2 < best2) best2 = d2, best = {q, 0.0, f};
            }
            continue;
        }
        // Push the farther child first, so the nearer one is explored (and tightens the bound) first.
        const std::uint32_t l = node.first, r = node.first + 1;
        const double dl = box_distance2(nodes_[l].lo, nodes_[l].hi, p), dr = box_distance2(nodes_[r].lo, nodes_[r].hi, p);
        if (dl <= dr) {
            if (dr < best2) stack[top++] = r;
            if (dl < best2) stack[top++] = l;
        } else {
            if (dl < best2) stack[top++] = l;
            if (dr < best2) stack[top++] = r;
        }
    }
    best.distance = std::sqrt(best2);
    return best;
}

ClosestPoint closest_brute_force(std::span<const Vec3> pos, std::span<const std::array<std::uint32_t, 3>> tri, const Vec3& p) {
    ClosestPoint best;
    double best2 = std::numeric_limits<double>::infinity();
    for (std::uint32_t f = 0; f < tri.size(); ++f) {
        const Vec3 q = closest_point_on_triangle(p, pos[tri[f][0]], pos[tri[f][1]], pos[tri[f][2]]);
        const Vec3 d = q - p;
        if (dot(d, d) < best2) best2 = dot(d, d), best = {q, 0.0, f};
    }
    best.distance = std::sqrt(best2);
    return best;
}

RayHit Bvh::intersect(const Vec3& o, const Vec3& d, double t_min, double t_max, std::uint32_t ignore) const {
    return trace<false>(o, d, t_min, t_max, ignore);
}

bool Bvh::occluded(const Vec3& o, const Vec3& d, double t_min, double t_max, std::uint32_t ignore) const {
    return trace<true>(o, d, t_min, t_max, ignore).hit;
}

RayHit intersect_brute_force(std::span<const Vec3> p, std::span<const std::array<std::uint32_t, 3>> tri, const Vec3& o, const Vec3& d,
                             double t_min, double t_max, std::uint32_t ignore) {
    RayHit best;
    for (std::uint32_t f = 0; f < tri.size(); ++f) {
        if (f == ignore) continue;
        const Vec3 &a = p[tri[f][0]], &b = p[tri[f][1]], &c = p[tri[f][2]];
        double t, u, v;
        if (hit_triangle(o, d, a, b - a, c - a, t_min, t_max, t, u, v)) best = {true, t, f, u, v}, t_max = t;
    }
    return best;
}

}  // namespace dmw
