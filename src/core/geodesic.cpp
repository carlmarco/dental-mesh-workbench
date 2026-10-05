#include "core/geodesic.h"

#include <cmath>
#include <functional>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

#include "detail/vec.h"

namespace dmw {
namespace {

using namespace detail;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// Vertices reachable from the sources through mesh edges (breadth-first over one-rings).
std::vector<bool> reachable(const HalfEdgeMesh& m, std::span<const std::uint32_t> sources) {
    std::vector<bool> seen(m.positions.size(), false);
    std::vector<std::uint32_t> queue(sources.begin(), sources.end());
    for (std::uint32_t s : sources) seen[s] = true;
    for (std::size_t q = 0; q < queue.size(); ++q) {
        for (std::uint32_t w : one_ring(m, queue[q])) {
            if (!seen[w]) seen[w] = true, queue.push_back(w);
        }
    }
    return seen;
}

std::string validate(const HalfEdgeMesh& m, std::span<const std::uint32_t> sources) {
    if (sources.empty()) return "no source vertices";
    for (std::uint32_t s : sources) {
        if (s >= m.positions.size()) return "source vertex " + std::to_string(s) + " out of range";
    }
    return {};
}

}  // namespace

namespace {

// t = m h^2 (paper, Sec. 3.2.4): diffuse about one edge length. Larger m smooths.
double time_step_for(const DecOperators& ops, double time_factor) {
    return time_factor * ops.mean_edge_length * ops.mean_edge_length;
}

// *0 - t L, with a unit diagonal for vertices that have no faces (*0 = 0 would be singular;
// such a vertex is decoupled and never receives heat unless it is itself a source).
SparseMatrix heat_matrix(const DecOperators& ops, double t) {
    std::vector<double> diag = ops.star0;
    for (double& d : diag) d = d > 0.0 ? d : 1.0;
    return scaled_plus_diagonal(ops.laplacian, -t, diag);
}

// -L is singular (constants per component). Pin the lowest-index vertex of each component.
std::vector<std::uint32_t> one_vertex_per_component(const HalfEdgeMesh& m) {
    const std::size_t nv = m.positions.size();
    std::vector<std::uint32_t> pinned, queue;
    std::vector<bool> seen(nv, false);
    for (std::uint32_t v = 0; v < nv; ++v) {
        if (seen[v]) continue;
        pinned.push_back(v);  // first vertex reached in this component
        seen[v] = true;
        queue.assign(1, v);
        for (std::size_t q = 0; q < queue.size(); ++q) {
            for (std::uint32_t w : one_ring(m, queue[q])) {
                if (!seen[w]) seen[w] = true, queue.push_back(w);
            }
        }
    }
    return pinned;
}

// Replace each pinned vertex's row and column of -L by the identity. Dropping column p from
// the other rows assumes phi_p = 0, so distance() must also set rhs_p = 0. The system is
// consistent, so this yields an exact solution; distance() shifts it afterwards anyway.
SparseMatrix pinned_poisson(const DecOperators& ops, std::span<const std::uint32_t> pinned_ids) {
    const SparseMatrix& L = ops.laplacian;
    std::vector<bool> pinned(L.rows, false);
    for (std::uint32_t p : pinned_ids) pinned[p] = true;
    std::vector<Triplet> t;
    for (std::uint32_t r = 0; r < L.rows; ++r) {
        if (pinned[r]) {
            t.push_back({r, r, 1.0});
            continue;
        }
        for (std::uint32_t k = L.row_start[r]; k < L.row_start[r + 1]; ++k) {
            if (!pinned[L.col[k]]) t.push_back({r, L.col[k], -L.value[k]});
        }
    }
    return SparseMatrix::from_triplets(L.rows, L.cols, std::move(t));
}

}  // namespace

HeatGeodesics::HeatGeodesics(const HalfEdgeMesh& mesh, double time_factor)
    : mesh_(&mesh),
      ops_(build_dec(mesh)),
      time_step_(time_step_for(ops_, time_factor)),
      heat_(heat_matrix(ops_, time_step_)),
      pinned_(one_vertex_per_component(mesh)),
      poisson_(pinned_poisson(ops_, pinned_)) {
    if (!heat_.ok()) error_ = "heat matrix: " + heat_.error();
    else if (!poisson_.ok()) error_ = "Poisson matrix: " + poisson_.error();
}

GeodesicResult HeatGeodesics::distance(std::span<const std::uint32_t> sources) const {
    const HalfEdgeMesh& m = *mesh_;
    GeodesicResult res;
    res.time_step = time_step_;
    res.error = validate(m, sources);
    if (!res.ok()) return res;
    const std::size_t nv = m.positions.size();

    if (!ok()) {
        res.error = error_;
        return res;
    }

    // I. Heat step. The right-hand side is a Kronecker delta (integrated, like *0 u).
    std::vector<double> delta(nv, 0.0);
    for (std::uint32_t s : sources) delta[s] = 1.0;
    const std::vector<double> u = heat_.solve(delta);

    // II. Normalized negative gradient. Where heat never arrived (other components), grad u = 0
    // and X stays 0; those vertices are set to NaN below.
    std::vector<Vec3> x = face_gradient(m, u);
    for (Vec3& g : x) {
        const double len = norm(g);
        g = len > 0.0 ? (-1.0 / len) * g : Vec3{};
    }

    // III. Poisson: L phi = div X  <=>  (-L) phi = -div X. Consistent because the discrete
    // divergence sums to zero on each component. Pinned vertices get rhs 0 (phi = 0 there).
    std::vector<double> rhs = vertex_divergence(m, x);
    for (double& v : rhs) v = -v;
    for (std::uint32_t p : pinned_) rhs[p] = 0.0;
    const std::vector<double> phi = poisson_.solve(rhs);
    // phi is determined up to a constant per component: shift so the sources average to 0.
    double shift = 0.0;
    for (std::uint32_t s : sources) shift += phi[s];
    shift /= static_cast<double>(sources.size());
    const auto seen = reachable(m, sources);
    res.distance.assign(nv, kNaN);
    for (std::size_t v = 0; v < nv; ++v) {
        if (seen[v]) res.distance[v] = phi[v] - shift;
    }
    return res;
}

std::vector<double> dijkstra_distance(const HalfEdgeMesh& m, std::span<const std::uint32_t> sources) {
    const std::size_t nv = m.positions.size();
    std::vector<double> dist(nv, std::numeric_limits<double>::infinity());
    using Item = std::pair<double, std::uint32_t>;  // (tentative distance, vertex)
    std::priority_queue<Item, std::vector<Item>, std::greater<>> heap;  // min-heap
    for (std::uint32_t s : sources) {
        if (s < nv) dist[s] = 0.0, heap.push({0.0, s});
    }
    while (!heap.empty()) {
        const auto [d, v] = heap.top();
        heap.pop();
        if (d > dist[v]) continue;  // stale entry: v was already settled with a shorter path
        for (std::uint32_t w : one_ring(m, v)) {
            const double nd = d + norm(m.positions[w] - m.positions[v]);
            if (nd < dist[w]) dist[w] = nd, heap.push({nd, w});
        }
    }
    for (double& d : dist) {
        if (std::isinf(d)) d = kNaN;
    }
    return dist;
}

}  // namespace dmw
