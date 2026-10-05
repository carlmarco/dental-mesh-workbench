#include "core/handles.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

#include "detail/disjoint_sets.h"
#include "detail/vec.h"

namespace dmw {
namespace {

using detail::DisjointSets;

// One id per undirected edge: the smaller half-edge id of a twin pair, or the lone boundary one.
std::uint32_t canonical(const HalfEdgeMesh& m, std::uint32_t h) {
    const std::uint32_t t = m.twin[h];
    return t == kInvalid ? h : std::min(h, t);
}

}  // namespace

std::vector<HandleLoop> handle_loops(const HalfEdgeMesh& m) {
    const auto nv = static_cast<std::uint32_t>(m.positions.size());
    const auto nh = static_cast<std::uint32_t>(m.origin.size());
    const std::uint32_t nf = nh / 3;
    auto edge_length = [&](std::uint32_t h) {
        return detail::norm(detail::operator-(m.positions[dest(m, h)], m.positions[m.origin[h]]));
    };

    // --- T: Dijkstra shortest-path forest, one basepoint (lowest index) per component.
    // parent_he[v] = the half-edge whose edge reached v (kInvalid at a basepoint).
    std::vector<double> dist(nv, std::numeric_limits<double>::infinity());
    std::vector<std::uint32_t> parent(nv, kInvalid), parent_he(nv, kInvalid), depth(nv, 0);
    using Item = std::pair<double, std::uint32_t>;
    std::priority_queue<Item, std::vector<Item>, std::greater<>> heap;
    auto relax = [&](std::uint32_t v, std::uint32_t w, std::uint32_t h) {
        const double nd = dist[v] + edge_length(h);
        if (nd < dist[w]) {
            dist[w] = nd, parent[w] = v, parent_he[w] = h, depth[w] = depth[v] + 1;
            heap.push({nd, w});
        }
    };
    for (std::uint32_t root = 0; root < nv; ++root) {
        if (m.vertex_halfedge[root] == kInvalid || std::isfinite(dist[root])) continue;
        dist[root] = 0.0;
        heap.push({0.0, root});
        while (!heap.empty()) {
            const auto [d, v] = heap.top();
            heap.pop();
            if (d > dist[v]) continue;
            // Walk v's outgoing half-edges (rho = twin(prev(h)), as in one_ring); at a boundary the
            // last neighbor is reached through the incoming boundary half-edge prev(h).
            const std::uint32_t start = m.vertex_halfedge[v];
            for (std::uint32_t h = start;;) {
                relax(v, dest(m, h), h);
                const std::uint32_t in = prev(h), t = m.twin[in];
                if (t == kInvalid) {
                    relax(v, m.origin[in], in);
                    break;
                }
                if (t == start) break;
                h = t;
            }
        }
    }
    std::vector<bool> in_tree(nh, false);  // indexed by canonical half-edge
    for (std::uint32_t v = 0; v < nv; ++v) {
        if (parent_he[v] != kInvalid) in_tree[canonical(m, parent_he[v])] = true;
    }

    // --- Dual nodes: faces 0..nf-1, then one virtual node per boundary loop (capping the hole).
    std::vector<std::uint32_t> loop_of(nh, kInvalid);
    std::uint32_t num_loops = 0;
    for (std::uint32_t h = 0; h < nh; ++h) {
        if (m.twin[h] != kInvalid || loop_of[h] != kInvalid) continue;
        // Follow the boundary: the next boundary half-edge leaves dest(h); a boundary vertex stores
        // exactly that half-edge (half-edge invariant, D25).
        for (std::uint32_t b = h; loop_of[b] == kInvalid; b = m.vertex_halfedge[dest(m, b)]) loop_of[b] = num_loops;
        ++num_loops;
    }

    // --- C: maximum spanning forest of the dual graph over non-tree edges, by loop length sigma.
    // Kruskal in descending sigma: an edge joining two dual components goes into C; one that
    // would close a dual cycle is a generator. Taking long loops into C leaves the short ones.
    struct Candidate {
        double sigma;
        std::uint32_t h;
    };
    std::vector<Candidate> candidates;
    for (std::uint32_t h = 0; h < nh; ++h) {
        if (canonical(m, h) != h || in_tree[h]) continue;
        candidates.push_back({dist[m.origin[h]] + edge_length(h) + dist[dest(m, h)], h});
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const Candidate& a, const Candidate& b) { return a.sigma > b.sigma; });
    DisjointSets dual(nf + num_loops);
    std::vector<std::uint32_t> generators;
    for (const Candidate& c : candidates) {
        const std::uint32_t t = m.twin[c.h];
        const std::uint32_t other = t != kInvalid ? face(t) : nf + loop_of[c.h];
        if (!dual.unite(face(c.h), other)) generators.push_back(c.h);
    }

    // --- Each generator edge (u, v) closes u -> lca -> v; the stem from the basepoint is dropped.
    std::vector<HandleLoop> loops;
    for (std::uint32_t h : generators) {
        std::uint32_t a = m.origin[h], b = dest(m, h);
        std::vector<std::uint32_t> left{a}, right{b};  // walking up from each end
        while (a != b) {
            if (depth[a] >= depth[b]) a = parent[a], left.push_back(a);
            else b = parent[b], right.push_back(b);
        }
        right.pop_back();  // the lca is already the last entry of `left`
        HandleLoop loop;
        loop.vertices = std::move(left);  // u ... lca
        loop.vertices.insert(loop.vertices.end(), right.rbegin(), right.rend());  // ... down to v
        for (std::size_t i = 0; i < loop.vertices.size(); ++i) {
            const Vec3& p = m.positions[loop.vertices[i]];
            const Vec3& q = m.positions[loop.vertices[(i + 1) % loop.vertices.size()]];
            loop.length += detail::norm(detail::operator-(p, q));
        }
        loops.push_back(std::move(loop));
    }
    std::sort(loops.begin(), loops.end(), [](const HandleLoop& x, const HandleLoop& y) { return x.length < y.length; });
    return loops;
}

}  // namespace dmw
