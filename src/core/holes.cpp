#include "core/holes.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

#include "core/cholesky.h"
#include "core/sparse.h"
#include "detail/vec.h"

namespace dmw {
namespace {

using namespace detail;
using Tri = std::array<std::uint32_t, 3>;

Vec3 unit_normal(const Vec3& a, const Vec3& b, const Vec3& c) {
    const Vec3 n = cross(b - a, c - a);
    const double l = norm(n);
    return l > 0.0 ? (1.0 / l) * n : Vec3{0, 0, 0};
}
double angle_between(const Vec3& n1, const Vec3& n2) { return std::acos(std::clamp(dot(n1, n2), -1.0, 1.0)); }
double area(const Vec3& a, const Vec3& b, const Vec3& c) { return 0.5 * norm(cross(b - a, c - a)); }
std::uint64_t edge_key(std::uint32_t a, std::uint32_t b) { return (static_cast<std::uint64_t>(a) << 32) | b; }

// Lexicographic weight: (largest dihedral angle, total area); smaller is better.
struct Weight {
    double angle = 0.0, area = 0.0;
    bool operator<(const Weight& o) const { return angle < o.angle - 1e-12 || (angle <= o.angle + 1e-12 && area < o.area); }
};

}  // namespace

std::vector<std::vector<std::uint32_t>> boundary_loops(const TriMesh& m) {
    std::set<std::pair<std::uint32_t, std::uint32_t>> directed;
    for (const Tri& t : m.triangles)
        for (int k = 0; k < 3; ++k) directed.insert({t[static_cast<std::size_t>(k)], t[static_cast<std::size_t>((k + 1) % 3)]});
    std::multimap<std::uint32_t, std::uint32_t> out;  // boundary edges a -> b
    for (const auto& [a, b] : directed)
        if (!directed.count({b, a})) out.emplace(a, b);
    std::vector<std::vector<std::uint32_t>> loops;
    while (!out.empty()) {
        auto it = out.begin();
        const std::uint32_t start = it->first;
        std::vector<std::uint32_t> loop{start};
        std::uint32_t next = it->second;
        out.erase(it);
        while (next != start) {
            loop.push_back(next);
            const auto nx = out.find(next);
            if (nx == out.end()) break;  // broken boundary (non-manifold input): stop this loop
            const std::uint32_t after = nx->second;
            out.erase(nx);
            next = after;
        }
        if (next == start && loop.size() >= 3) loops.push_back(std::move(loop));
    }
    return loops;
}

std::pair<double, double> triangulation_weight(std::span<const Vec3> p, std::span<const Vec3> outside, std::span<const Tri> tris) {
    const auto n = static_cast<std::uint32_t>(p.size());
    // Faces: the patch triangles, then the rim faces (i, i+1, outside[i]) as indices n + i for their third vertex.
    std::vector<std::array<Vec3, 3>> faces;
    std::vector<bool> is_patch;
    std::map<std::pair<std::uint32_t, std::uint32_t>, std::vector<std::size_t>> by_edge;  // undirected edge -> faces
    auto add = [&](std::array<std::uint32_t, 3> ids, std::array<Vec3, 3> pos, bool patch) {
        const std::size_t f = faces.size();
        faces.push_back(pos), is_patch.push_back(patch);
        for (int k = 0; k < 3; ++k) {
            const std::uint32_t a = ids[static_cast<std::size_t>(k)], b = ids[static_cast<std::size_t>((k + 1) % 3)];
            by_edge[{std::min(a, b), std::max(a, b)}].push_back(f);
        }
    };
    double total = 0.0;
    for (const Tri& t : tris) {
        add(t, {p[t[0]], p[t[1]], p[t[2]]}, true);
        total += area(p[t[0]], p[t[1]], p[t[2]]);
    }
    for (std::uint32_t i = 0; i < n; ++i) add({i, (i + 1) % n, n + i}, {p[i], p[(i + 1) % n], outside[i]}, false);
    double worst = 0.0;
    for (const auto& [e, fs] : by_edge) {
        if (fs.size() != 2 || (!is_patch[fs[0]] && !is_patch[fs[1]])) continue;
        const auto& a = faces[fs[0]];
        const auto& b = faces[fs[1]];
        worst = std::max(worst, angle_between(unit_normal(a[0], a[1], a[2]), unit_normal(b[0], b[1], b[2])));
    }
    return {worst, total};
}

std::vector<Tri> triangulate_loop(std::span<const Vec3> p, std::span<const Vec3> outside, bool use_dihedral,
                                  const std::function<bool(std::uint32_t, std::uint32_t)>& chord_ok) {
    const std::size_t n = p.size();
    if (n < 3) return {};
    // W[i][j], best[i][j]: optimal triangulation of the polygon i, i+1, ..., j closed by the chord (i, j).
    std::vector<Weight> w(n * n);
    std::vector<std::uint32_t> best(n * n, 0);
    auto at = [n](std::size_t i, std::size_t j) { return i * n + j; };
    // Normal of the neighbour across edge (a, b) with a < b, from the side of the sub-polygon a..b: a rim face if the
    // edge is a loop edge, else the sub-solution's triangle (a, b, best[a][b]).
    auto neighbour = [&](std::size_t a, std::size_t b) {
        if (b == a + 1) return unit_normal(p[a], p[b], outside[a]);
        return unit_normal(p[a], p[b], p[best[at(a, b)]]);
    };
    for (std::size_t len = 2; len < n; ++len)
        for (std::size_t i = 0; i + len < n; ++i) {
            const std::size_t j = i + len;
            Weight bestw{std::numeric_limits<double>::infinity(), 0.0};
            for (std::size_t m = i + 1; m < j; ++m) {
                if (w[at(i, m)].angle == std::numeric_limits<double>::infinity() || w[at(m, j)].angle == std::numeric_limits<double>::infinity())
                    continue;  // a sub-polygon without a valid triangulation
                if (chord_ok) {
                    auto ok = [&](std::size_t a, std::size_t b) { return b == a + 1 || (a == 0 && b == n - 1) || chord_ok(std::uint32_t(a), std::uint32_t(b)); };
                    if (!ok(i, m) || !ok(m, j) || !ok(i, j)) continue;
                }
                const Vec3 nt = unit_normal(p[i], p[j], p[m]);  // triangle (i, j, m): consistent with the rim faces
                double ang = std::max({w[at(i, m)].angle, w[at(m, j)].angle, angle_between(nt, neighbour(i, m)), angle_between(nt, neighbour(m, j))});
                if (i == 0 && j == n - 1) ang = std::max(ang, angle_between(nt, unit_normal(p[n - 1], p[0], outside[n - 1])));
                if (!use_dihedral) ang = 0.0;
                const Weight cand{ang, w[at(i, m)].area + w[at(m, j)].area + area(p[i], p[j], p[m])};
                if (cand < bestw) bestw = cand, best[at(i, j)] = static_cast<std::uint32_t>(m);
            }
            w[at(i, j)] = bestw;
        }
    if (w[at(0, n - 1)].angle == std::numeric_limits<double>::infinity()) return {};  // every triangulation was forbidden
    std::vector<Tri> out;
    std::vector<std::pair<std::size_t, std::size_t>> todo{{0, n - 1}};
    while (!todo.empty()) {
        const auto [i, j] = todo.back();
        todo.pop_back();
        if (j < i + 2) continue;
        const std::size_t m = best[at(i, j)];
        out.push_back({static_cast<std::uint32_t>(i), static_cast<std::uint32_t>(j), static_cast<std::uint32_t>(m)});
        todo.push_back({i, m}), todo.push_back({m, j});
    }
    return out;
}

HoleFillResult fill_holes(const TriMesh& mesh, const HoleFillParams& prm) {
    HoleFillResult res;
    res.mesh = mesh;
    TriMesh& m = res.mesh;
    const auto loops = boundary_loops(mesh);
    res.loops = loops.size();
    // Third vertex of the existing face on each directed edge, and the mean incident edge length (density scale).
    std::unordered_map<std::uint64_t, std::uint32_t> third;
    std::vector<double> scale(mesh.positions.size(), 0.0);
    std::vector<std::uint32_t> count(mesh.positions.size(), 0);
    for (const Tri& t : mesh.triangles)
        for (int k = 0; k < 3; ++k) {
            const std::uint32_t a = t[static_cast<std::size_t>(k)], b = t[static_cast<std::size_t>((k + 1) % 3)], c = t[static_cast<std::size_t>((k + 2) % 3)];
            third[edge_key(a, b)] = c;
            const double l = norm(mesh.positions[b] - mesh.positions[a]);
            scale[a] += l, ++count[a];
        }
    for (std::size_t v = 0; v < scale.size(); ++v) scale[v] = count[v] ? scale[v] / count[v] : 0.0;
    // Every undirected edge of the mesh so far (input + finished patches): a patch must not duplicate any of them.
    std::unordered_set<std::uint64_t> edges;
    auto undirected = [](std::uint32_t a, std::uint32_t b) { return edge_key(std::min(a, b), std::max(a, b)); };
    for (const Tri& t : mesh.triangles)
        for (int k = 0; k < 3; ++k) edges.insert(undirected(t[static_cast<std::size_t>(k)], t[static_cast<std::size_t>((k + 1) % 3)]));

    for (const auto& loop : loops) {
        if (loop.size() > prm.max_loop) continue;
        std::set<std::uint32_t> uniq(loop.begin(), loop.end());
        if (uniq.size() != loop.size()) continue;  // a loop through a vertex twice (pinched boundary): skip
        std::vector<Vec3> pos, outside;
        for (std::size_t i = 0; i < loop.size(); ++i) {
            pos.push_back(mesh.positions[loop[i]]);
            outside.push_back(mesh.positions[third.at(edge_key(loop[i], loop[(i + 1) % loop.size()]))]);
        }
        std::vector<Tri> patch;
        const auto chord_ok = [&](std::uint32_t a, std::uint32_t b) { return !edges.count(undirected(loop[a], loop[b])); };
        const auto tris = triangulate_loop(pos, outside, prm.dihedral, chord_ok);
        if (tris.empty()) continue;  // every triangulation would duplicate an existing edge: leave this hole open
        for (const Tri& t : tris) patch.push_back({loop[t[0]], loop[t[1]], loop[t[2]]});

        if (prm.refine) {
            // Liepa's refinement: split at the centroid where it is far from every corner relative to the local scale,
            // then relax the patch by Delaunay flips; repeat until nothing splits.
            for (int round = 0; round < 30; ++round) {
                bool split = false;
                std::vector<Tri> next;
                for (const Tri& t : patch) {
                    const Vec3 c = (1.0 / 3.0) * (m.positions[t[0]] + m.positions[t[1]] + m.positions[t[2]]);
                    const double sc = (scale[t[0]] + scale[t[1]] + scale[t[2]]) / 3.0;
                    bool far = true;
                    for (auto v : t) {
                        const double d = std::sqrt(2.0) * norm(c - m.positions[v]);
                        far = far && d > sc && d > scale[v];
                    }
                    if (!far) {
                        next.push_back(t);
                        continue;
                    }
                    const auto id = static_cast<std::uint32_t>(m.positions.size());
                    m.positions.push_back(c), scale.push_back(sc), res.patch_vertices.push_back(id);
                    next.push_back({t[0], t[1], id}), next.push_back({t[1], t[2], id}), next.push_back({t[2], t[0], id});
                    split = true;
                }
                patch = std::move(next);
                // Relax: flip interior patch edges whose opposite angles sum to more than pi (not locally Delaunay).
                for (int pass = 0; pass < 50; ++pass) {
                    std::unordered_map<std::uint64_t, std::size_t> face_of;  // directed edge -> patch face
                    for (std::size_t f = 0; f < patch.size(); ++f)
                        for (int k = 0; k < 3; ++k) face_of[edge_key(patch[f][static_cast<std::size_t>(k)], patch[f][static_cast<std::size_t>((k + 1) % 3)])] = f;
                    bool flipped = false;
                    std::vector<bool> touched(patch.size(), false);
                    std::unordered_set<std::uint64_t> made;  // edges created by flips earlier in this pass (face_of is stale)
                    for (std::size_t f = 0; f < patch.size(); ++f) {
                        for (int k = 0; k < 3 && !touched[f]; ++k) {
                            const std::uint32_t a = patch[f][static_cast<std::size_t>(k)], b = patch[f][static_cast<std::size_t>((k + 1) % 3)],
                                                c = patch[f][static_cast<std::size_t>((k + 2) % 3)];
                            const auto g_it = face_of.find(edge_key(b, a));
                            if (g_it == face_of.end() || touched[g_it->second]) continue;  // rim edge or already changed
                            const std::size_t g = g_it->second;
                            std::uint32_t d = 0;
                            for (auto v : patch[g]) if (v != a && v != b) d = v;
                            const Vec3 &pa = m.positions[a], &pb = m.positions[b], &pc = m.positions[c], &pd = m.positions[d];
                            const double alpha = std::acos(std::clamp(dot((1.0 / norm(pa - pc)) * (pa - pc), (1.0 / norm(pb - pc)) * (pb - pc)), -1.0, 1.0));
                            const double beta = std::acos(std::clamp(dot((1.0 / norm(pa - pd)) * (pa - pd), (1.0 / norm(pb - pd)) * (pb - pd)), -1.0, 1.0));
                            if (alpha + beta <= std::acos(-1.0) + 1e-9) continue;
                            if (face_of.count(edge_key(c, d)) || face_of.count(edge_key(d, c)) || edges.count(undirected(c, d)) ||
                                made.count(undirected(c, d)))
                                continue;  // edge c-d exists already (in the patch, in the mesh, or from a flip this pass)
                            // The flip must keep both new triangles facing the same way as the old pair.
                            // Faces (a, b, c) and (b, a, d) bound the quad a -> d -> b -> c; the flip gives (d, b, c) and (c, a, d),
                            // and must keep both facing the same way as the old pair.
                            const Vec3 n_old = cross(pb - pa, pc - pa) + cross(pa - pb, pd - pb);
                            if (dot(cross(pb - pd, pc - pd), n_old) <= 0.0 || dot(cross(pa - pc, pd - pc), n_old) <= 0.0) continue;
                            patch[f] = {d, b, c}, patch[g] = {c, a, d};
                            made.insert(undirected(c, d));
                            touched[f] = touched[g] = true, flipped = true;
                        }
                    }
                    if (!flipped) break;
                }
                if (!split) break;
            }
        }
        for (const Tri& t : patch)
            for (int k = 0; k < 3; ++k) edges.insert(undirected(t[static_cast<std::size_t>(k)], t[static_cast<std::size_t>((k + 1) % 3)]));
        m.triangles.insert(m.triangles.end(), patch.begin(), patch.end());
        ++res.filled;
    }

    if (prm.fair && !res.patch_vertices.empty()) {
        // Thin-plate fairing of the new vertices: minimize |K x|^2_{D^-1} with every other vertex fixed.
        const std::size_t nv = m.positions.size();
        std::vector<std::set<std::uint32_t>> adj(nv);
        for (const Tri& t : m.triangles)
            for (int k = 0; k < 3; ++k) adj[t[static_cast<std::size_t>(k)]].insert(t[static_cast<std::size_t>((k + 1) % 3)]),
                                        adj[t[static_cast<std::size_t>((k + 1) % 3)]].insert(t[static_cast<std::size_t>(k)]);
        std::vector<std::int64_t> unknown(nv, -1);
        for (std::size_t q = 0; q < res.patch_vertices.size(); ++q) unknown[res.patch_vertices[q]] = static_cast<std::int64_t>(q);
        const auto nu = static_cast<std::uint32_t>(res.patch_vertices.size());
        // M = K D^-1 K, K = D - A: M_ij = sum_k K_ik K_kj / d_k. Rows of the unknowns only.
        std::vector<Triplet> trip;
        std::vector<std::array<double, 3>> rhs(nu, {0.0, 0.0, 0.0});
        auto kval = [&](std::uint32_t i, std::uint32_t j) { return i == j ? double(adj[i].size()) : (adj[i].count(j) ? -1.0 : 0.0); };
        for (std::uint32_t q = 0; q < nu; ++q) {
            const std::uint32_t i = res.patch_vertices[q];
            std::map<std::uint32_t, double> row;
            std::vector<std::uint32_t> ks(adj[i].begin(), adj[i].end());
            ks.push_back(i);
            for (std::uint32_t k : ks) {
                const double kik = kval(i, k) / double(adj[k].size());
                row[k] += kik * kval(k, k);
                for (std::uint32_t j : adj[k]) row[j] += kik * kval(k, j);
            }
            for (const auto& [j, v] : row) {
                if (unknown[j] >= 0) trip.push_back({q, static_cast<std::uint32_t>(unknown[j]), v});
                else rhs[q][0] -= v * m.positions[j].x, rhs[q][1] -= v * m.positions[j].y, rhs[q][2] -= v * m.positions[j].z;
            }
        }
        const EnvelopeLdlt ldlt(SparseMatrix::from_triplets(nu, nu, std::move(trip)));
        if (ldlt.ok()) {
            for (int c = 0; c < 3; ++c) {
                std::vector<double> b(nu);
                for (std::uint32_t q = 0; q < nu; ++q) b[q] = rhs[q][static_cast<std::size_t>(c)];
                const std::vector<double> x = ldlt.solve(b);
                for (std::uint32_t q = 0; q < nu; ++q) {
                    Vec3& p = m.positions[res.patch_vertices[q]];
                    (c == 0 ? p.x : c == 1 ? p.y : p.z) = x[q];
                }
            }
        }
    }
    return res;
}

}  // namespace dmw
