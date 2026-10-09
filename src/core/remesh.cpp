#include "core/remesh.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "core/bvh.h"
#include "detail/vec.h"

namespace dmw {
namespace {

using namespace detail;
using Tri = std::array<std::uint32_t, 3>;

// Editable triangle mesh: faces (dead ones kept as holes in the array) and per-vertex lists of live faces.
struct Editable {
    std::vector<Vec3> p;
    std::vector<Tri> f;
    std::vector<bool> dead;
    std::vector<std::vector<std::uint32_t>> vf;

    explicit Editable(const TriMesh& m) : p(m.positions), f(m.triangles), dead(m.triangles.size(), false), vf(m.positions.size()) {
        for (std::uint32_t i = 0; i < f.size(); ++i)
            for (auto v : f[i]) vf[v].push_back(i);
    }
    static bool has(const Tri& t, std::uint32_t v) { return t[0] == v || t[1] == v || t[2] == v; }
    std::vector<std::uint32_t> edge_faces(std::uint32_t a, std::uint32_t b) const {
        std::vector<std::uint32_t> out;
        for (auto fi : vf[a])
            if (has(f[fi], b)) out.push_back(fi);
        return out;
    }
    static std::uint32_t opposite(const Tri& t, std::uint32_t a, std::uint32_t b) {
        for (auto v : t)
            if (v != a && v != b) return v;
        return t[0];
    }
    std::set<std::uint32_t> neighbours(std::uint32_t v) const {
        std::set<std::uint32_t> s;
        for (auto fi : vf[v])
            for (auto w : f[fi])
                if (w != v) s.insert(w);
        return s;
    }
    bool boundary_edge(std::uint32_t a, std::uint32_t b) const { return edge_faces(a, b).size() == 1; }
    bool boundary_vertex(std::uint32_t v) const {
        for (auto w : neighbours(v))
            if (boundary_edge(v, w)) return true;
        return false;
    }
    std::vector<std::pair<std::uint32_t, std::uint32_t>> edges() const {
        std::set<std::pair<std::uint32_t, std::uint32_t>> s;
        for (std::uint32_t i = 0; i < f.size(); ++i)
            if (!dead[i])
                for (int k = 0; k < 3; ++k) {
                    const auto a = f[i][static_cast<std::size_t>(k)], b = f[i][static_cast<std::size_t>((k + 1) % 3)];
                    s.insert({std::min(a, b), std::max(a, b)});
                }
        return {s.begin(), s.end()};
    }
    void remove_face_from(std::uint32_t v, std::uint32_t fi) { std::erase(vf[v], fi); }
    Vec3 normal(const Tri& t) const { return cross(p[t[1]] - p[t[0]], p[t[2]] - p[t[0]]); }
};

double len(const Editable& m, std::uint32_t a, std::uint32_t b) { return norm(m.p[a] - m.p[b]); }

// Split edge (a, b) at its midpoint: each incident face (u, v, c) with {u, v} = {a, b} becomes (u, m, c), (m, v, c).
void split(Editable& m, std::uint32_t a, std::uint32_t b) {
    const auto mid = static_cast<std::uint32_t>(m.p.size());
    m.p.push_back(0.5 * (m.p[a] + m.p[b]));
    m.vf.emplace_back();
    for (std::uint32_t fi : m.edge_faces(a, b)) {
        Tri t = m.f[fi];
        while (!((t[0] == a && t[1] == b) || (t[0] == b && t[1] == a))) t = {t[1], t[2], t[0]};  // rotate: edge first
        const std::uint32_t u = t[0], v = t[1], c = t[2];
        m.f[fi] = {u, mid, c};
        const auto nf = static_cast<std::uint32_t>(m.f.size());
        m.f.push_back({mid, v, c});
        m.dead.push_back(false);
        m.remove_face_from(v, fi);
        m.vf[v].push_back(nf);
        m.vf[c].push_back(nf);
        m.vf[mid].push_back(fi), m.vf[mid].push_back(nf);
    }
}

// Collapse b into a (a moves to `x`), if safe. Returns true if collapsed.
bool collapse(Editable& m, std::uint32_t a, std::uint32_t b, const Vec3& x, double max_len) {
    const auto ef = m.edge_faces(a, b);
    if (ef.empty() || ef.size() > 2) return false;
    // Link condition: the common neighbours of a and b are exactly the opposite vertices of the edge's faces.
    std::set<std::uint32_t> opp;
    for (auto fi : ef) opp.insert(Editable::opposite(m.f[fi], a, b));
    const auto na = m.neighbours(a), nb = m.neighbours(b);
    std::set<std::uint32_t> common;
    std::set_intersection(na.begin(), na.end(), nb.begin(), nb.end(), std::inserter(common, common.begin()));
    if (common != opp) return false;
    for (auto o : opp)
        if (m.neighbours(o).size() <= 3) return false;  // would drop below valence 3
    for (auto w : nb)
        if (w != a && norm(m.p[w] - x) > max_len) return false;  // no new long edge
    for (auto w : na)
        if (w != b && norm(m.p[w] - x) > max_len) return false;
    // No face may flip (or degenerate) when a and b move to x.
    for (std::uint32_t v : {a, b})
        for (auto fi : m.vf[v]) {
            if (std::find(ef.begin(), ef.end(), fi) != ef.end()) continue;
            Tri t = m.f[fi];
            const Vec3 before = m.normal(t);
            for (auto& w : t)
                if (w == a || w == b) w = a;
            std::vector<Vec3> q{m.p[t[0]], m.p[t[1]], m.p[t[2]]};
            for (int k = 0; k < 3; ++k)
                if (t[static_cast<std::size_t>(k)] == a) q[static_cast<std::size_t>(k)] = x;
            const Vec3 after = cross(q[1] - q[0], q[2] - q[0]);
            if (dot(before, after) <= 0.0 || norm(after) < 1e-12 * norm(before)) return false;
        }
    // Apply.
    for (auto fi : ef) {
        m.dead[fi] = true;
        for (auto v : m.f[fi]) m.remove_face_from(v, fi);
    }
    for (auto fi : m.vf[b]) {
        for (auto& w : m.f[fi])
            if (w == b) w = a;
        m.vf[a].push_back(fi);
    }
    m.vf[b].clear();
    m.p[a] = x;
    return true;
}

}  // namespace

TriMesh remesh_isotropic(const TriMesh& mesh, const RemeshParams& prm) {
    Editable m(mesh);
    double target = prm.target_edge;
    if (target <= 0.0) {
        double s = 0.0;
        const auto es = m.edges();
        for (const auto& [a, b] : es) s += len(m, a, b);
        target = es.empty() ? 1.0 : s / double(es.size());
    }
    const double hi = 4.0 / 3.0 * target, lo = 4.0 / 5.0 * target;
    const Bvh reference(mesh.positions, mesh.triangles);
    for (int it = 0; it < prm.iterations; ++it) {
        // 1. Split long edges (repeat: a split can leave halves that are still long).
        for (int pass = 0; pass < 8; ++pass) {
            bool any = false;
            for (const auto& [a, b] : m.edges())
                if (len(m, a, b) > hi) split(m, a, b), any = true;
            if (!any) break;
        }
        // 2. Collapse short edges. Boundary vertices are never moved inward: a boundary vertex keeps its position; an
        // interior edge between two boundary vertices is not collapsed (it would pinch the surface).
        for (const auto& [a0, b0] : m.edges()) {
            std::uint32_t a = a0, b = b0;
            if (m.vf[a].empty() || m.vf[b].empty() || len(m, a, b) >= lo) continue;
            const bool ba = m.boundary_vertex(a), bb = m.boundary_vertex(b), be = m.boundary_edge(a, b);
            if (ba && bb && !be) continue;
            if (bb && !ba) std::swap(a, b);  // keep the boundary vertex: collapse the interior one into it
            const Vec3 x = (ba || bb) ? (be ? 0.5 * (m.p[a] + m.p[b]) : m.p[a]) : 0.5 * (m.p[a] + m.p[b]);
            collapse(m, a, b, x, hi);
        }
        // 3. Flip edges towards valence 6 (4 on the boundary).
        auto valence = [&](std::uint32_t v) { return static_cast<int>(m.neighbours(v).size()); };
        auto target_valence = [&](std::uint32_t v) { return m.boundary_vertex(v) ? 4 : 6; };
        for (const auto& [a, b] : m.edges()) {
            const auto ef = m.edge_faces(a, b);
            if (ef.size() != 2) continue;
            // Orient: face f0 contains a -> b, f1 contains b -> a.
            std::uint32_t f0 = ef[0], f1 = ef[1];
            auto has_dir = [&](std::uint32_t fi, std::uint32_t u, std::uint32_t v) {
                const Tri& t = m.f[fi];
                for (int k = 0; k < 3; ++k)
                    if (t[static_cast<std::size_t>(k)] == u && t[static_cast<std::size_t>((k + 1) % 3)] == v) return true;
                return false;
            };
            if (!has_dir(f0, a, b)) std::swap(f0, f1);
            if (!has_dir(f0, a, b) || !has_dir(f1, b, a)) continue;  // inconsistent orientation: leave it
            const std::uint32_t c = Editable::opposite(m.f[f0], a, b), d = Editable::opposite(m.f[f1], a, b);
            if (c == d || m.neighbours(c).count(d)) continue;  // edge c-d exists already
            const int va = valence(a), vb = valence(b), vc = valence(c), vd = valence(d);
            const int before = std::abs(va - target_valence(a)) + std::abs(vb - target_valence(b)) + std::abs(vc - target_valence(c)) +
                               std::abs(vd - target_valence(d));
            const int after = std::abs(va - 1 - target_valence(a)) + std::abs(vb - 1 - target_valence(b)) +
                              std::abs(vc + 1 - target_valence(c)) + std::abs(vd + 1 - target_valence(d));
            if (after >= before || va <= 3 || vb <= 3) continue;
            // New faces (d, b, c) and (c, a, d) must keep the orientation of the old pair.
            const Vec3 n_old = m.normal(m.f[f0]) + m.normal(m.f[f1]);
            const Tri t0{d, b, c}, t1{c, a, d};
            if (dot(m.normal(t0), n_old) <= 0.0 || dot(m.normal(t1), n_old) <= 0.0) continue;
            for (auto v : m.f[f0]) m.remove_face_from(v, f0);
            for (auto v : m.f[f1]) m.remove_face_from(v, f1);
            m.f[f0] = t0, m.f[f1] = t1;
            for (auto v : t0) m.vf[v].push_back(f0);
            for (auto v : t1) m.vf[v].push_back(f1);
        }
        // 4. Tangential relaxation of interior vertices, then projection onto the original surface.
        std::vector<Vec3> next = m.p;
        for (std::uint32_t v = 0; v < m.p.size(); ++v) {
            if (m.vf[v].empty() || m.boundary_vertex(v)) continue;
            Vec3 c{}, n{};
            const auto nb = m.neighbours(v);
            for (auto w : nb) c += m.p[w];
            c = (1.0 / double(nb.size())) * c;
            for (auto fi : m.vf[v]) n += m.normal(m.f[fi]);
            const double nl = norm(n);
            if (nl == 0.0) continue;
            n = (1.0 / nl) * n;
            const Vec3 step = c - m.p[v];
            next[v] = m.p[v] + prm.relaxation * (step - dot(step, n) * n);
            next[v] = reference.closest(next[v]).point;
        }
        m.p = std::move(next);
    }
    // Compact: drop dead faces and unused vertices.
    TriMesh out;
    std::vector<std::uint32_t> remap(m.p.size(), 0xFFFFFFFFu);
    for (std::uint32_t i = 0; i < m.f.size(); ++i) {
        if (m.dead[i]) continue;
        Tri t = m.f[i];
        for (auto& v : t) {
            if (remap[v] == 0xFFFFFFFFu) remap[v] = static_cast<std::uint32_t>(out.positions.size()), out.positions.push_back(m.p[v]);
            v = remap[v];
        }
        out.triangles.push_back(t);
    }
    return out;
}

}  // namespace dmw
