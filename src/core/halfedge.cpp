#include "core/halfedge.h"

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>


namespace dmw {

// h = 3f + k. Subtracting h % 3 gives the face's first half-edge 3f; adding
// (k+1) % 3 stays inside the face. (h + 1 would leak into the next face when k == 2.)
std::uint32_t next(std::uint32_t h) { return h - h % 3 + (h % 3 + 1) % 3; }
std::uint32_t prev(std::uint32_t h) { return h - h % 3 + (h % 3 + 2) % 3; }
std::uint32_t face(std::uint32_t h) { return h / 3; }
std::uint32_t dest(const HalfEdgeMesh& m, std::uint32_t h) { return m.origin[next(h)]; }

namespace {

// Directed edge (u, v) packed into one 64-bit key: u in the high word, v in the low.
// Distinct pairs give distinct keys, and (u,v) != (v,u), so direction is preserved.
std::uint64_t edge_key(std::uint32_t u, std::uint32_t v) {
    return (static_cast<std::uint64_t>(u) << 32) | v;
}

BuildResult fail(std::string why) { return {HalfEdgeMesh{}, std::move(why)}; }

}  // namespace

BuildResult build_halfedge(const TriMesh& mesh) {
    const std::size_t num_vertices = mesh.positions.size();
    const std::size_t num_faces = mesh.triangles.size();
    // Half-edge ids must fit in uint32 and never collide with the kInvalid sentinel.
    if (num_faces > (kInvalid - 1) / 3) return fail("too many faces for 32-bit half-edge ids");
    const auto num_halfedges = static_cast<std::uint32_t>(3 * num_faces);

    HalfEdgeMesh m;
    m.positions = mesh.positions;
    m.origin.resize(num_halfedges);
    m.twin.assign(num_halfedges, kInvalid);
    m.vertex_halfedge.assign(num_vertices, kInvalid);
    std::vector<std::uint32_t> out_degree(num_vertices, 0);  // |H_out(v)|, needed in step 6

    // Steps 1-2: validate each face, then record origins and tally out-degrees.
    for (std::size_t f = 0; f < num_faces; ++f) {
        const auto& tri = mesh.triangles[f];
        for (std::uint32_t v : tri) {
            if (v >= num_vertices) {
                return fail("face " + std::to_string(f) + ": vertex index " + std::to_string(v) +
                            " out of range");
            }
        }
        if (tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2]) {
            return fail("face " + std::to_string(f) + ": degenerate (repeated vertex)");
        }
        for (std::uint32_t k = 0; k < 3; ++k) {
            const auto h = static_cast<std::uint32_t>(3 * f) + k;
            m.origin[h] = tri[k];
            ++out_degree[tri[k]];
        }
    }

    // Step 3: index directed edges. A failed insert means delta(h) is not injective:
    // 3+ faces share the edge, or two faces traverse it in the same direction (D10).
    std::unordered_map<std::uint64_t, std::uint32_t> halfedge_of;
    halfedge_of.reserve(num_halfedges);
    for (std::uint32_t h = 0; h < num_halfedges; ++h) {
        const std::uint32_t u = m.origin[h], v = dest(m, h);
        const auto [it, inserted] = halfedge_of.try_emplace(edge_key(u, v), h);
        if (!inserted) {
            return fail("non-manifold or inconsistently oriented edge (" + std::to_string(u) + "," +
                        std::to_string(v) + ") in faces " + std::to_string(face(it->second)) +
                        " and " + std::to_string(face(h)));
        }
    }

    // Step 4: twin(h) = the half-edge traversing the reverse direction, if any.
    for (std::uint32_t h = 0; h < num_halfedges; ++h) {
        const auto it = halfedge_of.find(edge_key(dest(m, h), m.origin[h]));
        if (it != halfedge_of.end()) m.twin[h] = it->second;
    }

    // Step 5: one outgoing half-edge per vertex, preferring a boundary one (twin == kInvalid)
    // so that the CCW sweep starts at one end of the fan. Must come after step 4: boundary
    // status is only known once twins are.
    for (std::uint32_t h = 0; h < num_halfedges; ++h) {
        std::uint32_t& vh = m.vertex_halfedge[m.origin[h]];
        if (vh == kInvalid || (m.twin[h] == kInvalid && m.twin[vh] != kInvalid)) vh = h;
    }

    // Step 6: manifold-vertex check. rho(h) = twin(prev(h)) is injective (prev is a bijection,
    // twin an involution), so the orbit from vh either returns to vh or stops at a boundary;
    // it can't loop elsewhere. It covers every outgoing half-edge iff v's faces form a single
    // fan. If vh had a twin yet the walk stopped, next(twin(vh)) is an unvisited outgoing
    // half-edge, so the count comes up short and we reject: no special case needed.
    for (std::uint32_t v = 0; v < num_vertices; ++v) {
        const std::uint32_t start = m.vertex_halfedge[v];
        if (start == kInvalid) continue;  // isolated vertex: allowed
        std::uint32_t visited = 1;
        for (std::uint32_t h = start;;) {
            const std::uint32_t t = m.twin[prev(h)];
            if (t == kInvalid || t == start) break;
            h = t;
            ++visited;
        }
        if (visited != out_degree[v]) {
            return fail("vertex " + std::to_string(v) + ": non-manifold (faces form " +
                        "more than one fan)");
        }
    }

    return {std::move(m), {}};
}

// Step 7: walk the rho-orbit from vertex_halfedge[v], emitting dest(h). For a boundary vertex
// the walk stops at the incoming boundary half-edge prev(h), whose origin is the last neighbor.
std::vector<std::uint32_t> one_ring(const HalfEdgeMesh& m, std::uint32_t v) {
    std::vector<std::uint32_t> ring;
    const std::uint32_t start = m.vertex_halfedge[v];
    if (start == kInvalid) return ring;
    for (std::uint32_t h = start;;) {
        ring.push_back(dest(m, h));
        const std::uint32_t incoming = prev(h);
        const std::uint32_t t = m.twin[incoming];
        if (t == kInvalid) {
            ring.push_back(m.origin[incoming]);
            break;
        }
        if (t == start) break;
        h = t;
    }
    return ring;
}

}  // namespace dmw
