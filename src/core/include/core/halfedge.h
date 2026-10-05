#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "core/mesh.h"

namespace dmw {

inline constexpr std::uint32_t kInvalid = std::numeric_limits<std::uint32_t>::max();

// Index-based half-edge structure for triangle meshes (D25).
//
// Numbering is implicit: face f owns half-edges 3f, 3f+1, 3f+2, and half-edge 3f+k
// runs from triangles[f][k] to triangles[f][(k+1) % 3]. So next/prev/face are pure
// arithmetic and are not stored; only origin and twin are.
//
// Invariants after a successful build:
//   next(next(next(h))) == h
//   twin[h] == kInvalid            <=> h is a boundary half-edge (no face on the other side)
//   twin[h] != kInvalid  =>  twin[twin[h]] == h, twin[h] != h,
//                            origin[twin[h]] == dest(h), dest(twin[h]) == origin[h]
//   vertex_halfedge[v] == kInvalid <=> v is referenced by no triangle (isolated)
//   otherwise origin[vertex_halfedge[v]] == v, and if v is on the boundary,
//   vertex_halfedge[v] is v's outgoing boundary half-edge (twin == kInvalid)
struct HalfEdgeMesh {
    std::vector<Vec3> positions;
    std::vector<std::uint32_t> origin;           // per half-edge: start vertex
    std::vector<std::uint32_t> twin;             // per half-edge: opposite half-edge or kInvalid
    std::vector<std::uint32_t> vertex_halfedge;  // per vertex: one outgoing half-edge or kInvalid
};

// Pure index arithmetic (valid for any h < 3 * face count).
std::uint32_t next(std::uint32_t h);
std::uint32_t prev(std::uint32_t h);
std::uint32_t face(std::uint32_t h);
std::uint32_t dest(const HalfEdgeMesh& m, std::uint32_t h);  // end vertex = origin of next

struct BuildResult {
    HalfEdgeMesh mesh;
    std::string error;  // empty on success
    bool ok() const { return error.empty(); }
};

// Builds the half-edge structure from an indexed triangle mesh (typically after welding).
// Fails with an error (D10, D25) on:
//   - a triangle index out of range
//   - a degenerate triangle (repeated vertex index)
//   - a non-manifold edge: the same directed edge (u,v) occurring twice, which means
//     3+ faces share the edge or two faces disagree on orientation
//   - a non-manifold vertex: a vertex whose faces form more than one fan ("bowtie")
BuildResult build_halfedge(const TriMesh& mesh);

// Neighbors of v in counter-clockwise order seen from outside (right-hand rule about
// the face normals). Interior vertex: starts at dest(vertex_halfedge[v]), k neighbors for
// k faces. Boundary vertex: starts at the far end of v's outgoing boundary half-edge and
// ends at the other boundary neighbor, k + 1 neighbors for k faces. Isolated: empty.
std::vector<std::uint32_t> one_ring(const HalfEdgeMesh& m, std::uint32_t v);

}  // namespace dmw
