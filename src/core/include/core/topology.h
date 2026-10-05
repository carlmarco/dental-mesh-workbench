#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include "core/halfedge.h"
#include "core/mesh.h"

namespace dmw {

// Undirected edge with v0 < v1.
struct Edge {
    std::uint32_t v0;
    std::uint32_t v1;
};

enum class EdgeKind : std::uint8_t {
    Boundary,     // 1 incident face
    Manifold,     // 2 faces traversing it in opposite directions (consistent winding)
    Misoriented,  // 2 faces traversing it in the same direction (winding flips here)
    NonManifold,  // 3 or more faces
};

// Topology of one connected component (vertex-connectivity, D28).
struct ComponentTopology {
    std::uint32_t num_vertices = 0;
    std::uint32_t num_edges = 0;
    std::uint32_t num_faces = 0;
    std::int64_t euler_characteristic = 0;  // V - E + F
    std::uint32_t boundary_loops = 0;       // b; meaningful only when manifold
    bool manifold = true;                   // no non-manifold edges or vertices
    bool orientable = false;                // a consistent winding exists; false if !manifold
    bool consistently_oriented = true;      // no misoriented edges (winding already consistent)
    // Set iff manifold && orientable (D29). g = (2 - chi - b) / 2;
    // betti = (b0, b1, b2) = (1, b > 0 ? 2g + b - 1 : 2g, b == 0 ? 1 : 0).
    std::optional<std::uint32_t> genus;
    std::optional<std::array<std::uint32_t, 3>> betti;
};

// Combinatorial diagnostics of a raw (possibly broken) indexed mesh. Works on any TriMesh,
// including those build_halfedge rejects (D10). Index lists are sorted ascending.
struct TopologyReport {
    std::vector<Edge> edges;          // unique undirected edges of valid faces
    std::vector<EdgeKind> edge_kind;  // parallel to edges

    std::vector<std::uint32_t> invalid_faces;         // out-of-range or repeated vertex index
    std::vector<std::uint32_t> duplicate_faces;       // same vertex set as an earlier face
    std::vector<std::uint32_t> isolated_vertices;     // in no valid face
    std::vector<std::uint32_t> nonmanifold_vertices;  // faces around v form more than one fan

    std::vector<std::uint32_t> vertex_component;  // per vertex; kInvalid if isolated
    std::vector<std::uint32_t> face_component;    // per face; kInvalid if invalid/duplicate
    std::vector<ComponentTopology> components;    // numbered by lowest vertex index

    std::size_t count(EdgeKind kind) const;
};

// Invalid and duplicate faces are reported and excluded from everything else (D31), so the
// counts describe the underlying simplicial complex. O(F log F) worst case, O(F) typical.
TopologyReport analyze_topology(const TriMesh& mesh);

// The mesh without the report's invalid and duplicate faces (D63). Vertices are kept unchanged,
// so per-vertex results still line up with the input. Real scans contain collapsed triangles
// (all three corners the same vertex); excluding them lets the half-edge build succeed when the
// rest of the surface is manifold. Non-manifold edges and vertices are NOT touched.
TriMesh without_excluded_faces(const TriMesh& mesh, const TopologyReport& report);

// Analysis view for the half-edge-based algorithms on real scans (D69). Repeatedly excludes the
// report's invalid and duplicate faces plus every face touching a non-manifold or misoriented edge
// or a non-manifold vertex, until the remainder is a manifold the half-edge build accepts (or
// `max_passes` is reached). Vertex indexing is unchanged; the input is not modified. Typical real
// scan: ~0.05% of faces excluded instead of rejecting the whole scan.
struct AnalysisMesh {
    TriMesh mesh;
    HalfEdgeMesh halfedge;           // valid iff manifold: built once here, so callers don't rebuild
    std::size_t excluded_faces = 0;  // relative to the input
    int passes = 0;
    bool manifold = false;           // build_halfedge(mesh) succeeds
};
AnalysisMesh manifold_analysis_mesh(const TriMesh& mesh, int max_passes = 4);

}  // namespace dmw
