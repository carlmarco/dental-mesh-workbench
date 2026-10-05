#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

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

}  // namespace dmw
