#pragma once

#include <array>
#include <cstdint>
#include <limits>
#include <vector>

namespace dmw {

// "No such element" for uint32 indices (no twin, isolated vertex, excluded face, ...).
inline constexpr std::uint32_t kInvalid = std::numeric_limits<std::uint32_t>::max();

// double, not float: see DECISIONS.md D12.
struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

// Indexed triangle mesh ("face-vertex" form): the interchange format between
// loaders, diagnostics, and the half-edge builder.
//   - vertex i is positions[i]
//   - each triangle lists 3 indices into positions; winding order is preserved
//     from the file (counter-clockwise seen from outside = outward normal by the
//     right-hand rule, for well-formed input)
// uint32_t indices: see DECISIONS.md D13.
struct TriMesh {
    std::vector<Vec3> positions;
    std::vector<std::array<std::uint32_t, 3>> triangles;
};

}  // namespace dmw
