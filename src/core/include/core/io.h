#pragma once

#include <cstdint>
#include <span>
#include <string>
#include <string_view>

#include "core/mesh.h"

namespace dmw {

// Errors are returned, not thrown: see DECISIONS.md D15.
struct LoadResult {
    TriMesh mesh;
    std::string error;  // empty on success; otherwise "<location>: <reason>", e.g. "line 4: ..."
    bool ok() const { return error.empty(); }
};

// Parses Wavefront OBJ text (D14: text in, no file I/O).
// Supported: `v x y z [w]`, `f` with `i`, `i/t`, `i/t/n`, `i//n` corners,
// negative (relative) indices, polygons (fan-triangulated, D16), comments,
// CRLF line endings. Other directives (vt, vn, o, g, s, usemtl, mtllib) are ignored.
LoadResult parse_obj(std::string_view text);

// Parses STL, binary or ASCII (D14: bytes in, no file I/O).
// Format detection is by size, not by the "solid" prefix: a file is binary iff it is
// exactly 84 + 50 * N bytes, where N is the uint32 at offset 80 (D19).
// Output is an unwelded triangle soup: triangle i uses vertices 3i, 3i+1, 3i+2.
// Stored facet normals are ignored; orientation comes from vertex order (D19).
LoadResult parse_stl(std::span<const std::uint8_t> bytes);

}  // namespace dmw
