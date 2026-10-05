#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include "core/io.h"

namespace dmw {
namespace {

// Explicit little-endian writes, independent of the host's byte order (mirror of io_stl.cpp).
void put_u32(std::vector<std::uint8_t>& out, std::uint32_t v) {
    for (int i = 0; i < 4; ++i) out.push_back(static_cast<std::uint8_t>(v >> (8 * i)));
}
void put_f32(std::vector<std::uint8_t>& out, float f) {
    std::uint32_t bits = 0;
    std::memcpy(&bits, &f, sizeof bits);
    put_u32(out, bits);
}

}  // namespace

std::string write_obj(const TriMesh& mesh) {
    std::string out;
    out.reserve(mesh.positions.size() * 64 + mesh.triangles.size() * 24);
    char buf[128];
    for (const Vec3& p : mesh.positions) {
        // %.17g: 17 significant digits always round-trip an IEEE double through strtod.
        const int n = std::snprintf(buf, sizeof buf, "v %.17g %.17g %.17g\n", p.x, p.y, p.z);
        out.append(buf, static_cast<std::size_t>(n));
    }
    for (const auto& t : mesh.triangles) {
        const int n = std::snprintf(buf, sizeof buf, "f %u %u %u\n", t[0] + 1, t[1] + 1, t[2] + 1);
        out.append(buf, static_cast<std::size_t>(n));
    }
    return out;
}

std::vector<std::uint8_t> write_stl_binary(const TriMesh& mesh) {
    std::vector<std::uint8_t> out(80, 0);  // header: zeros (deliberately not "solid")
    out.reserve(84 + 50 * mesh.triangles.size());
    put_u32(out, static_cast<std::uint32_t>(mesh.triangles.size()));
    for (const auto& t : mesh.triangles) {
        for (int k = 0; k < 3; ++k) put_f32(out, 0.0f);  // normal: readers recompute from winding
        for (std::uint32_t v : t) {
            const Vec3& p = mesh.positions[v];
            put_f32(out, static_cast<float>(p.x));
            put_f32(out, static_cast<float>(p.y));
            put_f32(out, static_cast<float>(p.z));
        }
        out.push_back(0);  // attribute byte count
        out.push_back(0);
    }
    return out;
}

}  // namespace dmw
