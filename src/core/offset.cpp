#include "core/offset.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <utility>

#include "core/isosurface.h"
#include "detail/vec.h"

namespace dmw {

using namespace detail;

TriMesh offset_shell_sdf(const TriMesh& surface, double wall, double grid_h, std::size_t max_nodes) {
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (const auto& t : surface.triangles)
        for (auto v : t) {
            const Vec3& p = surface.positions[v];
            lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)}, hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
        }
    SignedHeatParams prm;
    const Vec3 span = hi - lo;
    const double longest = std::max({span.x, span.y, span.z});
    prm.padding = std::max(2.0 * wall, 0.15 * longest);
    prm.h = grid_h > 0.0 ? grid_h : (longest + 2.0 * prm.padding) / double(std::max<std::size_t>(max_nodes, 8) - 1);
    const Grid3 phi = signed_heat_distance(surface, prm);
    TriMesh shell = extract_isosurface(phi, 0.0);
    const TriMesh inner = extract_isosurface(phi, -wall);
    const auto offset = static_cast<std::uint32_t>(shell.positions.size());
    shell.positions.insert(shell.positions.end(), inner.positions.begin(), inner.positions.end());
    for (auto t : inner.triangles) shell.triangles.push_back({t[0] + offset, t[2] + offset, t[1] + offset});  // flipped: faces the hollow
    return shell;
}

TriMesh offset_shell_naive(const TriMesh& surface, double wall) {
    const auto n = static_cast<std::uint32_t>(surface.positions.size());
    TriMesh s = surface;
    std::vector<Vec3> normal(n, Vec3{});
    std::map<std::pair<std::uint32_t, std::uint32_t>, int> edge_use;  // directed edges of the input faces
    for (const auto& t : surface.triangles) {
        const Vec3 nn = cross(surface.positions[t[1]] - surface.positions[t[0]], surface.positions[t[2]] - surface.positions[t[0]]);
        for (auto v : t) normal[v] += nn;
        for (auto [p, q] : {std::pair{t[0], t[1]}, std::pair{t[1], t[2]}, std::pair{t[2], t[0]}}) ++edge_use[{p, q}];
    }
    for (std::uint32_t v = 0; v < n; ++v) {
        const double l = norm(normal[v]);
        s.positions.push_back(l > 0.0 ? surface.positions[v] - (wall / l) * normal[v] : surface.positions[v]);
    }
    for (const auto& t : surface.triangles) s.triangles.push_back({t[0] + n, t[2] + n, t[1] + n});  // inner, flipped
    for (const auto& [e, count] : edge_use) {
        if (edge_use.count({e.second, e.first})) continue;  // interior edge
        const std::uint32_t a = e.first, b = e.second;      // boundary edge a -> b of an input face
        s.triangles.push_back({b, a, a + n});
        s.triangles.push_back({b, a + n, b + n});
    }
    return s;
}

}  // namespace dmw
