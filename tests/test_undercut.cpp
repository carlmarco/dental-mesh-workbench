#include <cmath>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/bvh.h"
#include "core/generate.h"
#include "core/undercut.h"

using namespace dmw;
using Catch::Matchers::WithinAbs;

namespace {

Vec3 face_normal(const TriMesh& m, std::size_t f) {
    const Vec3 &a = m.positions[m.triangles[f][0]], &b = m.positions[m.triangles[f][1]], &c = m.positions[m.triangles[f][2]];
    const Vec3 u{b.x - a.x, b.y - a.y, b.z - a.z}, v{c.x - a.x, c.y - a.y, c.z - a.z};
    const Vec3 n{u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x};
    const double l = std::sqrt(n.x * n.x + n.y * n.y + n.z * n.z);
    return {n.x / l, n.y / l, n.z / l};
}

// Open frustum (side wall only) around the z axis: radius r0 at z = 0 shrinking to r1 at z = h (a prepared tooth
// with taper), then rotated by `tilt` radians about the x axis. Outward normals.
TriMesh frustum(double r0, double r1, double h, double tilt, int seg = 64, int rings = 12) {
    TriMesh m;
    for (int k = 0; k <= rings; ++k) {
        const double z = h * k / rings, r = r0 + (r1 - r0) * k / rings;
        for (int s = 0; s < seg; ++s) {
            const double a = 2.0 * std::acos(-1.0) * s / seg;
            const double x = r * std::cos(a), y = r * std::sin(a);
            m.positions.push_back({x, y * std::cos(tilt) - z * std::sin(tilt), y * std::sin(tilt) + z * std::cos(tilt)});
        }
    }
    for (int k = 0; k < rings; ++k)
        for (int s = 0; s < seg; ++s) {
            const auto a = static_cast<std::uint32_t>(k * seg + s), b = static_cast<std::uint32_t>(k * seg + (s + 1) % seg);
            const auto c = a + static_cast<std::uint32_t>(seg), d = b + static_cast<std::uint32_t>(seg);
            m.triangles.push_back({a, b, d});
            m.triangles.push_back({a, d, c});
        }
    return m;
}

}  // namespace

TEST_CASE("undercut: convex sphere - undercut is exactly the back-facing faces", "[undercut]") {
    const auto m = make_icosphere(3);
    const Bvh bvh(m.positions, m.triangles);
    const Vec3 d{0, 0, 1};
    const auto r = undercut_map(m, bvh, d);
    for (std::size_t f = 0; f < m.triangles.size(); ++f) {
        const Vec3 n = face_normal(m, f);
        CHECK(r.undercut[f] == (n.z < 0.0 ? 1 : 0));  // convex: no front face is hidden
    }
    CHECK_THAT(r.fraction(), WithinAbs(0.5, 0.03));
}

TEST_CASE("undercut: an overhang hides exactly its footprint", "[undercut]") {
    // Base plate 4 x 4 at z = 0, facing up; a 2 x 2 plate floating at z = 1 above its centre (also facing up).
    TriMesh m = make_grid(40, 40);
    for (auto& p : m.positions) p.x = 4.0 * p.x - 2.0, p.y = 4.0 * p.y - 2.0;
    const std::size_t base_faces = m.triangles.size();
    TriMesh top = make_grid(10, 10);
    for (auto& p : top.positions) p.x = 2.0 * p.x - 1.0, p.y = 2.0 * p.y - 1.0, p.z = 1.0;
    append(m, top);
    const Bvh bvh(m.positions, m.triangles);
    std::vector<std::uint8_t> region(m.triangles.size(), 0);
    for (std::size_t f = 0; f < base_faces; ++f) region[f] = 1;
    const auto r = undercut_map(m, bvh, {0, 0, 1}, region);
    std::size_t mismatches = 0;
    for (std::size_t f = 0; f < base_faces; ++f) {
        double cx = 0, cy = 0;
        for (auto v : m.triangles[f]) cx += m.positions[v].x / 3.0, cy += m.positions[v].y / 3.0;
        const bool under = std::abs(cx) < 1.0 && std::abs(cy) < 1.0;
        mismatches += r.undercut[f] != (under ? 1 : 0);
    }
    CHECK(mismatches == 0);
    CHECK_THAT(r.fraction(), WithinAbs(4.0 / 16.0, 1e-9));  // footprint 2 x 2 of a 4 x 4 plate
    for (std::size_t f = base_faces; f < m.triangles.size(); ++f) CHECK(r.undercut[f] == 0);  // outside the region
}

TEST_CASE("undercut: a tapered frustum has no undercut along its axis, and the search recovers that axis", "[undercut]") {
    const double taper = 6.0 * std::acos(-1.0) / 180.0, tilt = 20.0 * std::acos(-1.0) / 180.0;
    const double h = 5.0, r0 = 3.0, r1 = r0 - h * std::tan(taper);
    const auto m = frustum(r0, r1, h, tilt);
    const Bvh bvh(m.positions, m.triangles);
    const Vec3 axis{0, -std::sin(tilt), std::cos(tilt)};
    CHECK(undercut_map(m, bvh, axis).undercut_area == 0.0);
    CHECK(undercut_map(m, bvh, {0, 0, 1}).fraction() > 0.1);  // 20 degrees off a 6-degree taper: undercut
    const InsertionAxis best = best_insertion_axis(m, bvh, {}, {0, 0, 1}, 30.0);
    CHECK(best.result.undercut_area == 0.0);
    const double angle = std::acos(std::min(1.0, best.axis.x * axis.x + best.axis.y * axis.y + best.axis.z * axis.z)) * 180.0 / std::acos(-1.0);
    CHECK(angle <= 6.0);  // any axis within the taper angle is undercut-free; the search must land inside that cone
}

TEST_CASE("undercut: local refinement finds a narrow optimum the coarse spiral misses; the tilt limit is respected", "[undercut]") {
    const double deg = std::acos(-1.0) / 180.0;
    // 1-degree taper, tilted 13 degrees: undercut-free axes form a 1-degree cone, which 12 spiral samples miss.
    const auto narrow = frustum(3.0, 3.0 - 5.0 * std::tan(1.0 * deg), 5.0, 13.0 * deg, 96, 12);
    const Bvh nb(narrow.positions, narrow.triangles);
    const InsertionAxis coarse = best_insertion_axis(narrow, nb, {}, {0, 0, 1}, 30.0, 12);
    CHECK(coarse.result.undercut_area == 0.0);
    // Limit the search to 10 degrees around the hint while the zero-undercut cone sits at 14-26 degrees: the answer
    // must stay inside the limit (and so keep some undercut).
    const auto tilted = frustum(3.0, 3.0 - 5.0 * std::tan(6.0 * deg), 5.0, 20.0 * deg);
    const Bvh tb(tilted.positions, tilted.triangles);
    const InsertionAxis limited = best_insertion_axis(tilted, tb, {}, {0, 0, 1}, 10.0);
    CHECK(std::acos(std::min(1.0, limited.axis.z)) <= 10.0 * deg + 1e-9);
    CHECK(limited.result.undercut_area > 0.0);
}

TEST_CASE("undercut: threads give identical results; coarse-to-fine still finds an undercut-free axis", "[undercut]") {
    const double deg = std::acos(-1.0) / 180.0;
    const auto m = frustum(3.0, 3.0 - 5.0 * std::tan(3.0 * deg), 5.0, 15.0 * deg, 256, 24);  // 12,288 faces
    const Bvh bvh(m.positions, m.triangles);
    const auto serial = undercut_map(m, bvh, {0.1, 0.2, 1.0}, {}, 0.0, 1);
    const auto parallel = undercut_map(m, bvh, {0.1, 0.2, 1.0}, {}, 0.0, 4);
    CHECK(serial.undercut == parallel.undercut);
    // Same faces; the summed area differs only by floating-point summation order (chunks), deterministic per thread count.
    CHECK_THAT(parallel.undercut_area, WithinAbs(serial.undercut_area, 1e-12 * serial.region_area));
    // Withdrawing against the frustum's own axis, every face of the (upward-tapered) wall faces away: all must be
    // flagged, by every thread.
    const auto all_down = undercut_map(m, bvh, {0.0, std::sin(15.0 * deg), -std::cos(15.0 * deg)}, {}, 0.0, 4);  // -axis
    std::size_t flagged = 0;
    for (auto x : all_down.undercut) flagged += x;
    CHECK(flagged == m.triangles.size());
    // Coarse stage on ~1,000 of 12,288 faces, then a full-resolution polish: must still reach the 3-degree cone.
    const InsertionAxis best = best_insertion_axis(m, bvh, {}, {0, 0, 1}, 30.0, 60, 0.25, 0.0, 1000);
    CHECK(best.result.undercut_area == 0.0);
}

TEST_CASE("undercut: a misleadingly sparse coarse stage is corrected by the full-resolution polish", "[undercut]") {
    // 1-degree taper tilted 13 degrees, coarse stage on only ~24 faces: the coarse optimum (nearest the hint among
    // its apparently undercut-free directions) need not be undercut-free at full resolution.
    const double deg = std::acos(-1.0) / 180.0;
    const auto m = frustum(3.0, 3.0 - 5.0 * std::tan(1.0 * deg), 5.0, 13.0 * deg, 192, 16);
    const Bvh bvh(m.positions, m.triangles);
    const InsertionAxis best = best_insertion_axis(m, bvh, {}, {0, 0, 1}, 30.0, 60, 0.25, 0.0, 24);
    CHECK(best.result.undercut_area == 0.0);
}
