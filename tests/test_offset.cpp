#include <algorithm>
#include <cmath>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/bvh.h"
#include "core/generate.h"
#include "core/halfedge.h"
#include "core/offset.h"
#include "core/thickness.h"
#include "core/topology.h"

using namespace dmw;
using Catch::Matchers::WithinAbs;

TEST_CASE("offset: SDF shell of a sphere - two closed walls at radius 1 and 1 - w, thickness w", "[offset]") {
    const double w = 0.25;
    const TriMesh shell = offset_shell_sdf(make_icosphere(4), w, 0.05);
    const auto topo = analyze_topology(shell);
    REQUIRE(topo.components.size() == 2);
    for (const auto& c : topo.components) {
        CHECK(c.manifold);
        CHECK(c.boundary_loops == 0);
    }
    CHECK(build_halfedge(shell).ok());  // consistently oriented: outer outward, inner towards the hollow
    const Bvh bvh(shell.positions, shell.triangles);
    const auto t = wall_thickness(shell, bvh);
    std::vector<double> along;
    for (double x : t.along_normal)
        if (std::isfinite(x)) along.push_back(x);
    REQUIRE(along.size() > shell.positions.size() / 2);
    std::sort(along.begin(), along.end());
    CHECK_THAT(along[along.size() / 2], WithinAbs(w, 0.01));
}

TEST_CASE("offset: naive shell of an open hemisphere is closed and consistently oriented", "[offset]") {
    const auto sphere = make_icosphere(3);
    TriMesh hemi;
    hemi.positions = sphere.positions;
    for (const auto& t : sphere.triangles) {
        bool keep = true;
        for (auto v : t) keep = keep && sphere.positions[v].z >= -1e-12;
        if (keep) hemi.triangles.push_back(t);
    }
    // Drop vertices no face uses (the lower half), so the shell has no isolated vertices.
    std::vector<std::uint32_t> remap(hemi.positions.size(), 0xFFFFFFFFu);
    TriMesh compact;
    for (auto& t : hemi.triangles) {
        for (auto& v : t) {
            if (remap[v] == 0xFFFFFFFFu) remap[v] = static_cast<std::uint32_t>(compact.positions.size()), compact.positions.push_back(hemi.positions[v]);
            v = remap[v];
        }
        compact.triangles.push_back(t);
    }
    const TriMesh shell = offset_shell_naive(compact, 0.1);
    const auto topo = analyze_topology(shell);
    REQUIRE(topo.components.size() == 1);
    CHECK(topo.components[0].boundary_loops == 0);
    CHECK(topo.components[0].genus.value_or(99u) == 0);
    CHECK(build_halfedge(shell).ok());
}

TEST_CASE("offset: synthetic molar - the naive shell folds under the fissure, the SDF shell keeps the wall", "[offset]") {
    const TriMesh tooth = make_synthetic_tooth(100);
    auto occlusal_min = [](const TriMesh& shell) {  // thinnest wall over the occlusal table (r < 3 mm, z > 3 mm)
        const Bvh bvh(shell.positions, shell.triangles);
        const auto t = wall_thickness(shell, bvh);
        double mn = 1e300;
        for (std::size_t v = 0; v < shell.positions.size(); ++v) {
            const Vec3& p = shell.positions[v];
            if (std::hypot(p.x, p.y) < 3.0 && p.z > 3.0 && std::isfinite(t.cone_min[v])) mn = std::min(mn, t.cone_min[v]);
        }
        return mn;
    };
    const double naive = occlusal_min(offset_shell_naive(tooth, 1.0)), sdf = occlusal_min(offset_shell_sdf(tooth, 1.0, 0.0, 64));
    INFO("thinnest occlusal wall for a 1 mm shell: naive " << naive << " mm, SDF " << sdf << " mm");
    CHECK(naive < 0.5);  // measured ~0: the inner copy folds through itself under the 0.3 mm fissure
    CHECK(sdf > 0.85);   // measured ~0.91 at 72 nodes: the nominal wall within grid error
}
