#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/curvature.h"
#include "core/dec.h"
#include "core/generate.h"
#include "core/margin.h"

using namespace dmw;
using Catch::Matchers::WithinAbs;

namespace {

// Synthetic "crown in gingiva": a 12 x 12 mm patch (spacing 0.1 mm). Flat gingiva at z = 0; a crown of
// height 4 mm whose wall rises over [r0 - w, r0] with a smooth profile, so the wall's BASE (r = r0)
// is a concave crease: the true margin. A mid-height plane cut lands mid-wall, offset by ~w/2.
constexpr double kR0 = 3.0, kWall = 1.2, kHeight = 4.0;

TriMesh crown_patch() {
    auto m = make_grid(120, 120);
    for (auto& p : m.positions) {
        p.x = 12.0 * p.x - 6.0, p.y = 12.0 * p.y - 6.0;
        const double r = std::hypot(p.x, p.y);
        double t = std::clamp((kR0 - r) / kWall, 0.0, 1.0);  // 0 at the base, 1 at the top of the wall
        t = t * t * (3 - 2 * t);                               // smoothstep: zero slope at both ends
        p.z = kHeight * t + 0.3 * std::max(0.0, 1.0 - r * r / ((kR0 - kWall) * (kR0 - kWall))) * (r < kR0 - kWall);
    }
    return m;
}

// Ground truth: tooth = inside the wall's base circle.
std::vector<std::uint8_t> truth_labels(const HalfEdgeMesh& m) {
    std::vector<std::uint8_t> t(m.positions.size());
    for (std::size_t v = 0; v < t.size(); ++v) t[v] = std::hypot(m.positions[v].x, m.positions[v].y) < kR0 ? 1 : 0;
    return t;
}

HalfEdgeMesh he(const TriMesh& m) {
    auto r = build_halfedge(m);
    REQUIRE(r.ok());
    return std::move(r.mesh);
}

}  // namespace

TEST_CASE("margin: label boundary edges join differently labelled vertices", "[margin]") {
    const auto mesh = he(make_grid(4, 4));
    std::vector<std::uint8_t> lab(mesh.positions.size(), 0);
    lab[12] = 1;  // centre vertex of the 5 x 5 grid
    const auto edges = label_boundary_edges(mesh, lab);
    CHECK(edges.size() == one_ring(mesh, 12).size());  // one edge to each neighbour
}

TEST_CASE("margin: boundary metrics are zero for identical lines and measure a known offset", "[margin]") {
    std::vector<Vec3> a, b;
    for (int i = 0; i < 100; ++i) a.push_back({0.1 * i, 0, 0}), b.push_back({0.1 * i, 0.3, 0});
    const auto same = compare_boundaries(a, a);
    CHECK(same.assd == 0.0);
    CHECK(same.f1_025 == 1.0);
    const auto shifted = compare_boundaries(a, b);
    CHECK_THAT(shifted.assd, WithinAbs(0.3, 1e-12));
    CHECK_THAT(shifted.hd95, WithinAbs(0.3, 1e-12));
    CHECK(shifted.f1_025 == 0.0);
    CHECK(shifted.f1_050 == 1.0);
}

TEST_CASE("margin: area-weighted IoU", "[margin]") {
    const std::vector<std::uint8_t> p{1, 1, 0, 0}, t{1, 0, 1, 0};
    const std::vector<double> w{1.0, 2.0, 3.0, 4.0};
    CHECK_THAT(region_iou(p, t, w), WithinAbs(1.0 / 6.0, 1e-12));  // |{0}| / |{0,1,2}| by weight
}

TEST_CASE("margin: concavity-weighted Voronoi snaps to the crease; a plane cut lands mid-wall", "[margin]") {
    const auto mesh = he(crown_patch());
    const auto truth = truth_labels(mesh);
    const auto gt_line = edge_midpoints(mesh, label_boundary_edges(mesh, truth));
    const auto area = build_dec(mesh).star0;

    MarginParams plane;
    plane.method = MarginParams::Method::HeightPlane;
    plane.plane_quantile = 0.5;  // whatever height: any level cuts the wall, not its base
    MarginParams voronoi;        // defaults: valley-weighted, cusp seeds off below (flat-topped crown)
    voronoi.cusp_seeds = false;

    const auto rp = detect_margin(mesh, plane), rv = detect_margin(mesh, voronoi);
    const auto mp = compare_boundaries(edge_midpoints(mesh, rp.margin), gt_line);
    const auto mv = compare_boundaries(edge_midpoints(mesh, rv.margin), gt_line);
    CHECK(mv.assd < 0.15);           // within ~1.5 mesh spacings of the true crease
    CHECK(mv.assd < 0.5 * mp.assd);  // and much closer than the plane cut
    CHECK(region_iou(rv.tooth, truth, area) > 0.95);
}

TEST_CASE("margin: seed classifier hook - reject-all equals no cusp seeds, accept-all equals default", "[margin]") {
    // On the bumpy crown patch there are detected cusp seeds; a constant model fixes P(tooth).
    const auto mesh = he(crown_patch());
    const MarginInputs in = margin_inputs(mesh, true);
    REQUIRE_FALSE(in.cusp_tips.empty());  // otherwise both checks below would pass vacuously
    REQUIRE(in.seed_features.size() == in.cusp_tips.size() * kSeedFeatureCount);
    LogisticModel constant;
    constant.mean.assign(kSeedFeatureCount, 0.0);
    constant.scale.assign(kSeedFeatureCount, 1.0);
    constant.weights.assign(kSeedFeatureCount, 0.0);

    MarginParams base = margin_operating_point();
    MarginParams off = base;
    off.cusp_seeds = false;
    MarginParams reject = base, accept = base;
    constant.bias = -20.0;  // P ~ 0
    reject.seed_model = &constant;
    const auto rejected = margin_labels(mesh, in, {}, reject);
    CHECK(rejected == margin_labels(mesh, in, {}, off));
    LogisticModel always = constant;
    always.bias = 20.0;  // P ~ 1
    accept.seed_model = &always;
    CHECK(margin_labels(mesh, in, {}, accept) == margin_labels(mesh, in, {}, base));
}

TEST_CASE("margin: graph cut with zero smoothness equals the Voronoi labelling (unary = arrival order)", "[margin]") {
    const auto mesh = he(crown_patch());
    MarginParams voronoi;
    voronoi.cusp_seeds = false;
    MarginParams cut = voronoi;
    cut.method = MarginParams::Method::GraphCut;
    cut.cut_smoothness = 0.0;  // no boundary term: each vertex takes its cheaper unary, i.e. d_T < d_G
    const MarginInputs in = margin_inputs(mesh, false);
    const auto valley = valley_strength(build_dec(mesh), in.kmin, voronoi.curvature_scale);
    const auto a = margin_labels(mesh, in, valley, voronoi), b = margin_labels(mesh, in, valley, cut);
    std::size_t differ = 0;
    for (std::size_t v = 0; v < a.size(); ++v) differ += a[v] != b[v];
    CHECK(differ <= a.size() / 1000);  // only exact distance ties may break differently
}

TEST_CASE("margin: graph cut snaps to the crease like Voronoi", "[margin]") {
    const auto mesh = he(crown_patch());
    const auto truth = truth_labels(mesh);
    const auto gt_line = edge_midpoints(mesh, label_boundary_edges(mesh, truth));
    MarginParams cut;
    cut.method = MarginParams::Method::GraphCut;
    cut.cusp_seeds = false;
    const auto r = detect_margin(mesh, cut);
    CHECK(compare_boundaries(edge_midpoints(mesh, r.margin), gt_line).assd < 0.15);
    CHECK(region_iou(r.tooth, truth, build_dec(mesh).star0) > 0.95);
}

TEST_CASE("margin: a false tooth seed on flat gingiva floods a Voronoi region; the graph cut confines it", "[margin]") {
    // The D74/D77 failure: one tooth seed on gingiva. First arrival gives it every vertex it reaches
    // first (a fake tooth region); the cut pays mu per mm of boundary on flat gingiva, so the region
    // shrinks as mu grows (measured: 6.94 mm^2 Voronoi; 5.12, 1.70, 0.73, 0.056 at mu 0.3, 1, 3, 10).
    const auto mesh = he(crown_patch());
    const auto truth = truth_labels(mesh);
    const auto area = build_dec(mesh).star0;
    MarginInputs in = margin_inputs(mesh, false);
    std::uint32_t fake = 0;  // flat gingiva at (4.5, 0), between the crown and the patch edge
    for (std::uint32_t v = 0; v < mesh.positions.size(); ++v) {
        const Vec3& p = mesh.positions[v];
        if (std::hypot(p.x - 4.5, p.y) < std::hypot(mesh.positions[fake].x - 4.5, mesh.positions[fake].y)) fake = v;
    }
    in.cusp_tips = {fake};
    MarginParams voronoi;  // cusp seeds on (the injected one), no classifier
    // Gingiva seeds = the scan cut (patch border) only, as on a real arch: on this flat patch the 15%
    // height quantile would make every gingiva vertex a seed and leave nothing to flood.
    voronoi.gingiva_quantile = 0.0;
    MarginParams cut = voronoi;
    cut.method = MarginParams::Method::GraphCut;
    const auto valley = valley_strength(build_dec(mesh), in.kmin, voronoi.curvature_scale);
    auto false_tooth_area = [&](const std::vector<std::uint8_t>& lab) {
        double a = 0.0;
        for (std::size_t v = 0; v < lab.size(); ++v) a += (lab[v] && !truth[v]) ? area[v] : 0.0;
        return a;
    };
    const double av = false_tooth_area(margin_labels(mesh, in, valley, voronoi));
    std::vector<double> ac;
    for (double mu : {0.3, 1.0, 3.0, 10.0}) {
        cut.cut_smoothness = mu;
        ac.push_back(false_tooth_area(margin_labels(mesh, in, valley, cut)));
    }
    INFO("false tooth area: Voronoi " << av << " mm^2; cut " << ac[0] << ", " << ac[1] << ", " << ac[2] << ", " << ac[3]);
    CHECK(av > 0.5);  // mm^2: Voronoi grows a real fake region
    CHECK(ac[0] < av);
    CHECK(std::is_sorted(ac.rbegin(), ac.rend()));  // monotone in mu
    CHECK(ac[3] < 0.05 * av);                       // a speck around the hard seed
}

TEST_CASE("margin: cusp grouping splits crowns at the valley between them, keeps one crown's tips together", "[margin]") {
    // Two smooth domes side by side (z = max of the two): where they meet, the surface has a concave
    // valley (the max of two functions creases downward). Tips: each dome's apex, plus a second tip on dome A.
    auto m = make_grid(120, 60);
    auto dome = [](double x, double y, double cx) {
        const double r2 = (x - cx) * (x - cx) + y * y;
        return std::max(0.0, 4.0 - r2) / 4.0 * 3.0;  // radius 2 mm, height 3 mm
    };
    for (auto& p : m.positions) {
        p.x = 12.0 * p.x - 6.0, p.y = 6.0 * p.y - 3.0;
        p.z = std::max(dome(p.x, p.y, -1.6), dome(p.x, p.y, 1.6));  // centres 3.2 mm apart: domes overlap
    }
    const auto mesh = he(m);
    auto nearest = [&](double x, double y) {
        std::uint32_t best = 0;
        for (std::uint32_t v = 0; v < mesh.positions.size(); ++v)
            if (std::hypot(mesh.positions[v].x - x, mesh.positions[v].y - y) <
                std::hypot(mesh.positions[best].x - x, mesh.positions[best].y - y))
                best = v;
        return best;
    };
    const std::vector<std::uint32_t> tips{nearest(-1.6, 0.0), nearest(-1.2, 0.5), nearest(1.6, 0.0)};
    std::vector<std::uint8_t> region(mesh.positions.size(), 0);
    for (std::size_t v = 0; v < region.size(); ++v) region[v] = mesh.positions[v].z > 0.05;  // the domes
    std::vector<double> kmin(mesh.positions.size());
    const auto curv = compute_curvature(mesh);
    for (std::size_t v = 0; v < kmin.size(); ++v) kmin[v] = std::isfinite(curv.k2[v]) ? curv.k2[v] : 0.0;
    const auto g = group_cusps(mesh, tips, region, kmin, -2.0);
    REQUIRE(g.size() == 3);
    CHECK(g[0] >= 0);
    CHECK(g[0] == g[1]);  // same dome
    CHECK(g[0] != g[2]);  // across the valley
    // A threshold below the valley's concavity merges everything.
    const auto merged = group_cusps(mesh, tips, region, kmin, -1e9);
    CHECK((merged[0] == merged[1] && merged[1] == merged[2]));
}

TEST_CASE("margin: per-tooth cut with strip carving keeps gingiva between two touching crowns; the binary cut merges them", "[margin]") {
    // Two domes on a flat gingiva plane, overlapping slightly: where they meet, a valley runs up between
    // them (D84). Both crowns are "tooth" to a binary labelling, so it needs no boundary there.
    auto m = make_grid(120, 80);
    auto dome = [](double x, double y, double cx) {
        const double r2 = (x - cx) * (x - cx) + y * y;
        return std::max(0.0, 4.0 - r2) / 4.0 * 3.0;
    };
    for (auto& p : m.positions) {
        p.x = 12.0 * p.x - 6.0, p.y = 8.0 * p.y - 4.0;
        p.z = std::max(dome(p.x, p.y, -1.8), dome(p.x, p.y, 1.8));
    }
    const auto mesh = he(m);
    auto nearest = [&](double x, double y) {
        std::uint32_t best = 0;
        for (std::uint32_t v = 0; v < mesh.positions.size(); ++v)
            if (std::hypot(mesh.positions[v].x - x, mesh.positions[v].y - y) <
                std::hypot(mesh.positions[best].x - x, mesh.positions[best].y - y))
                best = v;
        return best;
    };
    MarginInputs in = margin_inputs(mesh, false);
    in.cusp_tips = {nearest(-1.8, 0.0), nearest(1.8, 0.0)};
    MarginParams binary = margin_operating_point();
    binary.gingiva_quantile = 0.0;  // gingiva seeds: the patch border (as on a real arch)
    MarginParams per_tooth = binary;
    per_tooth.method = MarginParams::Method::PerToothCut;
    per_tooth.strip_kmin = 0.0;  // the carving mechanism (off at the operating point, D85)
    const auto valley = valley_strength(build_dec(mesh), in.kmin, binary.curvature_scale);
    const auto a = margin_labels(mesh, in, valley, binary), b = margin_labels(mesh, in, valley, per_tooth);
    // Gingiva vertices high up in the valley between the crowns (|x| < 0.3 mm, z > 0.5 mm).
    auto strip = [&](const std::vector<std::uint8_t>& t) {
        std::size_t n = 0;
        for (std::size_t v = 0; v < t.size(); ++v)
            n += !t[v] && std::abs(mesh.positions[v].x) < 0.3 && mesh.positions[v].z > 0.5;
        return n;
    };
    CHECK(strip(a) == 0);  // binary: the crowns merge, no boundary between them
    CHECK(strip(b) > 0);   // per tooth: a strip of gingiva separates them
    // Elsewhere the two agree: the strip is thin.
    std::size_t differ = 0, crowns = 0;
    for (std::size_t v = 0; v < a.size(); ++v) differ += a[v] != b[v], crowns += a[v];
    CHECK(differ < crowns / 20);
}

TEST_CASE("margin: learned data term - a neutral model changes nothing, a confident one only grows the teeth", "[margin]") {
    const auto mesh = he(crown_patch());
    const MarginInputs in = margin_inputs(mesh, false);
    MarginParams base = margin_operating_point();
    base.method = MarginParams::Method::GraphCut;
    base.cusp_seeds = false;
    const auto valley = valley_strength(build_dec(mesh), in.kmin, base.curvature_scale);
    std::vector<double> feat;
    MarginParams with_out = base;
    with_out.vertex_features_out = &feat;
    const auto plain = margin_labels(mesh, in, valley, with_out);
    CHECK(feat.size() == mesh.positions.size() * kVertexFeatureCount);
    CHECK(vertex_feature_names().size() == kVertexFeatureCount);

    LogisticModel neutral;  // P = 0.5 everywhere: adds the same cost to both labels
    neutral.mean.assign(kVertexFeatureCount, 0.0);
    neutral.scale.assign(kVertexFeatureCount, 1.0);
    neutral.weights.assign(kVertexFeatureCount, 0.0);
    MarginParams n = base;
    n.vertex_model = &neutral;
    CHECK(margin_labels(mesh, in, valley, n) == plain);

    LogisticModel tooth = neutral;  // P ~ 1: tooth everywhere it is not hard-constrained
    tooth.bias = 10.0;
    MarginParams t = base;
    t.vertex_model = &tooth, t.vertex_weight = 1e3;
    const auto grown = margin_labels(mesh, in, valley, t);
    std::size_t lost = 0, gained = 0;
    for (std::size_t v = 0; v < plain.size(); ++v) lost += plain[v] && !grown[v], gained += !plain[v] && grown[v];
    CHECK(lost == 0);
    CHECK(gained > 0);
}

TEST_CASE("margin: the synthetic molar's tooth-gingiva boundary is found at its crease", "[margin]") {
    const auto mesh = he(make_synthetic_tooth());
    const MarginResult r = detect_margin(mesh, margin_operating_point());
    std::vector<std::uint8_t> truth(mesh.positions.size());
    for (std::size_t v = 0; v < truth.size(); ++v)
        truth[v] = std::hypot(mesh.positions[v].x, mesh.positions[v].y) < synthetic_tooth_radius ? 1 : 0;
    const auto m = compare_boundaries(edge_midpoints(mesh, r.margin), edge_midpoints(mesh, label_boundary_edges(mesh, truth)));
    CHECK(m.assd < 0.1);  // measured 0.021 mm
    CHECK(m.f1_050 > 0.95);
}
