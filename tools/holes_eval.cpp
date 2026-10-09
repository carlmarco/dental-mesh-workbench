// Hole filling on real scans (M10f, D98). Two measurements per scan:
//   1. the scan's own small holes (boundary loops up to max_loop edges; the long open cut is left alone): how many,
//      whether all close, whether the result is still a consistently oriented manifold;
//   2. accuracy against ground truth: a disk of radius `r` mm is punched out of a tooth crown (FDI labels pick a
//      crown vertex away from boundaries), filled, and every patch vertex is compared with the removed surface.
//   holes_eval <scan-dir> [stride] [radius_mm]
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "core/halfedge.h"
#include "core/holes.h"
#include "core/io.h"
#include "core/margin.h"
#include "dataset.h"

using namespace dmw;
using Clock = std::chrono::steady_clock;

namespace {
double median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}
}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: holes_eval <scan-dir> [stride] [radius_mm]\n");
        return 2;
    }
    const std::size_t stride = argc > 2 ? std::stoul(argv[2]) : 15;
    const double radius = argc > 3 ? std::stod(argv[3]) : 1.0;
    const auto objs = dataset::index_files({argv[1]}, ".obj");
    const auto labels = dataset::index_files({argv[1]}, ".json");
    std::mt19937 rng(98);
    std::size_t sampled = 0;
    std::size_t index = 0, scans = 0, small_loops = 0, filled = 0, still_valid = 0, punched = 0, punched_closed = 0;
    std::vector<double> loops_per_scan, fill_ms, mean_err, max_err;
    for (const auto& [stem, obj] : objs) {
        if (index++ % stride != 0) continue;
        ++sampled;
        const auto lab = labels.find(stem);
        if (lab == labels.end()) continue;
        const JsonResult j = parse_json(dataset::read_text(lab->second));
        const LoadResult r = parse_obj(dataset::read_text(obj));
        const JsonValue* arr = j.ok() ? j.value.find("labels") : nullptr;
        if (!r.ok() || !arr || arr->array.size() != r.mesh.positions.size()) continue;
        if (!build_halfedge(r.mesh).ok()) continue;  // the measurement needs a clean manifold start
        ++scans;
        // 1. The scan's own holes.
        HoleFillParams prm;
        const auto loops = boundary_loops(r.mesh);
        std::size_t small = 0;
        for (const auto& l : loops) small += l.size() <= prm.max_loop;
        loops_per_scan.push_back(double(loops.size()));
        const auto t0 = Clock::now();
        const HoleFillResult own = fill_holes(r.mesh, prm);
        fill_ms.push_back(std::chrono::duration<double, std::milli>(Clock::now() - t0).count());
        small_loops += small, filled += own.filled, still_valid += build_halfedge(own.mesh).ok();
        // 2. Punch a disk out of a crown and fill it.
        std::vector<std::uint32_t> crown;
        for (std::uint32_t v = 0; v < r.mesh.positions.size(); ++v)
            if (arr->array[v].number != 0.0) crown.push_back(v);
        if (crown.empty()) continue;
        const Vec3 c = r.mesh.positions[crown[std::uniform_int_distribution<std::size_t>(0, crown.size() - 1)(rng)]];
        auto inside = [&](const Vec3& p) {
            const double dx = p.x - c.x, dy = p.y - c.y, dz = p.z - c.z;
            return dx * dx + dy * dy + dz * dz < radius * radius;
        };
        // Punch the disk out of the scan whose own holes are already filled, so the disk is the only loop the next fill
        // touches and every patch vertex measured belongs to it (measuring the whole result counted the patch over the
        // scan's open base too, centimetres away).
        const TriMesh& base = own.mesh;
        const std::size_t base_loops = boundary_loops(base).size();
        TriMesh holed;
        holed.positions = base.positions;
        std::vector<Vec3> removed;  // ground truth: vertices of the removed region
        for (const auto& t : base.triangles) {
            const bool cut = inside(base.positions[t[0]]) && inside(base.positions[t[1]]) && inside(base.positions[t[2]]);
            if (cut) {
                for (auto v : t) removed.push_back(base.positions[v]);
            } else {
                holed.triangles.push_back(t);
            }
        }
        if (removed.empty() || boundary_loops(holed).size() != base_loops + 1) continue;  // exactly one new hole
        ++punched;
        // Default loop limit (500): fills the punched disk (a few dozen edges) and the scan's own small holes, never the
        // long open cut (an earlier version raised the limit and spent its time triangulating the cut, O(n^3)).
        const HoleFillResult fix = fill_holes(holed, HoleFillParams{});
        punched_closed += boundary_loops(fix.mesh).size() == base_loops;
        std::vector<Vec3> patch;
        for (auto v : fix.patch_vertices) patch.push_back(fix.mesh.positions[v]);
        if (patch.empty()) continue;
        const auto d = nearest_point_distances(patch, removed, 0.25);
        double s = 0.0, mx = 0.0;
        for (double x : d) s += x, mx = std::max(mx, x);
        mean_err.push_back(s / double(d.size())), max_err.push_back(mx);
        std::fprintf(stderr, "\r%zu scans", scans);
    }
    std::fprintf(stderr, "\n");
    std::printf("scans: %zu manifold and labelled (of %zu sampled; the others fail the half-edge build)\n", scans, sampled);
    std::printf("own boundary loops (on these scans usually just the open base cut): median %.0f per scan; %zu of <= 500 edges, %zu filled; result a valid oriented"
                " manifold on %zu of %zu scans; fill time median %.0f ms\n", median(loops_per_scan), small_loops, filled, still_valid,
                scans, median(fill_ms));
    std::printf("punched %.1f mm disks on crowns: %zu, closed %zu; patch vertex distance to the removed surface: mean %.3f mm"
                " (median over scans), max %.3f mm (median), worst scan max %.3f mm\n", radius, punched, punched_closed,
                median(mean_err), median(max_err), max_err.empty() ? 0.0 : *std::max_element(max_err.begin(), max_err.end()));
    return 0;
}
