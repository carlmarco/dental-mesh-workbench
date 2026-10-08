// Undercut analysis on real scans (M10b, D90): per tooth crown (faces whose three vertices carry the same FDI
// label in the Teeth3DS ground truth; labels only define the regions), the undercut fraction along the occlusal
// axis and along the best path of insertion within 25 degrees; plus one common path for all teeth of the arch
// (the removable-appliance / aligner case). Reports timing and ray throughput.
//   undercut_eval <scan-dir> [stride]
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "core/bvh.h"
#include "core/cusps.h"
#include "core/io.h"
#include "core/topology.h"
#include "core/undercut.h"
#include "dataset.h"

using namespace dmw;
using Clock = std::chrono::steady_clock;

namespace {
double ms_since(Clock::time_point t) { return std::chrono::duration<double, std::milli>(Clock::now() - t).count(); }
double median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}
const char* type_name(int unit) {
    switch (unit) {
        case 1: case 2: return "incisors";
        case 3: return "canines";
        case 4: case 5: return "premolars";
        default: return "molars";
    }
}
}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: undercut_eval <scan-dir> [stride]\n");
        return 2;
    }
    const std::size_t stride = argc > 2 ? std::stoul(argv[2]) : 15;
    const auto objs = dataset::index_files({argv[1]}, ".obj");
    const auto labels = dataset::index_files({argv[1]}, ".json");
    struct PerType {
        std::vector<double> occlusal, best, tilt, ms;
    };
    std::map<std::string, PerType> types;
    std::vector<double> arch_occlusal, arch_best, arch_tilt, build_ms, rays_per_s;
    std::size_t index = 0, scans = 0, teeth = 0;
    for (const auto& [stem, obj] : objs) {
        if (index++ % stride != 0) continue;
        const auto lab = labels.find(stem);
        if (lab == labels.end()) continue;
        const JsonResult j = parse_json(dataset::read_text(lab->second));
        const LoadResult r = parse_obj(dataset::read_text(obj));
        const JsonValue* arr = j.ok() ? j.value.find("labels") : nullptr;
        if (!r.ok() || !arr || arr->array.size() != r.mesh.positions.size()) continue;
        const AnalysisMesh a = manifold_analysis_mesh(r.mesh);
        if (!a.manifold) continue;
        const TriMesh& m = r.mesh;
        const Vec3 axis = occlusal_axis(a.halfedge, largest_component_mask(a.halfedge));
        auto t0 = Clock::now();
        const Bvh bvh(m.positions, m.triangles);
        build_ms.push_back(ms_since(t0));
        // Face regions per FDI tooth.
        std::map<int, std::vector<std::uint8_t>> region;
        std::vector<std::uint8_t> all_teeth(m.triangles.size(), 0);
        for (std::size_t f = 0; f < m.triangles.size(); ++f) {
            const int l0 = int(arr->array[m.triangles[f][0]].number), l1 = int(arr->array[m.triangles[f][1]].number),
                      l2 = int(arr->array[m.triangles[f][2]].number);
            if (l0 == 0 || l0 != l1 || l0 != l2) continue;
            auto& reg = region[l0];
            if (reg.empty()) reg.assign(m.triangles.size(), 0);
            reg[f] = 1, all_teeth[f] = 1;
        }
        for (const auto& [fdi, reg] : region) {
            t0 = Clock::now();
            const double occ = undercut_map(m, bvh, axis, reg).fraction();
            const InsertionAxis best = best_insertion_axis(m, bvh, reg, axis, 25.0);
            const double ms = ms_since(t0);
            std::size_t faces = 0;
            for (auto x : reg) faces += x;
            rays_per_s.push_back(double(faces) * double(best.evaluations + 1) / (ms / 1000.0));
            PerType& pt = types[type_name(fdi % 10)];
            pt.occlusal.push_back(100.0 * occ), pt.best.push_back(100.0 * best.result.fraction());
            pt.tilt.push_back(std::acos(std::clamp(best.axis.x * axis.x + best.axis.y * axis.y + best.axis.z * axis.z, -1.0, 1.0)) * 180.0 /
                              std::acos(-1.0));
            pt.ms.push_back(ms);
            ++teeth;
        }
        const InsertionAxis arch = best_insertion_axis(m, bvh, all_teeth, axis, 25.0);
        arch_occlusal.push_back(100.0 * undercut_map(m, bvh, axis, all_teeth).fraction());
        arch_best.push_back(100.0 * arch.result.fraction());
        arch_tilt.push_back(std::acos(std::clamp(arch.axis.x * axis.x + arch.axis.y * axis.y + arch.axis.z * axis.z, -1.0, 1.0)) * 180.0 /
                            std::acos(-1.0));
        ++scans;
        std::fprintf(stderr, "\r%zu scans, %zu teeth", scans, teeth);
    }
    std::fprintf(stderr, "\n");
    std::printf("scans: %zu, teeth: %zu; BVH build median %.0f ms; ray throughput median %.2f M rays/s (one thread)\n\n", scans, teeth,
                median(build_ms), median(rays_per_s) / 1e6);
    std::printf("| tooth type | teeth | undercut along occlusal axis (median %% of crown) | best axis within 25 deg (median %%) | median tilt (deg) | median ms per tooth |\n");
    std::printf("|---|---:|---:|---:|---:|---:|\n");
    for (const auto& [name, pt] : types)
        std::printf("| %s | %zu | %.1f | %.1f | %.1f | %.0f |\n", name.c_str(), pt.occlusal.size(), median(pt.occlusal), median(pt.best),
                    median(pt.tilt), median(pt.ms));
    std::printf("\nwhole arch, one common path for all crowns: undercut %.1f%% along the occlusal axis -> %.1f%% at the best axis"
                " (median tilt %.1f deg)\n", median(arch_occlusal), median(arch_best), median(arch_tilt));
    return 0;
}
