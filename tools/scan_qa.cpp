// Topology and timing QA over real intraoral scans (D62). Reads every .obj under the given
// directories, runs the viewer's analysis pipeline, and prints per-scan rows to a CSV plus an
// aggregate summary. Scans are licensed for local use only (D60): the CSV stays local.
// Usage: scan_qa <out.csv> <dir> [<dir> ...]
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <string>
#include <vector>

#include "core/curvature.h"
#include "core/halfedge.h"
#include "core/io.h"
#include "core/topology.h"

using namespace dmw;
namespace fs = std::filesystem;
using Clock = std::chrono::steady_clock;

namespace {

double ms_since(Clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

double percentile(std::vector<double> v, double q) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[static_cast<std::size_t>(q * static_cast<double>(v.size() - 1))];
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: scan_qa <out.csv> <dir> [<dir> ...]\n");
        return 2;
    }
    std::vector<fs::path> files;
    for (int i = 2; i < argc; ++i) {
        for (const auto& e : fs::recursive_directory_iterator(argv[i])) {
            if (e.is_regular_file() && e.path().extension() == ".obj") files.push_back(e.path());
        }
    }
    std::sort(files.begin(), files.end());
    std::ofstream csv(argv[1]);
    csv << "scan,vertices,faces,read_ms,parse_ms,topology_ms,halfedge_ms,curvature_ms,components,"
           "largest_V,largest_chi,largest_b,largest_genus,largest_manifold,nonmanifold_edges,"
           "nonmanifold_vertices,misoriented_edges,invalid_faces,duplicate_faces,isolated_vertices,"
           "halfedge_ok,zero_area_faces\n";

    std::vector<double> parse, topo, he, curv, verts;
    std::map<std::string, int> genus_hist, boundary_hist, comp_hist;
    int halfedge_fail = 0, any_nm = 0, any_dup = 0, any_invalid = 0, multi_comp = 0, n = 0;
    for (const auto& path : files) {
        auto t0 = Clock::now();
        std::ifstream in(path, std::ios::binary);
        const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        const double t_read = ms_since(t0);
        t0 = Clock::now();
        const LoadResult r = parse_obj(text);
        const double t_parse = ms_since(t0);
        if (!r.ok()) {
            std::fprintf(stderr, "%s: %s\n", path.filename().c_str(), r.error.c_str());
            continue;
        }
        t0 = Clock::now();
        const TopologyReport topo_r = analyze_topology(r.mesh);
        const double t_topo = ms_since(t0);
        t0 = Clock::now();
        const BuildResult hb = build_halfedge(without_excluded_faces(r.mesh, topo_r));  // D63
        const double t_he = ms_since(t0);
        double t_curv = 0.0;
        std::size_t zero_area = 0;
        if (hb.ok()) {
            t0 = Clock::now();
            zero_area = compute_curvature(hb.mesh).degenerate_faces.size();
            t_curv = ms_since(t0);
        }
        // "Largest" component by vertex count: the arch; small ones are scan debris.
        const ComponentTopology* big = nullptr;
        for (const auto& c : topo_r.components) {
            if (!big || c.num_vertices > big->num_vertices) big = &c;
        }
        const std::size_t nm_e = topo_r.count(EdgeKind::NonManifold), mis = topo_r.count(EdgeKind::Misoriented);
        csv << path.stem().string() << ',' << r.mesh.positions.size() << ',' << r.mesh.triangles.size() << ','
            << t_read << ',' << t_parse << ',' << t_topo << ',' << t_he << ',' << t_curv << ','
            << topo_r.components.size() << ',' << (big ? big->num_vertices : 0) << ','
            << (big ? big->euler_characteristic : 0) << ',' << (big ? big->boundary_loops : 0) << ','
            << (big && big->genus ? std::to_string(*big->genus) : "") << ',' << (big && big->manifold) << ','
            << nm_e << ',' << topo_r.nonmanifold_vertices.size() << ',' << mis << ','
            << topo_r.invalid_faces.size() << ',' << topo_r.duplicate_faces.size() << ','
            << topo_r.isolated_vertices.size() << ',' << hb.ok() << ',' << zero_area << '\n';

        ++n;
        parse.push_back(t_parse), topo.push_back(t_topo), he.push_back(t_he), verts.push_back(double(r.mesh.positions.size()));
        if (hb.ok()) curv.push_back(t_curv); else ++halfedge_fail;
        genus_hist[big && big->genus ? std::to_string(*big->genus) : "n/a"]++;
        boundary_hist[big ? std::to_string(std::min<std::uint32_t>(big->boundary_loops, 10)) : "n/a"]++;
        comp_hist[std::to_string(std::min<std::size_t>(topo_r.components.size(), 5))]++;
        multi_comp += topo_r.components.size() > 1;
        any_nm += (nm_e + topo_r.nonmanifold_vertices.size()) > 0;
        any_dup += !topo_r.duplicate_faces.empty();
        any_invalid += !topo_r.invalid_faces.empty();
    }
    auto print_hist = [](const char* name, const std::map<std::string, int>& h) {
        std::printf("%s:", name);
        for (const auto& [k, v] : h) std::printf("  %s->%d", k.c_str(), v);
        std::printf("\n");
    };
    std::printf("scans analysed: %d\n", n);
    std::printf("vertices: median %.0f, min %.0f, max %.0f\n", percentile(verts, 0.5), percentile(verts, 0.0), percentile(verts, 1.0));
    std::printf("timing ms (median / p95): parse %.0f / %.0f, topology %.0f / %.0f, half-edge %.0f / %.0f, curvature %.0f / %.0f\n",
                percentile(parse, 0.5), percentile(parse, 0.95), percentile(topo, 0.5), percentile(topo, 0.95),
                percentile(he, 0.5), percentile(he, 0.95), percentile(curv, 0.5), percentile(curv, 0.95));
    std::printf("half-edge build failed (non-manifold / misoriented input): %d\n", halfedge_fail);
    std::printf("scans with >1 component: %d, with non-manifold elements: %d, duplicate faces: %d, invalid faces: %d\n",
                multi_comp, any_nm, any_dup, any_invalid);
    print_hist("components per scan (5 = 5+)", comp_hist);
    print_hist("largest component boundary loops (10 = 10+)", boundary_hist);
    print_hist("largest component genus", genus_hist);
    return 0;
}
