// Benchmarks for the README (D57). Times each pipeline stage on icospheres of growing size,
// median of repeated runs, and writes the meshes as binary STL so the WASM build can be timed
// on identical inputs (web/src/bench.ts). Usage: dmw_bench <output-dir>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include "core/curvature.h"
#include "core/generate.h"
#include "core/geodesic.h"
#include "core/halfedge.h"
#include "core/io.h"
#include "core/topology.h"
#include "core/weld.h"

using namespace dmw;
using Clock = std::chrono::steady_clock;

namespace {

// Median wall time in milliseconds over `runs` repetitions (median: robust to OS noise).
double median_ms(int runs, const std::function<void()>& work) {
    std::vector<double> t;
    for (int i = 0; i < runs; ++i) {
        const auto start = Clock::now();
        work();
        t.push_back(std::chrono::duration<double, std::milli>(Clock::now() - start).count());
    }
    std::sort(t.begin(), t.end());
    return t[t.size() / 2];
}

std::string fmt(double ms) {
    char buf[32];
    std::snprintf(buf, sizeof buf, ms < 10 ? "%.2f" : ms < 100 ? "%.1f" : "%.0f", ms);
    return buf;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: dmw_bench <output-dir>\n");
        return 2;
    }
    const std::filesystem::path out = argv[1];
    std::filesystem::create_directories(out);
#ifdef NDEBUG
    const char* build = "Release (NDEBUG)";
#else
    const char* build = "NOT Release: numbers are not meaningful";
#endif
#if defined(__clang__)
    const char* compiler = "clang " __clang_version__;
#elif defined(__GNUC__)
    const char* compiler = "gcc " __VERSION__;
#else
    const char* compiler = "unknown compiler";
#endif
    std::printf("build: %s, compiler: %s\n\n", build, compiler);
    std::printf("| mesh | V | F | STL MB | load STL (parse+weld) | topology | half-edge | curvature | "
                "load+analysis | curvature (iDT) | heat setup | heat query | Dijkstra |\n");
    std::printf("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n");

    for (std::uint32_t s = 3; s <= 7; ++s) {
        const TriMesh mesh = make_icosphere(s);
        const auto stl = write_stl_binary(mesh);
        const std::string name = "icosphere_s" + std::to_string(s) + ".stl";
        std::ofstream(out / name, std::ios::binary).write(reinterpret_cast<const char*>(stl.data()),
                                                          static_cast<std::streamsize>(stl.size()));
        const int runs = mesh.triangles.size() > 100000 ? 3 : 7;

        TriMesh loaded;
        const double t_load = median_ms(runs, [&] { loaded = weld_vertices(parse_stl(stl).mesh); });
        const double t_topo = median_ms(runs, [&] { (void)analyze_topology(loaded); });
        HalfEdgeMesh he;
        const double t_he = median_ms(runs, [&] { he = build_halfedge(loaded).mesh; });
        const double t_curv = median_ms(runs, [&] { (void)compute_curvature(he); });
        std::uint32_t flips = 0;
        const double t_idt = median_ms(runs, [&] { flips = compute_curvature(he, true).intrinsic_flips; });
        const double total = t_load + t_topo + t_he + t_curv;  // what the viewer does on load

        // The envelope factorization grows ~n^1.5 in memory (D49): stop at s = 6 (41k vertices).
        std::string setup = "skipped (D49)", query = "-", dij = "-";
        if (s <= 6) {
            const std::uint32_t src[1] = {0};
            setup = fmt(median_ms(s == 6 ? 1 : runs, [&] { HeatGeodesics g(he); }));
            const HeatGeodesics g(he);
            query = fmt(median_ms(runs, [&] { (void)g.distance(src); }));
            dij = fmt(median_ms(runs, [&] { (void)dijkstra_distance(he, src); }));
        }
        std::printf("| icosphere s=%u | %zu | %zu | %.2f | %s | %s | %s | %s | **%s** | %s (%u flips) | %s | %s | %s |\n", s,
                    loaded.positions.size(), loaded.triangles.size(), static_cast<double>(stl.size()) / 1e6,
                    fmt(t_load).c_str(), fmt(t_topo).c_str(), fmt(t_he).c_str(), fmt(t_curv).c_str(),
                    fmt(total).c_str(), fmt(t_idt).c_str(), flips, setup.c_str(), query.c_str(), dij.c_str());
        std::fflush(stdout);
    }
    std::printf("\nAll times in milliseconds, median of 7 runs (3 for F > 100k; heat setup at s=6: 1 run).\n");
    return 0;
}
