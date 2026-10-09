// Writes the synthetic molar (D96) as an OBJ plus ground-truth labels in the Teeth3DS JSON format (FDI 36 for tooth
// vertices, 0 for gingiva), for demos and the scripted README recording (web/scripts/record-demo.mjs).
//   export_demo <out.obj> <out.json>
#include <cmath>
#include <cstdio>
#include <fstream>

#include "core/generate.h"
#include "core/io.h"

using namespace dmw;

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: export_demo <out.obj> <out.json>\n");
        return 2;
    }
    const TriMesh m = make_synthetic_tooth();
    std::ofstream(argv[1]) << write_obj(m);
    std::ofstream json(argv[2]);
    json << "{\"jaw\":\"lower\",\"labels\":[";
    for (std::size_t v = 0; v < m.positions.size(); ++v)
        json << (v ? "," : "") << (std::hypot(m.positions[v].x, m.positions[v].y) < synthetic_tooth_radius ? 36 : 0);
    json << "]}\n";
    std::printf("wrote %s (%zu vertices, %zu faces) and %s\n", argv[1], m.positions.size(), m.triangles.size(), argv[2]);
    return 0;
}
