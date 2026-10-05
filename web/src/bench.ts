// WASM benchmark on the same STL files dmw_bench writes (run it first). Run: npm run bench
// Times the full viewer path through the module: copy bytes into WASM memory, parse, weld,
// topology, half-edge, curvature, copy results out (Dmw.loadFile), then geodesic queries.
import { readFileSync } from "node:fs";
import { Dmw } from "./dmw.ts";

const dir = process.argv[2] ?? "../build-release/bench";
const median = (xs: number[]) => [...xs].sort((a, b) => a - b)[Math.floor(xs.length / 2)];
const fmt = (ms: number) => (ms < 10 ? ms.toFixed(2) : ms < 100 ? ms.toFixed(1) : ms.toFixed(0));

const dmw = await Dmw.create();
console.log(`node ${process.version}, WASM (emsdk 6.0.11, CMake Release: -O3 compile and link)\n`);
console.log("| mesh | V | load+analysis (WASM) | heat setup+first query | heat query |");
console.log("|---|---:|---:|---:|---:|");
for (const s of [3, 4, 5, 6, 7]) {
  const bytes = new Uint8Array(readFileSync(`${dir}/icosphere_s${s}.stl`));
  const runs = s >= 7 ? 3 : 7;
  const loads: number[] = [];
  let vertices = 0;
  for (let i = 0; i < runs; ++i) {
    const d = dmw.loadFile(bytes, "stl");
    loads.push(d.millis);
    vertices = d.stats.vertices;
  }
  let first = "skipped (D49)", query = "-";
  if (s <= 6) {
    first = fmt(dmw.geodesic(0).millis); // factorization happens here, once per mesh
    query = fmt(median(Array.from({ length: runs }, () => dmw.geodesic(0).millis)));
  }
  console.log(`| icosphere s=${s} | ${vertices} | **${fmt(median(loads))}** | ${first} | ${query} |`);
}
dmw.dispose();
