// Cross-target check: the WASM build must give the same answers as the native tests.
// Run: npm run smoke
import assert from "node:assert/strict";
import { Dmw } from "./dmw.ts";

const dmw = await Dmw.create();
const close = (a: number, b: number, tol = 1e-9) => assert.ok(Math.abs(a - b) < tol, `${a} != ${b}`);

// Topology and Gauss-Bonnet through the generators.
const torus = dmw.generate("torus").stats;
assert.equal(torus.components[0].genus, 1);
close(torus.curvature!.totalAngleDefect, 0);

const sphere = dmw.generate("icosphere").stats;
assert.equal(sphere.components[0].genus, 0);
close(sphere.curvature!.totalAngleDefect, 4 * Math.PI);

const handleData = dmw.generate("plate_handle");
assert.equal(handleData.stats.handleLoops.count, 2); // M8b: 2g loops, the boundary is not one
assert.equal(handleData.handleEdges.length % 2, 0);
assert.equal(dmw.generate("torus").stats.handleLoops.count, 2);
const handle = handleData.stats;
assert.deepEqual([handle.components[0].chi, handle.components[0].b, handle.components[0].genus], [-1, 1, 1]);

const defects = dmw.generate("defects");
assert.equal(defects.stats.edgeKinds.nonmanifold, 1);
assert.deepEqual(Array.from(defects.nonmanifoldVertices), [24]);
// D69: faces at the defects are excluded from the half-edge analysis view, so curvature still runs;
// the topology report above still describes the full input.
assert.equal(defects.stats.halfedge.ok, true);
assert.ok(defects.stats.excludedFaces > 0 && defects.stats.excludedFaces < 20);
assert.ok(defects.stats.curvature !== null);

// Heat-method geodesics through WASM: icosahedron vertices 0 and 3 are antipodal on the
// unit sphere, so their distance is pi (measured heat error at this resolution: < 3e-2).
dmw.generate("icosphere");
const geo = dmw.geodesic(0);
close(geo.distance[0], 0, 1e-12);
close(geo.distance[3], Math.PI, 5e-2);
// Cusp detection through WASM runs and returns vertex ids within range.
const sph = dmw.generate("icosphere");
const tips = dmw.detectCusps().vertices;
assert.ok(tips.every((v) => v < sph.stats.vertices));

// Intrinsic Delaunay (M6b) through WASM: the brick grid has obtuse triangles, needs flips, and
// geodesics stay accurate; Gauss-Bonnet is untouched (angle defects are intrinsic).
const brickCot = dmw.generate("brick");
const brickIdt = dmw.setIntrinsicDelaunay(true);
assert.ok(brickIdt.stats.flips > 0);
close(brickIdt.stats.curvature!.totalAngleDefect, brickCot.stats.curvature!.totalAngleDefect, 1e-9);
const center = 15 * 31 + 15, corner = 0; // 30x30 grid: center vertex and corner (0, 0)
const dBrick = dmw.geodesic(center).distance;
const p = brickIdt.positions;
close(dBrick[corner], Math.hypot(p[3 * corner] - p[3 * center], p[3 * corner + 1] - p[3 * center + 1]), 3e-2);
dmw.setIntrinsicDelaunay(false);

// The shared-memory input path: a binary STL tetrahedron built byte by byte.
const tris = [
  [0, 0, 0, 0, 1, 0, 1, 0, 0],
  [0, 0, 0, 1, 0, 0, 0, 0, 1],
  [0, 0, 0, 0, 0, 1, 0, 1, 0],
  [1, 0, 0, 0, 1, 0, 0, 0, 1],
];
const stl = new DataView(new ArrayBuffer(84 + 50 * tris.length));
stl.setUint32(80, tris.length, true); // little-endian
tris.forEach((t, i) => t.forEach((c, k) => stl.setFloat32(84 + 50 * i + 12 + 4 * k, c, true)));
const tet = dmw.loadFile(new Uint8Array(stl.buffer), "stl");
assert.equal(tet.stats.vertices, 4); // 12 soup corners welded to 4
assert.equal(tet.stats.components[0].genus, 0);
close(tet.stats.curvature!.totalAngleDefect, 4 * Math.PI);

// OBJ text path, and a loader error surfacing as an exception.
const obj = dmw.loadFile(new TextEncoder().encode("v 0 0 0\nv 1 0 0\nv 0 1 0\nf 1 2 3\n"), "obj");
assert.equal(obj.stats.components[0].b, 1);
assert.throws(() => dmw.loadFile(new TextEncoder().encode("f 1 2 3\n"), "obj"), /line 1/);

dmw.dispose();
console.log("smoke ok: WASM topology, curvature and STL/OBJ input match the native results; geodesics, intrinsic Delaunay, cusp detection OK");
