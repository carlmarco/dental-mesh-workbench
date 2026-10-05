// Typed wrapper over the embind module: the only file that touches WASM memory directly.
import createDmw from "./wasm/dmw.mjs";
import type { MainModule, Session } from "./wasm/dmw.mjs";

export interface ComponentStats {
  V: number;
  E: number;
  F: number;
  chi: number;
  b: number;
  manifold: boolean;
  orientable: boolean;
  consistent: boolean;
  genus: number | null;
  betti: [number, number, number] | null;
}

export interface Stats {
  vertices: number;
  faces: number;
  edges: number;
  edgeKinds: { boundary: number; manifold: number; misoriented: number; nonmanifold: number };
  invalidFaces: number;
  duplicateFaces: number;
  isolatedVertices: number;
  nonmanifoldVertices: number;
  components: ComponentStats[];
  intrinsicDelaunay: boolean;
  flips: number; // intrinsic edge flips (when intrinsicDelaunay)
  halfedge: { ok: boolean; error: string };
  curvature: null | {
    totalAngleDefect: number;
    meanRange: [number, number] | null;
    gaussianRange: [number, number] | null;
    degenerateFaces: number;
  };
}

// JS-owned copies of everything the viewer needs (safe across later WASM calls).
export interface MeshData {
  positions: Float32Array; // xyz per vertex
  indices: Uint32Array; // 3 per face
  mean: Float32Array; // per vertex, NaN = no data
  gaussian: Float32Array; // per vertex, NaN = no data
  faceComponent: Uint32Array; // per face, 0xFFFFFFFF = excluded
  boundaryEdges: Uint32Array; // vertex pairs
  nonmanifoldEdges: Uint32Array;
  misorientedEdges: Uint32Array;
  nonmanifoldVertices: Uint32Array;
  stats: Stats;
  millis: number; // wall time of load/generate + analysis + copies, measured here
}

export const PRESETS = [
  "plate_handle",
  "grid_holes",
  "defects",
  "brick",
  "icosphere",
  "torus",
  "cylinder",
  "grid",
  "mobius",
] as const;
export type Preset = (typeof PRESETS)[number];

export class Dmw {
  private readonly module: MainModule;
  private readonly session: Session;

  private constructor(module: MainModule) {
    this.module = module;
    this.session = new module.Session();
  }

  static async create(): Promise<Dmw> {
    return new Dmw(await createDmw());
  }

  generate(name: Preset): MeshData {
    return this.run(() => this.session.generate(name));
  }

  // Shared-memory input (D41): allocate inside WASM, copy the bytes in, hand C++ the offset.
  loadFile(bytes: Uint8Array, format: "stl" | "obj", weldEpsilon = 0): MeshData {
    return this.run(() => {
      const ptr = this.module._malloc(bytes.length);
      try {
        // Read HEAPU8 *after* malloc: if malloc grew memory, the old view is detached.
        this.module.HEAPU8.set(bytes, ptr);
        return this.session.load(ptr, bytes.length, format, weldEpsilon);
      } finally {
        this.module._free(ptr);
      }
    });
  }

  // Switch curvature and geodesics to the intrinsic Delaunay Laplacian (M6b) or back.
  setIntrinsicDelaunay(on: boolean): MeshData {
    return this.run(() => (this.session.setIntrinsicDelaunay(on), ""));
  }

  // Heat-method distance from one vertex. The first call on a mesh also builds and factors
  // the solver, so it is slower; later calls reuse the factorization.
  geodesic(source: number): { distance: Float32Array; millis: number } {
    const t0 = performance.now();
    const error = this.session.geodesic(source);
    if (error) throw new Error(error);
    const distance = (this.session.distance() as Float32Array).slice();
    return { distance, millis: performance.now() - t0 };
  }

  // embind objects live in C++ memory and are not garbage-collected.
  dispose(): void {
    this.session.delete();
  }

  private run(call: () => string): MeshData {
    const t0 = performance.now();
    const error = call();
    if (error) throw new Error(error);
    const s = this.session;
    // Each accessor returns a view aliasing WASM memory; slice() copies it out at once.
    const data: MeshData = {
      positions: (s.positions() as Float32Array).slice(),
      indices: (s.indices() as Uint32Array).slice(),
      mean: (s.mean() as Float32Array).slice(),
      gaussian: (s.gaussian() as Float32Array).slice(),
      faceComponent: (s.faceComponent() as Uint32Array).slice(),
      boundaryEdges: (s.boundaryEdges() as Uint32Array).slice(),
      nonmanifoldEdges: (s.nonmanifoldEdges() as Uint32Array).slice(),
      misorientedEdges: (s.misorientedEdges() as Uint32Array).slice(),
      nonmanifoldVertices: (s.nonmanifoldVertices() as Uint32Array).slice(),
      stats: JSON.parse(s.stats()) as Stats,
      millis: 0,
    };
    data.millis = performance.now() - t0;
    return data;
  }
}
