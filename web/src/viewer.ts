import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { LineMaterial } from "three/addons/lines/LineMaterial.js";
import { LineSegments2 } from "three/addons/lines/LineSegments2.js";
import { LineSegmentsGeometry } from "three/addons/lines/LineSegmentsGeometry.js";

import { categorical, diverging, NO_DATA, robustRange, type RGB, sequential } from "./colormap.ts";
import type { MeshData } from "./dmw.ts";

export type Overlay = "shaded" | "mean" | "gaussian" | "components" | "geodesic";

export interface Layers {
  boundary: boolean;
  nonmanifold: boolean;
  misoriented: boolean;
  vertices: boolean;
  handles: boolean;
  cusps: boolean;
  margin: boolean;
  wireframe: boolean;
}

export const LAYER_COLORS = {
  boundary: 0x2f80ed,
  nonmanifold: 0xe0342b,
  misoriented: 0xf2a900,
  vertices: 0xc026d3,
  handles: 0x16a34a,
  cusps: 0xf97316,
  truth: 0x22d3ee,
  margin: 0xfacc15,
} as const;

const EXCLUDED = 0xffffffff;
const SHADED: RGB = [0.72, 0.75, 0.8];
export const ISOLINES = 20; // isolines drawn at multiples of (max distance / ISOLINES)

// 1D texture: sequential colors with a dark band at each isoline. Used with a texture
// coordinate u = distance / max: u interpolates linearly across each triangle and the texture
// is sampled per pixel, so isolines stay crisp even on coarse meshes (per-vertex colors would
// smear them across whole triangles).
function isolineTexture(): THREE.DataTexture {
  const width = 1024;
  const data = new Uint8Array(width * 4);
  for (let i = 0; i < width; ++i) {
    const u = i / (width - 1);
    const line = Math.abs(u * ISOLINES - Math.round(u * ISOLINES)) < 0.06;
    const [r, g, b] = sequential(u).map((c) => (line ? c * 0.35 : c));
    data.set([r * 255, g * 255, b * 255, 255], 4 * i);
  }
  const tex = new THREE.DataTexture(data, width, 1);
  tex.colorSpace = THREE.SRGBColorSpace;
  tex.magFilter = THREE.LinearFilter;
  tex.minFilter = THREE.LinearFilter;
  tex.needsUpdate = true;
  return tex;
}

export class Viewer {
  private readonly renderer: THREE.WebGLRenderer;
  private readonly scene = new THREE.Scene();
  private readonly camera = new THREE.PerspectiveCamera(40, 1, 0.001, 1000);
  private readonly controls: OrbitControls;
  private readonly content = new THREE.Group();
  private readonly container: HTMLElement;
  private lineMaterials: LineMaterial[] = [];
  private data: MeshData | null = null;
  private distance: Float32Array | null = null;
  private cuspPoints: Float32Array | null = null;  // detected tips, xyz
  private truthPoints: Float32Array | null = null; // loaded landmarks, xyz
  private marginEdges: Uint32Array | null = null;
  private truthMarginEdges: Uint32Array | null = null;
  private source: number | null = null;
  private surface: THREE.Mesh | null = null;
  private readonly isolines = isolineTexture();
  onPick: ((vertex: number) => void) | null = null;

  constructor(container: HTMLElement) {
    this.container = container;
    this.renderer = new THREE.WebGLRenderer({ antialias: true });
    this.renderer.setPixelRatio(window.devicePixelRatio);
    container.appendChild(this.renderer.domElement);

    this.camera.up.set(0, 0, 1); // our generators are z-up
    this.scene.add(this.camera, this.content);
    this.scene.add(new THREE.HemisphereLight(0xffffff, 0x445566, 1.6));
    const headlight = new THREE.DirectionalLight(0xffffff, 1.6); // follows the camera
    headlight.position.set(1, 1, 2);
    this.camera.add(headlight);

    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = true;

    // Click (press and release without dragging) picks the nearest vertex of the hit face.
    let down: { x: number; y: number } | null = null;
    const canvas = this.renderer.domElement;
    canvas.addEventListener("pointerdown", (e) => (down = { x: e.clientX, y: e.clientY }));
    canvas.addEventListener("pointerup", (e) => {
      if (down && Math.hypot(e.clientX - down.x, e.clientY - down.y) < 5) this.pick(e);
      down = null;
    });

    new ResizeObserver(() => this.resize()).observe(container);
    this.resize();
    this.renderer.setAnimationLoop(() => {
      this.controls.update();
      this.renderer.render(this.scene, this.camera);
    });
  }

  setBackground(color: string): void {
    this.scene.background = new THREE.Color(color);
  }

  // New mesh: rebuild and re-frame the camera. Returns the color range used (curvature modes).
  show(data: MeshData, overlay: Overlay, layers: Layers): number | null {
    this.data = data;
    const range = this.rebuild(overlay, layers);
    this.frame();
    return range;
  }

  setCusps(detected: Float32Array | null, truth: Float32Array | null): void {
    this.cuspPoints = detected;
    this.truthPoints = truth;
  }

  setMargin(predicted: Uint32Array | null, truth: Uint32Array | null): void {
    this.marginEdges = predicted;
    this.truthMarginEdges = truth;
  }

  // Distance field for the "geodesic" overlay (null clears it).
  setGeodesic(distance: Float32Array | null, source: number | null): void {
    this.distance = distance;
    this.source = source;
  }

  // Same mesh, different overlay/layers: keep the camera where the user put it.
  restyle(overlay: Overlay, layers: Layers): number | null {
    return this.data ? this.rebuild(overlay, layers) : null;
  }

  private rebuild(overlay: Overlay, layers: Layers): number | null {
    this.clear();
    const d = this.data!;
    let range: number | null = null;

    // Mesh. Per-face coloring (components) needs non-indexed geometry, since an indexed
    // vertex is shared by faces of different colors; per-vertex fields stay indexed.
    const geometry = new THREE.BufferGeometry();
    let map: THREE.Texture | null = null;
    if (overlay === "geodesic" && this.distance) {
      let max = 0;
      for (const x of this.distance) if (Number.isFinite(x)) max = Math.max(max, x);
      range = max;
      const uv = new Float32Array((d.positions.length / 3) * 2);
      this.distance.forEach((x, v) => uv.set([Number.isFinite(x) && max > 0 ? Math.max(0, x) / max : 0, 0.5], 2 * v));
      geometry.setAttribute("position", new THREE.BufferAttribute(d.positions, 3));
      geometry.setAttribute("uv", new THREE.BufferAttribute(uv, 2));
      geometry.setIndex(new THREE.BufferAttribute(d.indices, 1));
      map = this.isolines;
    } else if (overlay === "components") {
      const nf = d.indices.length / 3;
      const pos = new Float32Array(nf * 9);
      const col = new Float32Array(nf * 9);
      for (let f = 0; f < nf; ++f) {
        const c = d.faceComponent[f] === EXCLUDED ? NO_DATA : categorical(d.faceComponent[f]);
        for (let k = 0; k < 3; ++k) {
          const v = d.indices[3 * f + k];
          pos.set(d.positions.subarray(3 * v, 3 * v + 3), 9 * f + 3 * k);
          col.set(c, 9 * f + 3 * k);
        }
      }
      geometry.setAttribute("position", new THREE.BufferAttribute(pos, 3));
      geometry.setAttribute("color", new THREE.BufferAttribute(col, 3));
    } else {
      const nv = d.positions.length / 3;
      const col = new Float32Array(nv * 3);
      const field = overlay === "mean" ? d.mean : overlay === "gaussian" ? d.gaussian : null;
      if (field) range = robustRange(field);
      for (let v = 0; v < nv; ++v) {
        const x = field ? field[v] : 0;
        col.set(!field ? SHADED : Number.isFinite(x) ? diverging(x / range!) : NO_DATA, 3 * v);
      }
      geometry.setAttribute("position", new THREE.BufferAttribute(d.positions, 3));
      geometry.setAttribute("color", new THREE.BufferAttribute(col, 3));
      geometry.setIndex(new THREE.BufferAttribute(d.indices, 1));
    }
    geometry.computeVertexNormals();
    const material = new THREE.MeshStandardMaterial({
      vertexColors: map === null,
      map,
      side: THREE.DoubleSide, // open meshes show their back faces
      roughness: 0.65,
      metalness: 0.0,
      polygonOffset: true, // push faces back so edge overlays don't z-fight
      polygonOffsetFactor: 1,
      polygonOffsetUnits: 1,
    });
    this.surface = new THREE.Mesh(geometry, material);
    this.content.add(this.surface);

    if (layers.wireframe) {
      const wire = new THREE.LineSegments(
        new THREE.WireframeGeometry(geometry),
        new THREE.LineBasicMaterial({ color: 0x000000, transparent: true, opacity: 0.18 }),
      );
      this.content.add(wire);
    }
    // Defect overlays: fat lines (WebGL ignores lineWidth > 1, hence LineSegments2).
    // Defects that must never be hidden are drawn on top of everything (no depth test).
    if (layers.boundary) this.addEdges(d.boundaryEdges, LAYER_COLORS.boundary, 2.5, false);
    if (layers.misoriented) this.addEdges(d.misorientedEdges, LAYER_COLORS.misoriented, 4, true);
    if (layers.nonmanifold) this.addEdges(d.nonmanifoldEdges, LAYER_COLORS.nonmanifold, 5, true);
    if (layers.handles) this.addEdges(d.handleEdges, LAYER_COLORS.handles, 4, true);
    if (layers.margin && this.truthMarginEdges?.length) this.addEdges(this.truthMarginEdges, LAYER_COLORS.truth, 2.5, true);
    if (layers.margin && this.marginEdges?.length) this.addEdges(this.marginEdges, LAYER_COLORS.margin, 3.5, true);
    if (layers.cusps && this.cuspPoints?.length) this.addXyzPoints(this.cuspPoints, LAYER_COLORS.cusps, 11);
    if (layers.cusps && this.truthPoints?.length) this.addXyzPoints(this.truthPoints, LAYER_COLORS.truth, 7);
    if (overlay === "geodesic" && this.source !== null) {
      this.addPoints(Uint32Array.of(this.source), 0xffffff, 14);
    }
    if (layers.vertices && d.nonmanifoldVertices.length > 0) {
      const pts = new Float32Array(d.nonmanifoldVertices.length * 3);
      d.nonmanifoldVertices.forEach((v, i) => pts.set(d.positions.subarray(3 * v, 3 * v + 3), 3 * i));
      const g = new THREE.BufferGeometry();
      g.setAttribute("position", new THREE.BufferAttribute(pts, 3));
      const m = new THREE.PointsMaterial({ color: LAYER_COLORS.vertices, size: 12, sizeAttenuation: false });
      m.depthTest = false; // always visible: a single vertex is easy to hide behind faces
      this.content.add(new THREE.Points(g, m));
    }
    return range;
  }

  private addPoints(ids: Uint32Array, color: number, size: number): void {
    const d = this.data!;
    const pts = new Float32Array(ids.length * 3);
    ids.forEach((v, i) => pts.set(d.positions.subarray(3 * v, 3 * v + 3), 3 * i));
    const g = new THREE.BufferGeometry();
    g.setAttribute("position", new THREE.BufferAttribute(pts, 3));
    const m = new THREE.PointsMaterial({ color, size, sizeAttenuation: false });
    m.depthTest = false;
    this.content.add(new THREE.Points(g, m));
  }

  private addXyzPoints(xyz: Float32Array, color: number, size: number): void {
    const g = new THREE.BufferGeometry();
    g.setAttribute("position", new THREE.BufferAttribute(xyz, 3));
    const m = new THREE.PointsMaterial({ color, size, sizeAttenuation: false });
    m.depthTest = false;
    this.content.add(new THREE.Points(g, m));
  }

  private pick(e: PointerEvent): void {
    if (!this.surface || !this.data || !this.onPick) return;
    const rect = this.renderer.domElement.getBoundingClientRect();
    const ndc = new THREE.Vector2(((e.clientX - rect.left) / rect.width) * 2 - 1, -((e.clientY - rect.top) / rect.height) * 2 + 1);
    const ray = new THREE.Raycaster();
    ray.setFromCamera(ndc, this.camera);
    const hit = ray.intersectObject(this.surface)[0];
    if (!hit?.face) return;
    // Non-indexed geometry (components view) numbers corners, not vertices: map back.
    const corner = (c: number) => (this.surface!.geometry.index ? c : this.data!.indices[c]);
    let best = -1, bestDist = Infinity;
    for (const c of [hit.face.a, hit.face.b, hit.face.c]) {
      const v = corner(c);
      const p = this.data.positions;
      const dist = hit.point.distanceToSquared(new THREE.Vector3(p[3 * v], p[3 * v + 1], p[3 * v + 2]));
      if (dist < bestDist) (best = v), (bestDist = dist);
    }
    if (best >= 0) this.onPick(best);
  }

  private addEdges(pairs: Uint32Array, color: number, width: number, onTop: boolean): void {
    if (pairs.length === 0) return;
    const d = this.data!;
    const xyz = new Float32Array(pairs.length * 3);
    pairs.forEach((v, i) => xyz.set(d.positions.subarray(3 * v, 3 * v + 3), 3 * i));
    const geometry = new LineSegmentsGeometry();
    geometry.setPositions(xyz);
    const material = new LineMaterial({ color, linewidth: width }); // width in pixels
    material.resolution.set(this.container.clientWidth, this.container.clientHeight);
    material.depthTest = !onTop;
    this.lineMaterials.push(material);
    const lines = new LineSegments2(geometry, material);
    if (onTop) lines.renderOrder = 1; // draw after the mesh
    this.content.add(lines);
  }

  private frame(): void {
    const box = new THREE.Box3().setFromObject(this.content);
    const sphere = box.getBoundingSphere(new THREE.Sphere());
    const r = Math.max(sphere.radius, 1e-6);
    const dist = r / Math.sin(THREE.MathUtils.degToRad(this.camera.fov / 2));
    this.controls.target.copy(sphere.center);
    this.camera.position.copy(sphere.center).add(new THREE.Vector3(0.35, -1, 0.75).normalize().multiplyScalar(dist));
    this.camera.near = dist / 100;
    this.camera.far = dist * 100;
    this.camera.updateProjectionMatrix();
  }

  private clear(): void {
    for (const obj of [...this.content.children]) {
      this.content.remove(obj);
      obj.traverse((o) => {
        const any = o as THREE.Mesh;
        any.geometry?.dispose(); // GPU buffers are not garbage-collected either
        const mat = any.material as THREE.Material | THREE.Material[] | undefined;
        if (Array.isArray(mat)) mat.forEach((m) => m.dispose());
        else mat?.dispose(); // (disposing a material does not dispose its shared map texture)
      });
    }
    this.lineMaterials = [];
    this.surface = null;
  }

  private resize(): void {
    const w = this.container.clientWidth, h = this.container.clientHeight;
    this.renderer.setSize(w, h, false);
    this.renderer.domElement.style.width = "100%";
    this.renderer.domElement.style.height = "100%";
    this.camera.aspect = w / Math.max(h, 1);
    this.camera.updateProjectionMatrix();
    for (const m of this.lineMaterials) m.resolution.set(w, h);
  }
}
