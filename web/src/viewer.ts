import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";
import { LineMaterial } from "three/addons/lines/LineMaterial.js";
import { LineSegments2 } from "three/addons/lines/LineSegments2.js";
import { LineSegmentsGeometry } from "three/addons/lines/LineSegmentsGeometry.js";

import { categorical, diverging, NO_DATA, robustRange, type RGB } from "./colormap.ts";
import type { MeshData } from "./dmw.ts";

export type Overlay = "shaded" | "mean" | "gaussian" | "components";

export interface Layers {
  boundary: boolean;
  nonmanifold: boolean;
  misoriented: boolean;
  vertices: boolean;
  wireframe: boolean;
}

export const LAYER_COLORS = {
  boundary: 0x2f80ed,
  nonmanifold: 0xe0342b,
  misoriented: 0xf2a900,
  vertices: 0xc026d3,
} as const;

const EXCLUDED = 0xffffffff;
const SHADED: RGB = [0.72, 0.75, 0.8];

export class Viewer {
  private readonly renderer: THREE.WebGLRenderer;
  private readonly scene = new THREE.Scene();
  private readonly camera = new THREE.PerspectiveCamera(40, 1, 0.001, 1000);
  private readonly controls: OrbitControls;
  private readonly content = new THREE.Group();
  private readonly container: HTMLElement;
  private lineMaterials: LineMaterial[] = [];
  private data: MeshData | null = null;

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
    if (overlay === "components") {
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
      vertexColors: true,
      side: THREE.DoubleSide, // open meshes show their back faces
      roughness: 0.65,
      metalness: 0.0,
      polygonOffset: true, // push faces back so edge overlays don't z-fight
      polygonOffsetFactor: 1,
      polygonOffsetUnits: 1,
    });
    this.content.add(new THREE.Mesh(geometry, material));

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
        else mat?.dispose();
      });
    }
    this.lineMaterials = [];
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
