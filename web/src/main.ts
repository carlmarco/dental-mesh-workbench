import "./style.css";

import { categorical } from "./colormap.ts";
import { matchPoints, parseCuspLandmarks } from "./cusps.ts";
import { Dmw, type MeshData, type Preset, PRESETS, type Stats } from "./dmw.ts";
import { ISOLINES, LAYER_COLORS, type Layers, type Overlay, Viewer } from "./viewer.ts";

const LABELS: Record<Preset, string> = {
  plate_handle: "Plate with a handle (genus 1)",
  grid_holes: "Plate with two holes (genus 0)",
  defects: "Defect showcase",
  brick: "Obtuse \"brick\" grid (cotan fails)",
  icosphere: "Icosphere",
  torus: "Torus",
  cylinder: "Open cylinder",
  grid: "Flat grid",
  mobius: "Möbius strip",
};

const byId = <T extends HTMLElement>(id: string) => document.getElementById(id) as T;
const preset = byId<HTMLSelectElement>("preset");
const fileInput = byId<HTMLInputElement>("file");
const weld = byId<HTMLInputElement>("weld");
const errorBox = byId<HTMLParagraphElement>("error");
const view = byId<HTMLElement>("view");
const idt = byId<HTMLInputElement>("idt");

const dmw = await Dmw.create();
const viewer = new Viewer(view);
let data: MeshData | null = null;
let geodesicInfo = ""; // legend text for the current distance field
let geodesicSource = 0;
let cuspXyz: Float32Array | null = null;
let truthXyz: Float32Array | null = null;
let cuspInfo = "";
let cuspMillis: number | null = null;
let marginInfo = "";

const hex = (c: number) => `#${c.toString(16).padStart(6, "0")}`;
const rgb = ([r, g, b]: number[]) => `rgb(${r * 255}, ${g * 255}, ${b * 255})`;
const esc = (s: string) => s.replace(/[&<>"]/g, (ch) => `&#${ch.charCodeAt(0)};`);
const fmt = (x: number) => (Math.abs(x) >= 1e4 || (x !== 0 && Math.abs(x) < 1e-3) ? x.toExponential(2) : x.toFixed(3));

for (const el of document.querySelectorAll<HTMLElement>("i[data-swatch]")) {
  el.style.background = hex(LAYER_COLORS[el.dataset.swatch as keyof typeof LAYER_COLORS]);
}
for (const p of PRESETS) preset.add(new Option(LABELS[p], p));

function overlay(): Overlay {
  return document.querySelector<HTMLInputElement>('input[name="overlay"]:checked')!.value as Overlay;
}
function layers(): Layers {
  const on = (v: string) => document.querySelector<HTMLInputElement>(`#layers input[value="${v}"]`)!.checked;
  return { boundary: on("boundary"), nonmanifold: on("nonmanifold"), misoriented: on("misoriented"), vertices: on("vertices"), handles: on("handles"), cusps: on("cusps"), margin: on("margin"), wireframe: on("wireframe") };
}

// Heat-method distance from `source`; the first query on a mesh also factors the solver.
function computeGeodesic(source: number): boolean {
  try {
    const first = geodesicInfo === "";
    const g = dmw.geodesic(source);
    geodesicSource = source;
    viewer.setGeodesic(g.distance, source);
    geodesicInfo = `Source: vertex ${source}. ${first ? "Factor + solve" : "Solve (reusing factorization)"}: ${g.millis.toFixed(1)} ms.`;
    return true;
  } catch (e) {
    errorBox.textContent = (e as Error).message;
    return false;
  }
}

function display(d: MeshData, keepCamera: boolean): void {
  if (d !== data) {
    viewer.setGeodesic(null, null); // a new mesh invalidates the distance field
    geodesicInfo = "";
    if (!keepCamera) {
      cuspXyz = truthXyz = null;
      cuspInfo = "";
      cuspMillis = null;
      marginInfo = "";
      undercutInfo = "";
      viewer.setUndercut(null);
      viewer.setMargin(null, null);
      viewer.setCusps(null, null);
    }
  }
  data = d;
  errorBox.textContent = "";
  if (overlay() === "geodesic" && geodesicInfo === "" && !computeGeodesic(0)) {
    document.querySelector<HTMLInputElement>('input[name="overlay"][value="shaded"]')!.checked = true;
  }
  const range = keepCamera ? viewer.restyle(overlay(), layers()) : viewer.show(d, overlay(), layers());
  renderLegend(range);
  renderStats(d.stats, d.millis);
}

function run(load: () => MeshData): void {
  try {
    display(load(), false);
  } catch (e) {
    errorBox.textContent = `Could not load: ${(e as Error).message}`;
  }
}

function renderLegend(range: number | null): void {
  const legend = byId<HTMLDivElement>("legend");
  if (range === null) {
    legend.innerHTML = "";
    return;
  }
  if (overlay() === "geodesic") {
    legend.innerHTML = `<div class="ramp seq"></div>
      <div class="ramp-labels"><span>0</span><span>${fmt(range)}</span></div>
      <div>${ISOLINES} isolines, spacing ${fmt(range / ISOLINES)}. ${esc(geodesicInfo)} Click the surface to move the source.</div>`;
    return;
  }
  legend.innerHTML = `<div class="ramp"></div>
    <div class="ramp-labels"><span>${fmt(-range)}</span><span>0</span><span>${fmt(range)}</span></div>
    <div>Range: 98th percentile of |value|. Gray: no pointwise value (boundary).</div>`;
}

function renderStats(s: Stats, millis: number): void {
  const flags = (c: Stats["components"][number]) =>
    [!c.manifold && "non-manifold", c.manifold && !c.orientable && "non-orientable", !c.consistent && c.orientable && "mis-wound"]
      .filter(Boolean)
      .map((f) => `<span class="flag">${f}</span>`)
      .join(" ") || '<span class="ok">ok</span>';
  const rows = s.components
    .map(
      (c, i) => `<tr><td><span class="dot" style="background:${rgb(categorical(i))}"></span>${i}</td>
        <td>${c.V}</td><td>${c.F}</td><td>${c.chi}</td><td>${c.b}</td>
        <td>${c.genus ?? "–"}</td><td>${c.betti ? c.betti.join(", ") : "–"}</td></tr>
        <tr><td></td><td colspan="6" style="text-align:left">${flags(c)}</td></tr>`,
    )
    .join("");
  const totalChi = s.components.reduce((a, c) => a + c.chi, 0);
  const defect = (label: string, n: number, color?: number) =>
    `<dt>${color !== undefined ? `<span class="dot" style="background:${hex(color)}"></span>` : ""}${label}</dt><dd class="${n ? "flag" : ""}">${n}</dd>`;
  const curv = s.curvature
    ? `<dl>
        <dt>Σ angle defect / 2π</dt><dd>${fmt(s.curvature.totalAngleDefect / (2 * Math.PI))}</dd>
        <dt>Σ χ (Gauss-Bonnet target)</dt><dd>${totalChi}</dd>
        <dt>H range</dt><dd>${s.curvature.meanRange ? s.curvature.meanRange.map(fmt).join(" … ") : "–"}</dd>
        <dt>K range</dt><dd>${s.curvature.gaussianRange ? s.curvature.gaussianRange.map(fmt).join(" … ") : "–"}</dd>
        <dt>Zero-area faces</dt><dd>${s.curvature.degenerateFaces}</dd>
        <dt>Laplacian</dt><dd>${s.intrinsicDelaunay ? `intrinsic Delaunay (${s.flips} flips)` : "cotan"}</dd>
      </dl>`
    : `<p class="flag">Not computed: the half-edge build rejected this mesh (${esc(s.halfedge.error)}).
        The curvature code needs a manifold with consistent orientation. Gray in the curvature views means
        "no data", not zero.</p>`;

  byId<HTMLDivElement>("stats").innerHTML = `
    <h2>Mesh</h2>
    <dl><dt>Vertices</dt><dd>${s.vertices}</dd><dt>Edges</dt><dd>${s.edges}</dd><dt>Faces</dt><dd>${s.faces}</dd>
        <dt>Load + analysis (this browser)</dt><dd>${millis.toFixed(1)} ms</dd></dl>
    <h2>Components</h2>
    <table><tr><th>#</th><th>V</th><th>F</th><th>χ</th><th>b</th><th>g</th><th>β</th></tr>${rows}</table>
    <h2>Defects</h2>
    <dl>
      ${defect("Boundary edges", s.edgeKinds.boundary, LAYER_COLORS.boundary)}
      ${defect("Non-manifold edges", s.edgeKinds.nonmanifold, LAYER_COLORS.nonmanifold)}
      ${defect("Misoriented edges", s.edgeKinds.misoriented, LAYER_COLORS.misoriented)}
      ${defect("Non-manifold vertices", s.nonmanifoldVertices, LAYER_COLORS.vertices)}
      ${defect("Handle loops (2 per handle)", s.handleLoops.count, LAYER_COLORS.handles)}
      ${defect("Invalid faces", s.invalidFaces)}
      ${defect("Duplicate faces", s.duplicateFaces)}
      ${defect("Isolated vertices", s.isolatedVertices)}
    </dl>
    ${s.handleLoops.count ? `<p class="hint">Handle loops: shortest ${s.handleLoops.lengths.slice(0, 5).map(fmt).join(", ")}${s.handleLoops.count > 5 ? ", …" : ""} (mesh units; upper bounds on handle size).</p>` : ""}
    ${s.excludedFaces ? `<p class="hint">${s.excludedFaces} faces at defects excluded from the half-edge analysis view (D69).</p>` : ""}
    ${cuspInfo ? `<h2>Cusps</h2><p>${cuspInfo}</p>` : ""}
    ${marginInfo ? `<h2>Margin</h2><p>${marginInfo}</p>` : ""}
    ${undercutInfo ? `<h2>Undercut</h2><p>${undercutInfo}</p>` : ""}
    <h2>Curvature</h2>${curv}`;
}

async function openFile(file: File): Promise<void> {
  idt.checked = false;
  const ext = file.name.toLowerCase().split(".").pop();
  if (ext !== "stl" && ext !== "obj") {
    errorBox.textContent = "Only .stl and .obj files are supported.";
    return;
  }
  const bytes = new Uint8Array(await file.arrayBuffer());
  run(() => dmw.loadFile(bytes, ext, Math.max(0, Number(weld.value) || 0)));
}

preset.addEventListener("change", () => {
  idt.checked = false; // each new mesh starts with the plain cotan Laplacian
  run(() => dmw.generate(preset.value as Preset));
});
idt.addEventListener("change", () => {
  if (!data) return;
  const keepSource = geodesicInfo !== "" ? geodesicSource : null;
  try {
    const d = dmw.setIntrinsicDelaunay(idt.checked);
    viewer.setGeodesic(null, null);
    geodesicInfo = "";
    data = d;
    if (overlay() === "geodesic" && keepSource !== null) computeGeodesic(keepSource);
    display(d, true);
  } catch (e) {
    errorBox.textContent = (e as Error).message;
  }
});
// A real <button> is focusable and keyboard-operable; it forwards to the hidden file input.
byId<HTMLButtonElement>("open").addEventListener("click", () => fileInput.click());
fileInput.addEventListener("change", () => fileInput.files?.[0] && openFile(fileInput.files[0]));
view.addEventListener("dragover", (e) => e.preventDefault());
view.addEventListener("drop", (e) => {
  e.preventDefault();
  const f = e.dataTransfer?.files[0];
  if (f) void openFile(f);
});
function updateCuspInfo(): void {
  const parts: string[] = [];
  if (cuspXyz) parts.push(`${cuspXyz.length / 3} tips detected${cuspMillis !== null ? ` in ${cuspMillis.toFixed(0)} ms` : ""}`);
  if (truthXyz) parts.push(`${truthXyz.length / 3} landmarks loaded`);
  if (cuspXyz && truthXyz) {
    const m = matchPoints(cuspXyz, truthXyz, 1.0);
    const p = m.tp / Math.max(m.detections, 1), r = m.tp / Math.max(m.truth, 1);
    parts.push(`at 1 mm: precision ${p.toFixed(2)}, recall ${r.toFixed(2)}, F1 ${(p + r > 0 ? (2 * p * r) / (p + r) : 0).toFixed(2)}, median error ${m.medianError.toFixed(2)} mm`);
  }
  cuspInfo = parts.join(" · ");
  viewer.setCusps(cuspXyz, truthXyz);
  if (data) display(data, true);
}

byId<HTMLButtonElement>("detect").addEventListener("click", () => {
  if (!data) return;
  try {
    const c = dmw.detectCusps();
    cuspXyz = new Float32Array(c.vertices.length * 3);
    c.vertices.forEach((v, i) => cuspXyz!.set(data!.positions.subarray(3 * v, 3 * v + 3), 3 * i));
    cuspMillis = c.millis;
    updateCuspInfo();
  } catch (e) {
    errorBox.textContent = (e as Error).message;
  }
});
const landmarkInput = byId<HTMLInputElement>("landmarks");
byId<HTMLButtonElement>("open-landmarks").addEventListener("click", () => landmarkInput.click());
landmarkInput.addEventListener("change", async () => {
  const f = landmarkInput.files?.[0];
  if (!f) return;
  try {
    truthXyz = parseCuspLandmarks(await f.text());
    updateCuspInfo();
  } catch (e) {
    errorBox.textContent = `Could not read landmarks: ${(e as Error).message}`;
  }
});

let marginEdges: Uint32Array | null = null;
byId<HTMLButtonElement>("detect-margin").addEventListener("click", () => {
  if (!data) return;
  try {
    const r = dmw.detectMargin();
    marginEdges = r.edges;
    marginInfo = `${r.edges.length / 2} margin edges in ${r.millis.toFixed(0)} ms`;
    viewer.setMargin(marginEdges, null);
    display(data, true);
  } catch (e) {
    errorBox.textContent = (e as Error).message;
  }
});
const labelInput = byId<HTMLInputElement>("labels");
byId<HTMLButtonElement>("open-labels").addEventListener("click", () => labelInput.click());
labelInput.addEventListener("change", async () => {
  const f = labelInput.files?.[0];
  if (!f || !data) return;
  try {
    if (!marginEdges) throw new Error("detect the margin first");
    // Teeth3DS labels: FDI tooth number per vertex, 0 = gingiva.
    const labels = (JSON.parse(await f.text()) as { labels: number[] }).labels;
    const flags = Uint8Array.from(labels, (l) => (l !== 0 ? 1 : 0));
    const { metrics: m, truthEdges } = dmw.compareMargin(flags);
    marginInfo = `${marginEdges.length / 2} margin edges · vs labels: ASSD ${m.assd.toFixed(3)} mm, HD95 ${m.hd95.toFixed(2)} mm, boundary F1@0.5 mm ${m.f1_050.toFixed(2)}, tooth IoU ${m.iou.toFixed(3)}`;
    viewer.setMargin(marginEdges, truthEdges);
    display(data, true);
  } catch (e) {
    errorBox.textContent = `Could not compare: ${(e as Error).message}`;
  }
});

let undercutInfo = "";
function runUndercut(best: boolean): void {
  if (!data) return;
  try {
    const u = dmw.computeUndercut(best);
    viewer.setUndercut(u.faces);
    undercutInfo = `${best ? "Best path of insertion" : "Occlusal axis"}: ${(100 * u.fraction).toFixed(1)}% of the ${u.teethOnly ? "detected teeth" : "mesh"} undercut, tilt ${u.tilt.toFixed(1)}° from the occlusal axis (${u.millis.toFixed(0)} ms)`;
    document.querySelector<HTMLInputElement>('input[name="overlay"][value="undercut"]')!.checked = true;
    display(data, true);
  } catch (e) {
    errorBox.textContent = (e as Error).message;
  }
}
byId<HTMLButtonElement>("undercut").addEventListener("click", () => runUndercut(false));
byId<HTMLButtonElement>("best-axis").addEventListener("click", () => runUndercut(true));

viewer.onPick = (v) => {
  if (overlay() !== "geodesic" || !data) return;
  if (computeGeodesic(v)) display(data, true);
};
for (const input of document.querySelectorAll<HTMLInputElement>("#overlay input, #layers input")) {
  input.addEventListener("change", () => data && display(data, true));
}
const applyTheme = () => viewer.setBackground(getComputedStyle(document.documentElement).getPropertyValue("--view-bg").trim());
window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", applyTheme);
applyTheme();

// Local inspection (dev server only, D84): ?scan=<abs path .obj>&labels=<abs path .json>&focus=x,y,z[&dist=mm]
// loads a scan from the git-ignored data folder, detects the margin, overlays the labelled boundary and frames
// the given point. Vite serves files outside the web root under /@fs (allowed in vite.config.ts).
async function inspectFromUrl(): Promise<boolean> {
  if (!import.meta.env.DEV) return false; // local scans are served by the dev server only; never on the public site
  const q = new URLSearchParams(location.search);
  const scan = q.get("scan");
  if (!scan) return false;
  try {
    const bytes = new Uint8Array(await (await fetch(`/@fs${scan}`)).arrayBuffer());
    display(dmw.loadFile(bytes, "obj"), false);
    const r = dmw.detectMargin();
    marginEdges = r.edges;
    marginInfo = `${r.edges.length / 2} margin edges in ${r.millis.toFixed(0)} ms`;
    viewer.setMargin(marginEdges, null);
    const labelsPath = q.get("labels");
    if (labelsPath) {
      const labels = ((await (await fetch(`/@fs${labelsPath}`)).json()) as { labels: number[] }).labels;
      const { metrics: m, truthEdges } = dmw.compareMargin(Uint8Array.from(labels, (l) => (l !== 0 ? 1 : 0)));
      marginInfo += ` · vs labels: ASSD ${m.assd.toFixed(3)} mm, HD95 ${m.hd95.toFixed(2)} mm, F1@0.5 ${m.f1_050.toFixed(2)}, IoU ${m.iou.toFixed(3)}`;
      viewer.setMargin(marginEdges, truthEdges);
    }
    const ov = q.get("overlay");
    if (ov) document.querySelector<HTMLInputElement>(`input[name="overlay"][value="${ov}"]`)!.checked = true;
    display(data!, true);
    const f = q.get("focus")?.split(",").map(Number);
    if (f && f.length === 3) viewer.focus([f[0], f[1], f[2]], Number(q.get("dist") ?? 12));
    return true;
  } catch (e) {
    errorBox.textContent = `Could not inspect: ${(e as Error).message}`;
    return false;
  }
}

if (!(await inspectFromUrl())) run(() => dmw.generate("plate_handle"));
