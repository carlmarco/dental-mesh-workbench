import "./style.css";

import { categorical } from "./colormap.ts";
import { Dmw, type MeshData, type Preset, PRESETS, type Stats } from "./dmw.ts";
import { LAYER_COLORS, type Layers, type Overlay, Viewer } from "./viewer.ts";

const LABELS: Record<Preset, string> = {
  plate_handle: "Plate with a handle (genus 1)",
  grid_holes: "Plate with two holes (genus 0)",
  defects: "Defect showcase",
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

const dmw = await Dmw.create();
const viewer = new Viewer(view);
let data: MeshData | null = null;

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
  return { boundary: on("boundary"), nonmanifold: on("nonmanifold"), misoriented: on("misoriented"), vertices: on("vertices"), wireframe: on("wireframe") };
}

function display(d: MeshData, keepCamera: boolean): void {
  data = d;
  errorBox.textContent = "";
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
      </dl>`
    : `<p class="flag">Unavailable: the half-edge build rejected this mesh (${esc(s.halfedge.error)}). Curvature needs a manifold.</p>`;

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
      ${defect("Invalid faces", s.invalidFaces)}
      ${defect("Duplicate faces", s.duplicateFaces)}
      ${defect("Isolated vertices", s.isolatedVertices)}
    </dl>
    <h2>Curvature</h2>${curv}`;
}

async function openFile(file: File): Promise<void> {
  const ext = file.name.toLowerCase().split(".").pop();
  if (ext !== "stl" && ext !== "obj") {
    errorBox.textContent = "Only .stl and .obj files are supported.";
    return;
  }
  const bytes = new Uint8Array(await file.arrayBuffer());
  run(() => dmw.loadFile(bytes, ext, Math.max(0, Number(weld.value) || 0)));
}

preset.addEventListener("change", () => run(() => dmw.generate(preset.value as Preset)));
fileInput.addEventListener("change", () => fileInput.files?.[0] && openFile(fileInput.files[0]));
view.addEventListener("dragover", (e) => e.preventDefault());
view.addEventListener("drop", (e) => {
  e.preventDefault();
  const f = e.dataTransfer?.files[0];
  if (f) void openFile(f);
});
for (const input of document.querySelectorAll<HTMLInputElement>("#overlay input, #layers input")) {
  input.addEventListener("change", () => data && display(data, true));
}
const applyTheme = () => viewer.setBackground(getComputedStyle(document.documentElement).getPropertyValue("--view-bg").trim());
window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", applyTheme);
applyTheme();

run(() => dmw.generate("plate_handle"));
