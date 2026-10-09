// Scripted README recording (D97): drives the built viewer in headless Chromium through an OBJ upload, rendering,
// overlays, margin detection scored against known labels, and naive vs signed-distance offset shells with the
// thickness map, then encodes the screenshots as a GIF (gifenc, inside the browser: no ffmpeg, no PNG decoder).
//
//   npm run build && node scripts/record-demo.mjs <mesh.obj> <labels.json> [out.gif]
//
// The mesh is the synthetic molar from tools/export_demo (licensed scans must not appear in public media).
import { spawn } from "node:child_process";
import { readFileSync, writeFileSync } from "node:fs";
import { createRequire } from "node:module";
import { resolve } from "node:path";
import { chromium } from "playwright-core";

const [objPath, labelsPath, outPath = "../docs/demo.gif"] = process.argv.slice(2);
if (!objPath || !labelsPath) {
  console.error("usage: node scripts/record-demo.mjs <mesh.obj> <labels.json> [out.gif]");
  process.exit(2);
}
const WIDTH = 1120, HEIGHT = 700, FPS = 4, SCALE = 0.7; // captured size; the GIF is scaled by SCALE
const require = createRequire(import.meta.url);
const gifencSource = readFileSync(require.resolve("gifenc").replace("gifenc.js", "gifenc.esm.js"), "utf8");

// Serve the production build (vite preview) and wait until it answers.
const server = spawn("npx", ["vite", "preview", "--port", "4173", "--strictPort"], { stdio: "ignore" });
const url = "http://localhost:4173/";
for (let i = 0; i < 100; ++i) {
  try {
    if ((await fetch(url)).ok) break;
  } catch {}
  await new Promise((r) => setTimeout(r, 200));
}

// SwiftShader gives headless Chromium a software WebGL context.
const browser = await chromium.launch({ args: ["--use-angle=swiftshader", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist"] });
const page = await browser.newPage({ viewport: { width: WIDTH, height: HEIGHT }, deviceScaleFactor: 1 });
page.on("pageerror", (e) => console.error("page error:", e.message));
const frames = []; // { jpeg: base64, delay: ms }

async function hold(ms, fps = FPS) {
  const n = Math.max(1, Math.round((ms / 1000) * fps));
  for (let i = 0; i < n; ++i) {
    frames.push({ jpeg: (await page.screenshot({ type: "jpeg", quality: 88 })).toString("base64"), delay: Math.round(1000 / fps) });
  }
}
async function caption(text) {
  console.log("step:", text.slice(0, 60));
  // A small caption bar over the 3D view, so the GIF explains itself without the README.
  await page.evaluate((t) => {
    let el = document.getElementById("demo-caption");
    if (!el) {
      el = document.createElement("div");
      el.id = "demo-caption";
      Object.assign(el.style, {
        position: "fixed", right: "16px", bottom: "16px", maxWidth: "560px", padding: "8px 12px", borderRadius: "6px",
        background: "rgba(20,22,28,0.82)", color: "#f2f2f2", font: "500 15px/1.35 system-ui, sans-serif", zIndex: 1000,
      });
      document.body.appendChild(el);
    }
    el.textContent = t;
  }, text);
}
// Text the viewer displays, so captions quote the measured numbers rather than claims.
const pageText = (sel) => page.evaluate((q) => document.querySelector(q)?.textContent ?? "", sel);
async function legendLine() {
  return (await page.evaluate(() => document.getElementById("legend")?.innerText ?? "")).split("\n").find((l) => l.startsWith("Thinnest")) ?? "";
}
async function orbit(dx, steps = 24, dy = 0) {
  const box = await page.locator("#view canvas").boundingBox();
  const cx = box.x + box.width / 2, cy = box.y + box.height / 2;
  await page.mouse.move(cx, cy);
  await page.mouse.down();
  for (let i = 1; i <= steps; ++i) {
    await page.mouse.move(cx + (dx * i) / steps, cy + ((dy || -dx * 0.15) * i) / steps);
    if (i % 3 === 0) await hold(160);
  }
  await page.mouse.up();
}
const overlay = (value) => page.locator(`input[name="overlay"][value="${value}"]`).check();
const click = (id) => page.locator(`#${id}`).click();

process.on("unhandledRejection", async (e) => {
  console.error("failed:", e.message.split("\n")[0]);
  console.error("stats:", (await page.evaluate(() => document.querySelector("#stats")?.textContent ?? "")).replace(/\s+/g, " ").slice(0, 200));
  console.error("error box:", await page.evaluate(() => document.getElementById("error")?.textContent));
  await page.screenshot({ path: "record-failure.png" });
  process.exit(1);
});
await page.goto(url);
await page.waitForFunction(() => document.querySelector("#stats")?.textContent.includes("Vertices"));
await caption("Mesh Inspection Workbench: C++ geometry core compiled to WebAssembly");
await hold(1500);

await caption("Open an OBJ (a synthetic molar, generated in code): parsed, welded and analysed in the browser");
await page.locator("#file").setInputFiles(resolve(objPath));
await page.waitForFunction(() => document.querySelector("#stats")?.textContent.replace(/[,\s\u202f]/g, "").includes("19881"));
const loadMs = (await pageText("#stats")).match(/Load \+ analysis \(this browser\)\s*([\d.]+)\s*ms/)?.[1];
await caption(`OBJ opened: 19,881 vertices, 39,200 faces; parsed, welded and analysed in the browser in ${loadMs} ms`);
await hold(2000);
await orbit(220);

await caption("Minimum principal curvature: the concave crease at the crown base is the tooth-gingiva margin");
await overlay("kmin");
await hold(2200);

await caption("Margin detection (graph cut), scored against the known boundary");
await overlay("shaded");
await click("detect-margin");
await page.locator("#labels").setInputFiles(resolve(labelsPath));
await page.waitForFunction(() => document.body.textContent.includes("vs labels"));
const scored = (await pageText("body")).match(/vs labels: ASSD ([\d.]+) mm, HD95 ([\d.]+) mm, boundary F1@0\.5 mm ([\d.]+)/);
await caption(`Margin detection (graph cut) vs the known boundary: ASSD ${scored?.[1]} mm, boundary F1 ${scored?.[3]} at 0.5 mm`);
await hold(3000);

await caption("Offset shell, 1 mm wall, by moving vertices along normals; thickness map (red = thinner than 1 mm)");
await page.locator("#wall").fill("1");
await click("shell-naive");
await click("thickness");
await page.waitForFunction(() => document.getElementById("legend")?.textContent.includes("Thinnest"));
await orbit(0, 18, 170); // look down onto the occlusal table and the fissure
// Quote only the thinnest wall: it comes from the fold under the fissure (checked on the occlusal table alone in
// tests/test_offset.cpp); the share below 95% also counts vertices at the patch rim, where cone rays reach the
// stitched side wall.
const naiveMin = (await legendLine()).match(/Thinnest wall ([^,]+)/)?.[1];
await caption(`Naive offset: thinnest wall ${naiveMin} mm, the inner copy folds through itself under the fissure`);
await hold(3200);

await caption("Same 1 mm wall from the signed distance field (signed heat method)");
await page.locator("#file").setInputFiles(resolve(objPath)); // back to the original surface
await page.waitForFunction(() => !document.getElementById("legend")?.textContent.includes("Thinnest"));
await page.locator("#wall").fill("1");
await click("shell-sdf");
await click("thickness");
await page.waitForFunction(() => document.getElementById("legend")?.textContent.includes("Thinnest"));
await orbit(0, 18, 170);
const sdf = (await legendLine()).match(/Thinnest wall ([^,]+), median ([^;]+)/);
await caption(`Signed-distance offset: thinnest wall ${sdf?.[1]} mm, median ${sdf?.[2]} mm: the nominal 1 mm within grid error`);
await hold(3000);
await orbit(-160, 18, 40);
await hold(1200);

// Encode in the page: decode JPEGs with createImageBitmap, one shared 256-colour palette, gifenc frames.
console.log(`${frames.length} frames; encoding...`);
if (process.env.DEMO_FRAMES_DIR) {  // optional: a few frames as JPEGs, for checking the recording by eye
  for (const f of [0.1, 0.3, 0.5, 0.7, 0.95]) {
    const i = Math.floor(f * frames.length);
    writeFileSync(`${process.env.DEMO_FRAMES_DIR}/frame-${String(i).padStart(3, "0")}.jpg`, Buffer.from(frames[i].jpeg, "base64"));
  }
}
const encoder = await browser.newPage();
const gifBase64 = await encoder.evaluate(
  async ({ src, frames, w, h, scale }) => {
    const { GIFEncoder, quantize, applyPalette } = await import("data:text/javascript;base64," + btoa(unescape(encodeURIComponent(src))));
    const ow = Math.round(w * scale), oh = Math.round(h * scale);
    const canvas = new OffscreenCanvas(ow, oh);
    const ctx = canvas.getContext("2d", { willReadFrequently: true });
    ctx.imageSmoothingQuality = "high";
    const rgba = async (b64) => {
      const bytes = Uint8Array.from(atob(b64), (c) => c.charCodeAt(0));
      ctx.drawImage(await createImageBitmap(new Blob([bytes], { type: "image/jpeg" })), 0, 0, ow, oh);
      return ctx.getImageData(0, 0, ow, oh).data;
    };
    // Palette from a sample of frames spread over the recording (consistent colours, smaller file).
    const sample = [];
    for (const i of [0, 0.25, 0.5, 0.75, 0.999].map((f) => Math.floor(f * frames.length))) sample.push(await rgba(frames[i].jpeg));
    const joined = new Uint8ClampedArray(sample.reduce((n, a) => n + a.length, 0));
    sample.reduce((off, a) => (joined.set(a, off), off + a.length), 0);
    const palette = quantize(joined, 256);
    const gif = GIFEncoder();
    for (const f of frames) gif.writeFrame(applyPalette(await rgba(f.jpeg), palette), ow, oh, { palette, delay: f.delay });
    gif.finish();
    const out = gif.bytes();
    let s = "";
    for (let i = 0; i < out.length; i += 0x8000) s += String.fromCharCode(...out.subarray(i, i + 0x8000));
    return btoa(s);
  },
  { src: gifencSource, frames, w: WIDTH, h: HEIGHT, scale: SCALE },
);
writeFileSync(outPath, Buffer.from(gifBase64, "base64"));
console.log(`wrote ${outPath} (${(Buffer.byteLength(gifBase64, "base64") / 1e6).toFixed(1)} MB)`);
await browser.close();
server.kill();
