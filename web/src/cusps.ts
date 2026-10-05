// Cusp landmark loading and matching for the viewer (mirrors core/metrics: greedy one-to-one).
export type Points = Float32Array; // xyz triples

// 3DTeethLand __kpt.json: {"objects": [{"class": "Cusp", "coord": [x, y, z]}, ...]}
export function parseCuspLandmarks(text: string): Points {
  const json = JSON.parse(text) as { objects?: { class?: string; coord?: number[] }[] };
  const xyz: number[] = [];
  for (const o of json.objects ?? []) {
    if (o.class === "Cusp" && o.coord?.length === 3) xyz.push(...o.coord);
  }
  return Float32Array.from(xyz);
}

export interface Match {
  tp: number;
  detections: number;
  truth: number;
  medianError: number;
}

export function matchPoints(det: Points, gt: Points, tolerance: number): Match {
  const pairs: [number, number, number][] = [];
  for (let i = 0; i < det.length / 3; ++i) {
    for (let j = 0; j < gt.length / 3; ++j) {
      const d = Math.hypot(det[3 * i] - gt[3 * j], det[3 * i + 1] - gt[3 * j + 1], det[3 * i + 2] - gt[3 * j + 2]);
      if (d <= tolerance) pairs.push([d, i, j]);
    }
  }
  pairs.sort((a, b) => a[0] - b[0]);
  const usedD = new Set<number>(), usedG = new Set<number>(), dist: number[] = [];
  for (const [d, i, j] of pairs) {
    if (usedD.has(i) || usedG.has(j)) continue;
    usedD.add(i), usedG.add(j), dist.push(d);
  }
  dist.sort((a, b) => a - b);
  return { tp: dist.length, detections: det.length / 3, truth: gt.length / 3, medianError: dist.length ? dist[dist.length >> 1] : NaN };
}
