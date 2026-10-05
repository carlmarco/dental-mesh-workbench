export type RGB = [number, number, number];

export const NO_DATA: RGB = [0.5, 0.5, 0.52]; // boundary / non-manifold: no pointwise curvature

// Diverging map for signed quantities: blue (negative), near-white (0), red (positive).
// Endpoints follow Moreland's "cool to warm" map; linear interpolation through white.
const COOL: RGB = [0.23, 0.299, 0.754];
const MID: RGB = [0.92, 0.92, 0.92];
const WARM: RGB = [0.706, 0.016, 0.15];

export function diverging(t: number): RGB {
  const s = Math.max(-1, Math.min(1, t));
  const [a, b, u] = s < 0 ? [MID, COOL, -s] : [MID, WARM, s];
  return [a[0] + (b[0] - a[0]) * u, a[1] + (b[1] - a[1]) * u, a[2] + (b[2] - a[2]) * u];
}

// Symmetric color range [-R, R] with R = the q-quantile of |value| over finite values.
// Curvature has heavy-tailed outliers (slivers, creases); a max-based range would wash
// everything else out to white.
export function robustRange(values: Float32Array, q = 0.98): number {
  const mags: number[] = [];
  for (const v of values) if (Number.isFinite(v)) mags.push(Math.abs(v));
  if (mags.length === 0) return 1;
  mags.sort((x, y) => x - y);
  const r = mags[Math.min(mags.length - 1, Math.floor(q * mags.length))];
  return r > 0 ? r : 1;
}

// Distinct colors for component ids.
const PALETTE: RGB[] = [
  [0.27, 0.51, 0.71],
  [0.87, 0.52, 0.16],
  [0.33, 0.66, 0.41],
  [0.77, 0.31, 0.32],
  [0.51, 0.45, 0.7],
  [0.58, 0.46, 0.33],
  [0.85, 0.52, 0.71],
  [0.55, 0.55, 0.3],
];
export function categorical(i: number): RGB {
  return PALETTE[i % PALETTE.length];
}
