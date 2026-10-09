# Case study: training-free tooth-gingiva segmentation of intraoral scans

**Result.** On 299 held-out Teeth3DS scans, the detected tooth-gingiva boundary is **0.285 mm** from the labelled one on
average (ASSD; median 0.219 mm), boundary F1 **0.906** at 0.5 mm, tooth IoU **0.936**. No neural network: geodesic
distances, curvature, and graph cuts, with a handful of parameters tuned on training scans. The first working
version measured 0.551 mm and F1 0.815. This page is the path between the two; every number links to an entry in
[DECISIONS.md](../DECISIONS.md).

## Setup and protocol

- **Data.** Teeth3DS+ (CC BY-NC-ND 4.0), used locally, never redistributed. Part 1 (300 training scans) split by
  index into sweep, classifier-training and validation folds. Test parts 5 and 6 (the official test split, 300 + 300).
- **Metrics.** Boundary-to-boundary distances between edge midpoints of the predicted and labelled tooth/gingiva
  boundaries: ASSD, HD95, boundary F1 at 0.25 and 0.5 mm; area-weighted tooth IoU (D72).
- **Discipline.** Tune on sweep scans, confirm on validation, then test. Every test run is counted and disclosed.
  Changes are judged by paired per-scan statistics (mean, SE, t, scans better/worse, size of the worst loss), not
  only by the mean. Part 6 was held back until the end and evaluated once per final method (D81, D85).

## The path

| Step | What changed | Test ASSD | F1 @ 0.5 |
|---|---|---|---|
| Baseline | Plane cut at a height quantile | 1.663 mm | 0.208 |
| D73 | Two-label geodesic Voronoi from cusp and gingiva seeds, edge cost weighted by concavity | 0.551 | 0.815 |
| D77 | No height-band tooth seeds | 0.520 | 0.819 |
| D78 | Graph cut: same seeds, boundary cost cheap along creases (max-flow written from scratch) | 0.338 | 0.883 |
| D85 | Per-tooth multi-label cut (cusp grouping + alpha-expansion) | **0.285** (part 6) | **0.906** |

### 1. First arrival leaks (D74-D77)
The error was a tail: median ~0.09 mm per scan, but some scans far worse. Decomposing false boundary by distance to
the nearest true tooth showed fake tooth regions on flat gingiva, correlated +0.88 with per-scan error. A sliver
hypothesis was rejected (same triangle quality on correct and false edges). The mechanism: Voronoi is first arrival,
so a tooth front that slips through one gap in a crease floods flat gingiva. Removing height-band tooth seeds helped
modestly; it treated the symptom.

### 2. A cut that pays for its boundary (D78)
A minimum s-t cut over the energy *area-weighted distance evidence + μ × (boundary length / (1 + β · crease
strength))* cannot flood cheaply: a leaked region pays for its whole perimeter on flat gingiva. Two sweeps ended at
the edge of the grid and were widened until the optimum was interior. Validation predicted the test gain
(−0.174 vs −0.182 mm). The first sweep looked catastrophic; the cause was the evaluation tool computing the crease
field only for the old method, which a consistency check in test mode would have caught.

### 3. Three hypotheses, measured and rejected (D80-D84)
F1 was recall-limited. The largest loss: interdental papillae labelled tooth, bridging two crowns.
- *The cut's length term shortcuts the papilla*: rejected, Voronoi (no length term) bridges even more.
- *The scanner cannot see between teeth, so there is no crease*: rejected, the missed boundary sits on the strongest
  creases of the arch.
- *The cut follows a competing crease*: rejected, almost no predicted boundary lies elsewhere.
- Looking at the worst cases showed the answer: the labels keep a gingiva strip ~0.17 mm wide (median) between
  touching crowns. A two-label model cannot express "tooth A | tooth B"; merging two crowns costs nothing. A
  representation limit, which is why no tuning helped.

### 4. Labels per tooth (D85)
Cusp tips are grouped into teeth (connected through vertices without deep creases; 86% of teeth come out clean
against the FDI labels), then a multi-label Potts energy is minimized by alpha-expansion, each move an exact cut.
Two surprises, both measured: carving the thin strip the labels contain made results worse (the multi-label cut
alone helped), and slightly over-splitting teeth beat the "correct" grouping. Profiling then made it 2.6× faster.

## What did not help, and why it is reported
Seed discs, island removal, strip carving, a geodesic star-shape prior, and a learned per-vertex data term were
each implemented, tested, measured, and rejected (D80, D85, D87, D88). The learned term's ROC AUC was below that of
one of its own inputs, the distance evidence the cut already uses, which predicted the null result.

## How it compares
Supervised networks trained on ~1,400 labelled scans report gingiva IoU around 0.96 per point (other splits); this
method reaches 0.911 per vertex (D79). Per-vertex IoU also understates boundary methods, since scans are densest at
teeth and margins (D80). A head-to-head with a published checkpoint is prepared on the 84 test scans that neither
method has seen (D81); a split audit found the checkpoint had been trained on 216 of the 300 "test" scans of part 6.

## Engineering
C++20 core with no I/O, compiled natively and to WebAssembly; Three.js viewer. Max-flow, alpha-expansion, BVH, the
heat method and the sparse factorization are written from scratch and checked against brute force or closed forms,
with mutation testing (which found gaps in the tests and once a bug in the mutation harness itself).

## Limits
"Margin" here is the tooth-gingiva boundary on unprepared arches, not the finish line of a crown preparation.
HD95 is still ~2 mm: a tail remains, mostly fake regions around wrong seeds and missed teeth.
