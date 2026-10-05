# Design Decisions

Each entry: the choice, the alternatives considered, and the reason.

## D1. C++20, extensions off
- **Choice:** `CMAKE_CXX_STANDARD 20`, `CMAKE_CXX_EXTENSIONS OFF`.
- **Alternatives:** C++17; GNU extensions on.
- **Reason:** C++20 gives `std::span` (safe views over mesh buffers) and concepts at no
  real cost; Emscripten's clang supports it. Extensions off keeps us on the standard.

## D2. Strict warnings as errors, on our targets only
- **Choice:** `-Wall -Wextra -Wpedantic -Wconversion -Wshadow -Werror` on `dmw_core`.
- **Alternatives:** global flags; no `-Werror`.
- **Reason:** geometry code is full of int/size_t/float/double conversions; `-Wconversion`
  catches silent precision and sign bugs. Scoped per target so third-party code is unaffected.

## D3. Catch2 v3.16.0 via FetchContent (pinned tag)
- **Choice:** Catch2, fetched at configure time at tag `v3.16.0`.
- **Alternatives:** GoogleTest; vendoring the source; Homebrew install.
- **Reason:** Catch2 syntax (`REQUIRE(a == b)`, `Approx`/matchers for floats) reads well for
  geometry tests. FetchContent + a pinned tag keeps the build reproducible without a
  system-installed dependency. Release date of v3.16.0: 2026-08-25.

## D4. Build with Make, not Ninja
- **Choice:** CMake's default Unix Makefiles generator.
- **Alternatives:** Ninja.
- **Reason:** one less tool to install; project is small. Revisit if build time becomes a problem.

## D5. Emscripten via emsdk, pinned to 6.0.11
- **Choice:** emsdk cloned to `~/emsdk`, version `6.0.11` (current `latest` in emsdk's release
  manifest as of 2026-10-04).
- **Alternatives:** Homebrew `emscripten` formula.
- **Reason:** emsdk lets us pin the exact toolchain for reproducible builds.
- **Status:** confirmed. 6.0.11 builds `dmw_core` plus the embind module under `-Werror`, and
  `add` is callable from Node (verified 2026-10-04).

## D6. Core is a pure C++ library
- **Choice:** `dmw_core` has no I/O, no Emscripten, no JS dependencies; WASM bindings will be a
  separate thin target.
- **Reason:** the same code builds natively (for tests/benchmarks) and to WASM, and is testable
  without a browser.

## D7. Web tooling: TypeScript 7.0.2, Vite 8.3.2, exact-pinned
- **Choice:** `typescript@7.0.2`, `vite@8.3.2`, `@types/node@26.6.4`, no `^` ranges; Node 26.10.0 / npm 11.19.1 for
  the web side (emsdk's bundled Node 24 is used only by Emscripten itself).
- **Alternatives:** webpack/esbuild-only; caret ranges; ts-node/tsx for the smoke test.
- **Reason:** Vite handles the `.wasm` asset and ES-module glue with no config. Exact pins for
  reproducibility. The smoke test runs `.ts` directly via Node's built-in type stripping, so no extra
  runner dependency (hence `erasableSyntaxOnly` in tsconfig).

## D8. Hand-written type declaration for the WASM module
- **Choice:** `web/src/dmw-module.d.ts` declares the embind surface by hand.
- **Alternatives:** `emcc --emit-tsd` (auto-generated typings); untyped `any`.
- **Reason:** simplest for one function. Known cost: it must be kept in sync with `bindings.cpp`.
  Revisit `--emit-tsd` when the API surface grows (Milestone 5).

## D9. Known issue: `node:module` warning in Vite build
- `-sENVIRONMENT=web,node` leaves a Node code path in the glue, so Vite warns that `node:module` was
  externalized for the browser. The build succeeds, and the page renders correctly in a browser
  (verified manually, M1).
  Candidate fix: build a web-only variant for Vite and keep the node variant for tests.

## D10. Half-edge builder rejects non-manifold input
- **Choice:** construction fails with an error if any directed edge occurs twice (an edge with 3+ faces, or
  inconsistent orientation). Topology diagnostics (M3) run on the raw indexed mesh instead.
- **Alternatives:** a tolerant half-edge structure that records non-manifold edges.
- **Reason:** keeps the half-edge invariants strict (`twin` is an involution) and the traversal code simple;
  broken meshes never need half-edge form to be diagnosed.

## D11. Division of work
- Loaders, welding and plumbing are written by Claude; half-edge, cotan weights, curvature and the heat
  method are written by the author.
- **Amendment (M2):** at the author's explicit request, Claude wrote `halfedge.cpp`, following the formal
  algorithm reviewed beforehand. Tests were verified to catch planted bugs (see D26).
- **Amendment (M3 onward):** the author chose to have Claude implement all remaining code, including the
  curvature and heat-method core, with every decision explained and recorded for later review.

## D12. `double` positions in the core
- **Alternatives:** `float`.
- **Reason:** cotangent weights and angle-defect sums amplify rounding error; `float` is used only at the
  GPU boundary.

## D13. `uint32_t` vertex indices
- **Alternatives:** `size_t`, `int`.
- **Reason:** half the memory of `size_t`, and maps directly to WebGL `Uint32Array` index buffers.

## D14. Parsers take text/bytes, not file paths
- **Reason:** follows from D6; the browser provides bytes, not paths. Tests use inline fixtures.

## D15. Errors as return values, not exceptions
- **Alternatives:** exceptions; `std::expected` (C++23).
- **Reason:** Emscripten's default `DISABLE_EXCEPTION_CATCHING=1` (checked in emsdk 6.0.11) makes a throw
  abort the module. A small `LoadResult { mesh; error; ok() }` works on both targets.

## D16. OBJ polygons are fan-triangulated
- **Reason:** simple and preserves winding.
- **Limitation:** correct only for convex polygons (or star-shaped from the first vertex). Ear clipping
  would be the general fix; out of scope.

## D18. `strtod` (not `std::from_chars`) for parsing floats
- **Finding:** `std::from_chars` for `double` compiles and is correct under Emscripten 6.0.11, but is a
  deleted overload in Apple's system libc++ (Apple clang 17). Integer `from_chars` works on both.
- **Choice:** `std::strtod` on a NUL-terminated copy of each token, rejecting partial parses and non-finite
  values.
- **Alternatives:** `from_chars` behind a feature check (two code paths that behave differently between
  native and WASM); a third-party parser such as fast_float (extra dependency).
- **Caveat:** `strtod` honors the C locale's decimal separator. We never call `setlocale`, so it stays `"C"`
  (decimal point `.`). An application embedding this core that changes the locale could break parsing.

## D19. STL: detect binary by size; ignore stored normals
- **Choice:** binary iff `size == 84 + 50*N` (N = uint32 at offset 80, computed in 64-bit); otherwise parse as
  ASCII. Facet normals are read past but not used.
- **Alternatives:** sniff for a leading `"solid"`; trust stored normals.
- **Reason:** many binary exporters put `"solid"` in the 80-byte header, so prefix sniffing misclassifies them.
  Stored normals are frequently zero or inconsistent; orientation is defined by vertex winding, which is what
  the half-edge structure uses.

## D20. Vertex welding: exact by default, optional epsilon (grid hash)
- **Choice:** `weld_vertices(soup, eps)`. `eps == 0`: exact coordinate match (with -0.0 normalized to +0.0).
  `eps > 0`: uniform grid with cell size eps, checking the 27 neighboring cells; first-representative rule.
- **Alternatives:** exact only; epsilon only; union-find clustering (transitive, can chain-merge distant points).
- **Reason:** exact can never merge distinct vertices and fits well-formed STL; epsilon repairs noisy exports.
  The first-representative rule is deterministic and bounds drift: every merged point is within eps of its
  representative, which union-find does not guarantee.

## D21. Welding never drops triangles
- **Reason:** silently deleting geometry hides defects. Collapsed (degenerate) triangles are kept and reported by
  the topology diagnostics (M3).

## Proposals (not yet decided)
- **P1.** Adopted, see D22.
- **P2. Intrinsic Delaunay cotan Laplacian** (Sharp & Crane 2020, "A Laplacian for Nonmanifold Triangle Meshes",
  CGF / SGP 2020) as an M6 extension addressing the obtuse-triangle caveat.
- **P3. Signed heat method** (Feng & Crane 2024, "A Heat Method for Generalized Signed Distance", ACM TOG 43(4))
  as a far stretch after M6.
- **P4.** Adopted, see D23.
- **P5. Handle loops** via tree-cotree decomposition (Eppstein 2003, not yet read) for visualizing H1 generators.
- **P6. Explicit repair** (hole filling, handle cutting/patching): rejected for now as mesh modification adjacent to
  the out-of-scope remeshing; implicit repair via P3 preferred.

## D22. Genus and Betti numbers per component in M3
- **Choice:** M3 reports, per connected component, chi, boundary loop count b, genus g = (2 - chi - b)/2 and
  beta1, flagging non-integer or negative g as "not a manifold surface".
- **Reason:** nearly free given chi, components and boundary loops; detects spurious handles, a real scan defect.

## D23. Test data: synthetic plus CC0 only
- **Choice:** generated meshes, plus selected CC0 models from Keenan Crane's model repository, each with its
  license recorded in the repo. No dental scan data.
- **Reason:** unambiguous licensing.

## D24. Shared private number parser
- `detail/parse_number.h` (not under `include/`, so not part of the public API) holds `parse_double` for both text
  parsers, so OBJ and ASCII STL cannot drift apart in how they read numbers.

## D25. Half-edge layout: index-based, implicit next/face; bowtie vertices rejected
- **Choice:** half-edge `3f+k` runs `tri[f][k] -> tri[f][(k+1)%3]`; `next`, `prev`, `face` are arithmetic. Stored per
  half-edge: `origin`, `twin`; per vertex: one outgoing half-edge (the boundary one for boundary vertices).
  Structure-of-arrays of `uint32_t`, `kInvalid = UINT32_MAX` for "none".
- **Alternatives:** pointer-based half-edge objects; storing `next`/`face` explicitly (needed for general polygon
  meshes); OpenMesh / CGAL / geometry-central.
- **Reason:** triangle-only (all loaders triangulate), so implicit `next` removes a stored array and a class of
  consistency bugs. Indices instead of pointers survive vector reallocation, are half the size of 64-bit pointers,
  and can be passed to JS / WASM memory as plain arrays (M5). Writing it ourselves is the point of the project.
- **Bowtie vertices** (faces around a vertex forming several fans) are rejected along with non-manifold edges (D10),
  so every vertex's star is a single fan and one-ring traversal is complete.

## D26. Mutation check for geometry tests
- **Practice:** after a routine passes, plant the known classic bugs and confirm the tests fail.
- **M2 result:** `next(h) = h + 1` -> 6 of 30 tests fail; clockwise rotation in `one_ring` -> 3 of 30 fail.
- **Reason:** a passing suite only means something if it can fail.

## D27. Synthetic mesh generators live in the core
- `make_grid`, `make_grid_with_holes`, `make_torus`, `make_icosphere`, `make_mobius`, `append`. Each tested
  against closed-form V/E/F counts and, for closed meshes, outward winding via positive signed volume.
- **Reason:** the viewer (M5), curvature convergence tests (M4) and benchmarks (M7) reuse them.

## D28. Components by vertex connectivity
- **Alternatives:** edge (face-adjacency) connectivity.
- **Reason:** beta0 counts path components of the space; two triangles meeting at a single vertex are
  path-connected. Edge connectivity would report a bowtie as two components. Bowties are flagged separately as
  non-manifold vertices.

## D29. Genus and Betti numbers only for manifold, orientable components
- g = (2 - chi - b)/2, beta = (1, b > 0 ? 2g + b - 1 : 2g, b == 0 ? 1 : 0). Consistent winding is NOT required:
  chi and b are combinatorial and unchanged by flipping faces, so a mis-wound but orientable surface gets a genus
  and a separate "inconsistently oriented" flag.
- **Reason:** the formula comes from the classification of compact orientable surfaces; outside its hypotheses
  it produces meaningless numbers.

## D30. Orientability by breadth-first orientation propagation
- Choose flip[f] in {0,1} per face so every two-face edge becomes consistently wound; a contradiction means
  non-orientable (e.g. Moebius). O(F).
- **Reason:** distinguishes fixable winding errors from genuinely non-orientable surfaces, which need different
  remedies.

## D31. Invalid and duplicate faces are reported and excluded
- Out-of-range / repeated-index faces and later copies of a face's vertex set are listed, then ignored, so all
  counts describe the underlying simplicial complex.

## D32. Geometric degeneracy deferred
- Zero-area and sliver triangles need a tolerance; they are handled where they matter (curvature, M4).

## D33. Shared private helpers
- `kInvalid` moved to `mesh.h`; `detail/hash.h` (Hash3, pack_pair) and `detail/disjoint_sets.h` (union-find) are
  shared by welding and topology.
- **Mutation check (D26) for M3:** orientability ignoring the flip rule -> 1/50 fail (only the Moebius test:
  thin coverage, noted); counting every boundary vertex as a loop -> 3/50; inverted misorientation test -> 3/50.

## D34. Curvature per Meyer, Desbrun, Schroeder, Barr (2003)
- Mixed Voronoi area, cotan mean-curvature normal K(x_i) = (1/2A) sum (cot a + cot b)(x_i - x_j) = 2Hn, Gaussian
  curvature as angle defect / A. (Author order corrected from the original brief: Meyer, Desbrun, Schroeder, Barr.)
- **Implementation:** per-face scatter instead of per-vertex one-ring gather; each edge's two opposite-angle
  cotangents arrive from its two faces, so the sums are identical. Angles via atan2(|u x w|, u.w), cotangents as
  u.w / |u x w|: no acos (ill-conditioned near 0 and pi) and no trig for cot.

## D35. Pointwise curvature = integrated quantity / mixed area
- Discrete curvature is naturally an integral over a vertex's region; dividing by A_mixed gives a pointwise
  estimate, as in Meyer et al.

## D36. Boundary vertices: angle defect pi - sum(theta); no pointwise curvature
- The boundary term makes sum(defects) = 2 pi chi hold with boundary. Mean curvature normal is undefined at the
  boundary (one-sided ring), so H, K, k1, k2 are NaN there; viewers must treat NaN as "no data".

## D37. Sign of H from area-weighted face normals
- **Alternatives:** use K(x_i)'s direction (undefined where H = 0); angle-weighted normals.
- **Reason:** sign needs a normal independent of K(x_i) that follows the winding. Measured consequence: on the
  icosphere, all H error for s >= 2 comes from this normal deviating from radial (see D39).

## D38. Zero-area faces
- Listed in `degenerate_faces`; their angles (0, 0, pi) still count toward angle sums so Gauss-Bonnet holds, but
  their (infinite) cotangents are skipped. Near-degenerate slivers are NOT special-cased; they produce large
  cotangents, which is the motivation for M6b.

## D39. Measured accuracy (2026-10-05, Apple clang 17, -O2 probe; max error over interior vertices)
- Icosphere r = 2: H error 2.2e-16 at s=1 (exact by symmetry), then 5.7e-5, 1.4e-5, 3.6e-6, 9.0e-7, 2.3e-7 for
  s = 2..6 (order 2). K error 2.0e-2, 5.1e-3, 1.4e-3, 3.5e-4, 8.9e-5, 2.2e-5 for s = 1..6 (order 2).
- Cylinder r = 0.5: H = 1/(2r) and K = 0 to roundoff (1e-16 .. 2e-13) at every resolution: an identity, not
  convergence.
- Torus R = 1, r = 0.3, nu = 3 nv: H error 3.6e-2, 9.8e-3, 2.5e-3, 6.3e-4, 1.6e-4 and K error 1.9e-1, 5.0e-2,
  1.3e-2, 3.2e-3, 8.0e-4 for nv = 8..128 (order 2).
- **Limitation:** icosphere with vertices jittered by ~0.25h along the sphere: max H error ~0.28 and mean ~9e-3 at
  s = 4, 5, 6, i.e. no convergence. With 0.1h jitter the mean converges but the max stalls near 1e-2. Pinned by a
  characterization test as the baseline for M6b.

## D40. Viewer generators: plate with a handle, defect showcase
- `make_plate_with_handle`: two holes joined by an arched square tube (chi = -1, b = 1, genus 1). The tube's last
  ring is attached in mirrored corner order; the natural order was verified to produce misoriented edges.
- `make_defect_showcase`: one mesh with every defect class at known indices, tested, used as the viewer demo.
- The M1 placeholder `add()` and its test were removed.

## D41. JS <-> WASM data transfer
- **Input:** JS `_malloc`s inside WASM, copies bytes with `HEAPU8.set`, passes the integer offset; C++ reads via
  `std::span` in place; JS frees. One bulk copy. HEAPU8 is read after malloc because malloc may grow memory and
  detach older views.
- **Output:** C++ returns `typed_memory_view`s (zero-copy JS views into WASM memory); JS copies each with
  `.slice()` immediately, since the views die on the next call that grows memory.
- **Alternatives:** embind `std::vector` / `val` conversions (per-element), returning JSON for bulk data.
- Stats (small) are returned as hand-built JSON.

## D42. Generated TypeScript declarations (D8 revisited)
- `emcc --emit-tsd` generates `dmw.d.ts` from the embind registry; copied as `dmw.d.mts` next to `dmw.mjs`. Removes
  the hand-synced declaration file. Typed-array accessors come out as `any` and are typed in `dmw.ts`.

## D43. Three.js 0.186.1 viewer choices
- Per-vertex fields on indexed geometry; per-face colors (components) on non-indexed geometry, since an indexed
  vertex is shared by faces of different colors.
- Fat lines (`LineSegments2`) because WebGL ignores line widths > 1; polygon offset on faces to avoid z-fighting;
  non-manifold and misoriented edges and non-manifold vertices drawn without depth test so they are never hidden.
- Diverging cool-warm color map, symmetric range at the 98th percentile of |value| (curvature is heavy-tailed);
  NaN (boundary) drawn gray.
- GPU geometries and materials are disposed explicitly on every rebuild.
- **Known gaps (for M7):** 617 kB bundle (mostly Three.js); the file picker is a hidden input inside a label and is
  not keyboard-accessible.

## D44. Cross-target smoke test
- `npm run smoke` runs the WASM build under Node and checks genus, Gauss-Bonnet totals, defect detection, and the
  STL/OBJ byte paths against the values the native tests establish.
- **Verified in the browser (2026-10-05):** plate with handle (chi -1, b 1, g 1, only the outer rim highlighted),
  defect showcase overlays, torus Gaussian curvature sign pattern; no console errors.

## D45. DEC operators (M6)
- d0 (signed vertex-edge incidence), d1 (signed edge-face incidence), *0 (barycentric lumped areas), *1
  ((cot a + cot b)/2), L = -d0^T *1 d0 assembled as a weighted Gram product. Tests: d1 d0 = 0 exactly, L 1 = 0,
  L symmetric NSD, linear precision on flat irregular meshes, div(grad u) = L u exactly, and L x = -A K(x)
  matching M4's independent implementation to 1e-12.

## D46. Hand-written CSR + Jacobi-preconditioned CG
- As the original brief specified ("hand-written conjugate gradient first").

## D47. FINDING: CG cannot solve the heat step accurately enough at t = h^2 (measured 2026-10-05)
- Heat-method error grew under refinement with CG at relative tolerance 1e-10 (icosphere mean error 2.3e-2 at
  s=3 but 6.5e-1 at s=5, 1.09 at s=6; grid likewise beyond n=64), while Dijkstra stayed ~1.2e-1.
- Cause: one backward-Euler heat step decays roughly like exp(-d/sqrt(t)) = exp(-d/h); far vertices have u many
  orders of magnitude below CG's residual tolerance, so u there is solver noise, and step II normalizes noise into
  arbitrary directions.
- Evidence (icosphere s=4, same matrix): CG 1e-10 -> 841/2562 vertices with u <= 0, distance error mean 1.4e-1,
  max 1.2; CG 1e-14 -> 106 such vertices, mean 1.8e-2, max 3.9e-1; dense Cholesky -> 0 such vertices, mean
  1.3e-2, max 2.9e-2. Larger t (m = 16) also hides it (mean 5.1e-3 at s=5) but gives up t ~ h^2 convergence.
- Conclusion: the heat step needs a direct (factorization) solver. Decision on how: D49.

## D49. Direct solver: hand-written envelope LDL^T with reverse Cuthill-McKee ordering (author's choice)
- **Alternatives:** Eigen SimplicialLDLT (AMD ordering; fast at 40k+ vertices, MPL2 dependency); keep CG with a
  larger t (gives up convergence).
- **Reason:** no dependency and fully explainable; factor once, two triangular solves per query (the paper's
  prefactorization). Envelope fill is O(n * bandwidth), roughly n^1.5 on surface meshes: fine for viewer-size
  meshes, heavy beyond ~40k vertices. Eigen remains the upgrade path.
- The Poisson matrix -L is made SPD by pinning one vertex per connected component (phi = 0 there; exact, since
  the system is consistent and phi is shifted afterwards). Vertices without faces get a decoupled unit diagonal.

## D48. Heat method: Neumann boundary conditions, t = m h^2 with m = 1, h = mean edge length
- Neumann (natural) conditions come for free with the cotan Laplacian. The paper also discusses Dirichlet and
  an average of the two near boundaries; not implemented (open).

## D50. Measured accuracy of the heat method (2026-10-05, Apple clang 17 -O2 probe, source vertex 5)
- Icosphere r = 1, mean / max |error| vs great-circle distance, heat vs Dijkstra-on-edges:
  s=2: 4.5e-2 / 9.7e-2 vs 1.1e-1; s=3: 2.3e-2 / 5.1e-2 vs 1.2e-1; s=4: 1.3e-2 / 2.9e-2 vs 1.2e-1;
  s=5: 8.7e-3 / 1.9e-2 vs 1.2e-1; s=6: 6.5e-3 / 1.4e-2 vs 1.2e-1. Dijkstra does not converge.
- Observed heat convergence rate slows with refinement (~0.96, 0.81, 0.61, 0.42 per level). Cause not yet
  identified (candidates: source singularity, far-field heat precision). Open.
- Flat unit grid, source at center, vs Euclidean: heat mean 2.2e-2 (n=8) -> 4.4e-3 (n=128); Dijkstra ~7e-2.
- Jittered (0.25h) icosphere: heat mean 2.9e-2 (s=3) -> 1.0e-2 (s=6): irregularity no longer breaks it.
- Larger t at s=5: m=0.25 -> 1.7e-2, m=1 -> 8.7e-3, m=4 -> 5.8e-3, m=16 -> 5.1e-3 (mean); on this smooth closed
  surface, more smoothing helped.
- Cost (native): envelope factor entries / setup / query: s=4 0.35M / 78 ms / 5 ms; s=5 2.7M / 0.49 s / 28 ms;
  s=6 21.9M / 4.9 s / 172 ms. Browser (WASM), plate with handle, 381 vertices: first query (factor + solve)
  4.7 ms, later queries 0.3 ms.
- Mutation check (D26): pinned rhs not zeroed -> 1/95 fail (only the non-pinned-source test, added for it);
  X not negated -> 4/95; *1 wrong corner -> 7/95; d1 sign ignored -> 1/95; LDL^T overlap sum dropped -> 6/95.
  The harness now aborts if a mutant fails to build (an earlier run silently re-tested the old binary).

## D51. Geodesic distance in the viewer
- Session builds HeatGeodesics lazily on the first query and keeps the half-edge mesh on the heap (the solver
  holds a pointer to it). Click picks the nearest vertex of the hit face. Isolines via a 1D striped texture
  indexed by distance / max (interpolated per pixel, so lines stay crisp).
- Fixed: a uint32 -> iterator offset in SparseMatrix::at() broke only the wasm32 build (-Wsign-conversion,
  32-bit difference_type).
- Viewer message for meshes without half-edge structure now says "consistently oriented manifold" and that gray
  means no data (the Moebius strip was being read as "curvature ~ 0").

## D52. Length-only geometry: IntrinsicTriangulation and a single DEC builder
- Connectivity + one length per half-edge; areas (Kahan's stable Heron), angles and cotangents (law of cosines),
  gradient and divergence (per-face 2D layout) from lengths only. `build_dec(HalfEdgeMesh)` now delegates to
  `build_dec(intrinsic_from_mesh(mesh))`: one implementation; all M6 DEC tests pass through it.
- The heat method now owns an IntrinsicTriangulation (no pointer to the caller's mesh).

## D53. Intrinsic Delaunay by edge flipping (Fisher, Springborn, Bobenko, Schroeder 2007)
- Flip interior edges with cot a + cot b < -1e-12; new diagonal length by unfolding the two triangles. The operator
  is the intrinsic Delaunay Laplacian of Bobenko & Springborn (Discrete Comput. Geom. 2007). Sharp & Crane 2020
  ("A Laplacian for Nonmanifold Triangle Meshes") generalize it to non-manifold input; not implemented.
- **Limitation:** an intrinsic Delaunay triangulation may need self-loops or multi-edges (Fisher et al. note
  faces with only two distinct edges). Our vertex-triple connectivity cannot represent them, so such flips are
  skipped and counted (`skipped`). Measured: 0 skipped on all test meshes except 1 on a 0.45h-jittered s=5
  icosphere.
- Tests: flipped connectivity equals a from-scratch half-edge rebuild; total area and every cone angle preserved;
  edge count preserved; all weights >= 0 when nothing is skipped; regular icosphere needs no flips; linear
  precision on flat meshes; Euclidean lengths reproduce the extrinsic operators to 1e-12.

## D54. Intrinsic Delaunay option in curvature and heat geodesics
- `compute_curvature(mesh, true)`: K(x_i) = -(L_idt x)_i / A_i with intrinsic mixed area; angle defects unchanged.
- `HeatGeodesics(mesh, t, true)`: all operators from the intrinsic Delaunay triangulation.

## D55. Measured effect of intrinsic Delaunay (2026-10-05, native -O2 probes)
- **Maximum principle (flat "brick" grid, alternate rows shifted 0.45 cell, min cotan weight -0.90):** heat from
  a point source has 54 (n=20) and 446 (n=40) negative values with cotan, none with intrinsic Delaunay.
  Heat-geodesic mean error: cotan 4.3e-2 -> 1.07e-1 (gets worse under refinement); intrinsic Delaunay
  9.9e-3 -> 5.9e-3 (converges). At shift 0.30: no negative heat, but cotan stalls at ~2e-2 while intrinsic
  Delaunay converges 8.6e-3 -> 4.7e-3.
- **Jittered icosphere 0.45h (s=3,4,5):** geodesic mean error 4.3e-2/2.2e-2/1.9e-2 -> 2.5e-2/1.8e-2/1.6e-2.
  Max |H - 1|: 3.0/8.9/16 -> 0.79/1.1/2.2 (about 8x lower, but still growing); mean |H - 1| 2.5-3x lower.
- **Jittered icosphere 0.25h:** geodesics essentially unchanged (2.9e-2 -> 2.3e-2 at s=3, equal at s=4,5);
  max |H - 1| about 0.2 either way.
- **Conclusion:** intrinsic Delaunay fixes negative weights and their consequences (lost maximum principle,
  non-convergence of heat geodesics, curvature blow-ups). It does not make pointwise mean curvature converge on
  irregular meshes (the M4 limitation is reduced, not removed).

## D56. Intrinsic Delaunay in the viewer; mutation-harness hardening
- Session `setIntrinsicDelaunay(bool)` recomputes curvature and resets the geodesic solver; stats report the flip
  count. New preset `brick` (`make_brick_grid(30, 0.45)`). **Verified in the browser (2026-10-05):** geodesic
  distance from a corner of the brick grid: cotan gives swirling, broken isolines and max distance 0.505 (true
  value sqrt(2) = 1.414); intrinsic Delaunay (450 flips) gives concentric quarter-circles and max 1.436.
- Mutation results for the flip code (D26): unfolding d on c's side -> 6/106 fail; inverted Delaunay test ->
  infinite flipping, caught by per-test timeouts; vertex remap skipped -> 2/106; twin redirect skipped -> 7/106;
  duplicate-edge guard removed -> initially SURVIVED (no test exercised it), now 1/106 after adding the
  flattened-tetrahedron test.
- Harness lessons: per-test timeouts (mutants can hang); a killed run leaves the mutant on disk (restore from the
  saved copy); make's 1-second timestamp comparison can leave a stale object, so the harness must verify the
  mutated file was actually recompiled.

## D57. Benchmarks
- `tools/dmw_bench` (native, Release build directory, median of repeated runs) times each stage on icospheres
  s=3..7 and writes them as binary STL; `web/src/bench.ts` times the WASM build on the same files under Node.
  `write_obj` (17 significant digits, round-trips doubles exactly) and `write_stl_binary` added, round-trip
  tested. The full test suite also passes in the Release build (UB that only appears under -O3 would show).
- Results (2026-10-05, Apple M4 16 GB, Apple clang 17 -O3; WASM emsdk 6.0.11, Node 26.10) are in README.md.

## D58. FINDING (corrected): Emscripten link flags
- I first concluded the WASM module was linked without optimization (CMAKE_EXE_LINKER_FLAGS_RELEASE is empty) and
  added -O3 at link. The output was byte-identical: link.txt shows CMake already passes CMAKE_CXX_FLAGS_RELEASE
  (-O3) on the em++ link line. Change reverted. Also measured: WASM memory growth is not the cause of the 4.3x
  WASM/native gap at 164k vertices (later runs grow the heap less yet take the same time). Cause open.

## D59. Publishing: MIT license, GitHub Pages workflow (author's choices)
- `.github/workflows/ci.yml`: native Release build + ctest (clang on ubuntu), pinned emsdk 6.0.11 + WASM build +
  typecheck + Node smoke test + Vite build, then Pages deploy on pushes to main. Action versions: checkout@v6,
  setup-node@v6, configure-pages@v5, upload-pages-artifact@v4, deploy-pages@v4 (looked up 2026-10-05).
  Not run yet: nothing has been pushed.
- Vite `base: "./"` so the site works under /<repo>/; verified by serving the production build locally (WASM
  loads through the relative hashed asset path).
- File picker is now a real <button> forwarding to the hidden input (keyboard-accessible; M5 gap closed).
  Bundle size left as a documented limitation.

## D60. Data policy change: public dental scans allowed (author's decision, 2026-10-05)
- Supersedes D23: public, legally licensed dental scan data is allowed (the author's decision).
- **License findings (checked 2026-10-05):**
  - **Teeth3DS+** (1,800 intraoral scans, 900 patients; OBJ + JSON per-vertex FDI labels, 0 = gingiva; OSF
    osf.io/xctdy): **CC BY-NC-ND 4.0** per the official challenge repository (github.com/abenhamadou/
    3dteethseg22_challenge); its code is MIT. A ResearchGate listing saying CC BY 4.0 conflicts; the official
    source wins.
  - **3DTeethLand landmarks** (340 scans; Mesial, Distal, Cusp, Inner, Outer, Facial points as 3D coordinates
    in JSON; OSF osf.io/gvjah): **CC BY-NC-ND**. (The Zenodo record 10991302 is CC BY 4.0 but holds only the
    challenge description PDF.)
  - **Tsinghua pre/post-orthodontic dataset** (Zenodo 11392406, labelled CC0): access requires a signed data use
    agreement restricting use to academic purposes. **Not used** (the author works in industry).
  - **FDI 16** (single-tooth meshes, reported CC BY-NC-SA): details unverified; not used yet.
- **Consequences of NC-ND:** use locally for non-commercial evaluation and publish aggregate metrics with
  attribution; do NOT commit or host scans or any modified/derived meshes (decimated, converted, processed).
  Data lives outside the repo (`data/` is git-ignored); the public demo hosts no scans (users open their own
  files locally, as the viewer already allows).

## D61. Scope change: dental feature-line extraction; ML not blocked
- Tooth-gingiva margin and cusp/ridge detection by geometric feature-line extraction first. ML segmentation is no
  longer out of scope: if geometry plateaus and ML is the obvious next step, it proceeds (open question to
  resolve then: whether trained weights count as derivatives under NC-ND). Optimization for real scan sizes is in
  scope.

## D62. Real-scan QA (3DTeethLand test set, 100 scans = Teeth3DS+ data part 7; local only)
- `tools/scan_qa` runs the viewer pipeline over every OBJ in a directory; per-scan CSV stays local (D60).
- The downloaded part 7 is exactly the 3DTeethLand **test** set (all 100 IDs in its testing lists, none in the
  Teeth3DS splits) and contains **no labels or landmarks**: margin and cusp ground truth still need the
  labelled Teeth3DS parts and the 3DTeethLand landmark files.
- Measured (2026-10-05, native Release, Apple M4): median 106k vertices (75k-216k); median parse 23 ms,
  topology 65 ms, half-edge 67 ms, curvature 16 ms per scan. Before D63, **25/100** scans failed the half-edge
  build (24 because of collapsed faces, i.e. all three corners the same vertex, typically exactly 64 per scan);
  after D63, **3/100** (genuinely non-manifold). 97/100 scans have more than one component (92 have 5+: debris).
  Largest component genus: 0 in only 9 scans, up to 26; boundary loops 1-5. Exact-welding three scans (genus 0,
  4, 26) changed neither vertex count nor genus, so the handles are not unwelded-vertex artifacts.
- Debugging note: a probe crash turned out to be zsh not word-splitting an unquoted `$F` (three paths passed as
  one) plus the probe dereferencing a null "largest component" for an empty mesh. Sanitizer builds of the core
  on real scans were clean.

## D63. Exclude invalid and duplicate faces before building the half-edge structure
- `without_excluded_faces(mesh, topology_report)` drops the faces topology already reports (D31), keeping vertex
  indexing. Used by the WASM session and scan_qa. Non-manifold edges and vertices are still rejected (D10).
