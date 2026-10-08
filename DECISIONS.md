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

## D64. Handle loops (M8b): tree-cotree with the greedy shortest system of loops
- Eppstein (SODA 2003) tree-cotree decomposition; Erickson & Whittlesey (SODA 2005) greedy choice: T = Dijkstra
  shortest-path tree per component, C = maximum spanning dual tree on loop length sigma(e) = d(u) + l(e) + d(v);
  leftover edges close 2g loops, each trimmed to u -> lca -> v -> u. Boundary loops are capped by one virtual
  dual node each, so holes are never reported.
- Tests: none on genus-0 meshes even with holes (b = 2, 3); torus: 2 independent non-separating loops, shortest =
  the tube 8-gon exactly (1e-9); per-component counts; every loop a valid closed edge cycle.
- **Limitation (measured):** loops are valid generators but NOT shortest in their homology class. On the plate
  with a handle, the loop homologous to the tube's cross-section circles hole A at length 1.158 vs 0.333 for the
  tight cycle. Reported lengths are therefore upper bounds on handle size. Tightening (shortest homologous
  cycles) is future work.
- **Real scans (100, native Release):** 1,062 loops; count = 2 * sum(genus) on every scan (exact invariant on real
  data); median 24 ms per scan (p95 34). Loop length p10 0.47 mm, median 1.91 mm, p90 21.8 mm, max 126 mm. In the
  browser (WASM), scan 3TROMS4N_lower (93,610 vertices, 17.7 MB OBJ) loads and is fully analysed, including handle
  loops, in 382 ms; the loops visibly cluster at interproximal contacts (qualitative, consistent with fused
  neighbouring crowns).
- Viewer: "Handle loops" layer and stats; dev-server-only `server.fs.allow: ["../data"]` to open local scans.

## D65. Minimal JSON parser in the core
- Recursive descent, depth-limited, errors with byte offsets; text in, value out (D14/D15). Used for 3DTeethLand
  landmarks and Teeth3DS labels. Tested on the landmark shape, escapes, literals, malformed input.

## D66. Point-detection metrics
- Greedy one-to-one matching by ascending distance within a tolerance; precision, recall, F1, matched distances;
  micro-averaged across scans. (Not the 3DTeethLand challenge's official metric: not comparable to its leaderboard.)

## D67. Cusp detectors
- A raw Meyer mean curvature (baseline); B mean curvature diffused to scale sigma; C occlusal prominence
  p = h - diffuse(h, R) with h the height along the occlusal axis (smallest PCA axis of the arch, signed away from
  the gingival cut), gated by smoothed H > 0 and h above an arch-height quantile. One-ring local maxima, then NMS
  (spatial hash). Only the largest component (the arch) is searched. Synthetic tests: Gaussian "cusps" with 10 um
  roughness: C finds all tips with precision 1 at 0.3 mm, B all tips, A precision < 0.5.

## D68. Cusp evaluation protocol and results (2026-10-05)
- Tuning: 67 training scans with 3DTeethLand landmarks (66 in Teeth3DS part 1, 1 in the sample); grids widened
  over 4 sweeps until each optimum was interior (except B's sigma = 3.5 mm at the range edge, +0.007 F1 over 2.5:
  diminishing returns, stopped). Chosen by F1 at 1 mm on training: A thr 2.0/nms 2.0 (0.101); B sigma 3.5/thr
  0.15/nms 4.0 (0.596); C sigma 0.5/R 4.0/q 0.5/thr 1.3/nms 3.5 (0.663).
- **Held-out test (100 scans = 3DTeethLand test set, 2,343 Cusp landmarks), run once with the fixed points:**
  A: P 0.062 R 0.534 F1@1mm 0.112 (200 detections/scan); B: P 0.521 R 0.656 F1@1mm 0.581, median error 0.51 mm;
  **C: P 0.553 R 0.736 F1@0.5/1/2mm 0.368/0.631/0.744, median error 0.44 mm**, 31.2 detections/scan (truth: 23.4).
  Train-to-test gap for C: 0.663 -> 0.631. The library path (detect_cusps) reproduces the sweep path exactly.
- Visible failure mode: false positives on incisal edges of anterior teeth, which carry no Cusp landmarks.
  Tooth-type awareness (margin/segmentation, or ML) is the obvious next lever for precision.

## D69. Manifold analysis view (supersedes the reject-whole-scan behaviour of D10 for real data)
- `manifold_analysis_mesh`: iteratively excludes invalid/duplicate faces and faces touching non-manifold or
  misoriented edges or non-manifold vertices until the half-edge build succeeds; returns the built structure.
  Measured: training scans with half-edge failures 22/67 -> 0/67; test scans 3/100 -> 0/100. The topology report
  still describes the full input; the viewer reports the excluded-face count.

## D70. Cusp detector speed: solve for prominence directly, at a measured tolerance
- (*0 - tL) p = -t L h instead of reconstructing the diffused height field: at loose tolerance this is far more
  accurate (full-field at 1e-4 was off by 1.5 mm vs 5e-4 mm for the direct solve), because CG's tolerance is
  relative to the right-hand side, which is now small (same lesson as D47, from the other side).
- Tolerance 1e-4 for the detector's solves: max prominence error 5.2e-4 mm against a 1e-11 reference (threshold
  1.3 mm; scanner accuracy 10-90 um per the Teeth3DS+ paper).
- Measured on a 93.6k-vertex scan: native detect_cusps 1,525 -> 383 ms (4.0x); in the browser (WASM)
  3,685 -> 459 ms (8x). Held-out test metrics identical to 3 decimals after the change.

## D71. Margin detection: concavity-weighted geodesic Voronoi
- Binary tooth/gingiva labelling; the margin is its boundary (same definition as the ground truth). Seeds: teeth at
  detected cusp tips (D68) and the highest 10% of the arch; gingiva at the scan's cut boundary and the lowest 15%.
  Two-label multi-source Dijkstra with edge cost l * (1 + alpha * s), s = max(0, -kappa_min). Baseline: a plane
  cut at a height quantile along the occlusal axis. Synthetic test: crown with a concave crease at its base; the
  weighted Voronoi lands within 0.15 mm ASSD and at least 2x closer than a plane cut.
- Building blocks (margin_inputs, valley_strength, margin_labels) are shared by the library and the evaluation tool;
  operating points come from single functions (cusp_operating_point, margin_operating_point).

## D72. Margin metrics
- Margin edges sampled at midpoints; nearest-neighbour distances via a uniform grid with ring search. ASSD, HD95,
  Hausdorff, boundary F1 at 0.25 and 0.5 mm, area-weighted tooth IoU. Per-scan, averaged over scans.

## D73. Margin evaluation (2026-10-06)
- Data: Teeth3DS labels (FDI per vertex, 0 = gingiva). Range-probing the zip central directories (803 KB fetched)
  showed data parts 5 and 6 hold the whole Teeth3DS test split (300 each), parts 1-4 only training; part 5 was
  downloaded (author approved).
- Tuning: 60 training scans (every 5th of part 1), three sweeps until the optimum was interior or at a natural end:
  best alpha 2560 (converged: 1280 -> 2560 changes ASSD by 0.004 mm), sigma 0 (raw kappa_min; the end of the
  range: the cervical crease is narrow and smoothing blurs it), gingiva quantile 0.15, tooth quantile 0.9, cusp
  seeds on. Training ASSD 0.647 mm. Without cusp seeds: 2.54 / 1.91 mm (alpha 0 / 10).
- **Held-out test (300 scans of the Teeth3DS test split, run once):** plane cut ASSD 1.663 mm, F1@0.5 0.208, IoU
  0.670; Voronoi without concavity weighting ASSD 1.781 mm, F1@0.5 0.214, IoU 0.682; **concavity-weighted Voronoi
  ASSD 0.551 mm (median 0.450), HD95 3.27 mm, F1@0.25 0.761, F1@0.5 0.815, tooth IoU 0.841.** Library path
  reproduces all three exactly. Example test scan EJWZZZRF_lower in the browser: ASSD 0.257 mm, F1@0.5 0.92, 733 ms.
- Not yet analysed: which scans produce the HD95 tail (3.27 mm mean).

## D74. Margin failure analysis (validation set: 240 training scans never used in the sweeps)
- Tool `margin_analyze`: per-scan error decomposition (false vs missed margin, over- vs under-segmentation), error
  by tooth category, and scan characteristics. The held-out test set is not used for analysis or tuning.
- Findings: typical margin points are accurate (median missed-margin distance ~0.09 mm for every tooth category);
  error is a tail of localized gross failures (p95 2.6-3.6 mm; incisors mean 0.72, molars 0.61, canines 0.54,
  premolars 0.40 mm). Over-segmentation dominates: 14.9% vs 3.1% of true tooth area (gingiva labelled tooth vs
  the reverse); over-segmentation correlates most with per-scan ASSD (r = +0.53).
- Hypotheses tested: H1 cusp seeds on true gingiva: r = +0.38 with over-segmentation, +0.46 with ASSD; 179/240
  scans have >= 1 such seed (mean ASSD 0.648 vs 0.422 mm without). H2 resolution: rejected (mean edge length
  r = +0.09 / -0.07; vertex count correlates via scan extent). H3 handles: genus r = +0.48 with ASSD, +0.21 with
  over-segmentation (contributor, partly confounded). Second, separate mode seen visually: false margin loops at
  the distal ends behind the last molars on fragmented scans (HD95 up to 16 mm).

## D75. Seed-filtering experiment (validation) and oracle bound
- Label-free fix tried: keep cusp seeds only above an arch-height quantile. ASSD 0.590 (none), 0.585 (q 0.6),
  0.597 (0.7), 0.667 (0.8), 1.331 (0.9): height does not separate gingival false cusps from real ones. **Not
  adopted** (0.005 mm is noise; adopting it would overfit validation). `cusp_seed_quantile` defaults to 0.
- **Oracle** (seeds on true gingiva removed using labels, diagnostic only): ASSD 0.504 mm (-15%), HD95 3.13, IoU
  0.896. That is the ceiling for seed filtering, and it needs a better tooth/gingiva discriminator than height
  (local shape features or a small learned classifier, D61). Even the oracle leaves the other failure modes
  (distal ends, handles).

## D76. Cusp-seed classifier (logistic regression) and the second margin test evaluation
- Model: logistic regression trained by Newton/IRLS with L2 (core/learn), on six label-free features per seed:
  prominence, height quantile, vertex normal . occlusal axis, distance to the scan cut, smoothed mean and Gaussian
  curvature. Split of part 1 by index % 5: {0,1,2} training (180 scans, 5,409 seeds, 7.9% on gingiva), {3,4}
  validation (120 scans, 3,572 seeds, 6.8% on gingiva).
- Seed classification ROC AUC: training 0.904, validation 0.871. Strongest feature: normal . occlusal axis
  (single-feature AUC 0.849, standardized weight +1.49): real cusps face occlusally, gingival bumps sideways.
  Height quantile weight +0.07 (consistent with D75's failed height gate).
- Margin on validation (paired vs the operating point): threshold 0.5: -0.024 mm (t = -4.1, worse on 2/120);
  0.8: -0.030 (t = -3.4, worse 23); 0.85: -0.033 (t = -2.8, worse 34); 0.9: -0.004; 0.95: +0.161. **Chosen 0.5
  for reliability** (fewest regressions, strongest t; 0.85's mean advantage of 0.009 mm is not distinguishable);
  criterion stated explicitly because it was not pre-registered. Oracle on the same scans: -0.101 mm.
- **Second test evaluation (disclosed; the method was developed on validation only):** 300 scans: ASSD 0.551 ->
  **0.528 mm** (median 0.450 -> 0.428), HD95 3.27 -> 3.18, F1@0.5 0.815 -> 0.819, IoU 0.841 -> **0.869**; paired
  -0.024 mm, t = -5.3, better on 130, worse on 5 of 300. Baseline rows reproduce the first test run exactly.
- Licence handling: the model file (12 numbers + names) stays local in git-ignored data/models; the public build
  and viewer keep the classifier-free operating point until the author decides whether trained weights from
  CC BY-NC-ND data may be published.

## D77. Far-from-teeth false margins: height seeds leak, dropped (tooth_quantile 1.0)
- Decomposition (validation, 120 scans, classifier on, ASSD 0.560): false margin = predicted margin samples > 1 mm
  from the true margin, split by distance to the nearest true tooth vertex: "far" (> 2 mm, a spurious tooth region
  inside gingiva) is 8.3% of predicted samples on average, "near" (an offset line) 5.5%. Per scan, ASSD
  correlates +0.88 with the far share (HD95 +0.69): the tail of D74 is mostly false margin on gingiva.
- Sliver hypothesis rejected (LAFJJKAE_upper): correct margin edges have median min triangle angle 21.8 deg,
  false far edges 20.8 deg; median kappa_min -4.53 vs -0.93. The false lines sit on flat gingiva, not on bad
  triangles. Mechanism: first-arrival Voronoi; a tooth front that slips through a gap in the crease floods flat
  gingiva before the gingiva front arrives, leaving parallel stripes.
- Height seeds on true gingiva: median 0.1% but mean 2.3% of height seeds (a few scans, e.g. tilted or with high
  palatal gingiva, put large tooth seed sets on gingiva). Height-seed oracle (those removed with labels): 0.498.
- Experiment (validation, paired vs tq 0.9): tq 0.95 -0.031 mm (t -3.5, worse 33/120); 0.98 -0.047 (t -4.3, worse
  28); 1.0 -0.071 (t -3.9, worse 27). Both oracles (cusp + height seeds on gingiva removed): -0.158 (t -7.0).
- Confirmed on the 60 sweep scans (independent of the 120; same paired baseline 0.625): 0.98 -0.072 (t -3.4,
  worse 7/60); 1.0 -0.106 (t -3.4, worse 11/60), median 0.391, IoU 0.880; tq 1.0 without classifier 0.558 (worse
  23/60): the classifier still pays without height seeds.
- **Chosen tq 1.0** by the D76 criterion: regressions are close (validation 27 vs 28, sweep 11 vs 7; 38 vs 35 of
  180 combined), so the larger mean gain decides, and it is the simpler method (only the single highest arch
  vertex remains a height seed; cusp tips seed the teeth).
- Omission owned: the D73 sweep never tried tq above 0.9 (edge of the grid). The fix was available from the start.
- The principled anti-leak fix is a min-cut labelling (regional + boundary terms), not first arrival; not done.
- **Third test evaluation (disclosed, 2026-10-06; method fixed on validation + sweep scans before the run):** 300
  scans. Public point (no classifier) tq 0.9 -> 1.0: ASSD 0.551 -> **0.520 mm** (median 0.450 -> 0.409), HD95 3.27
  -> 2.95, IoU 0.841 -> 0.847; paired -0.032 mm, t = -3.6, better 178 / worse 122. With the classifier: 0.528 ->
  **0.491 mm** (median 0.392), HD95 3.18 -> **2.85**, F1@0.5 0.824, IoU **0.882**; paired -0.037, t = -4.0, better
  183 / worse 117. Classifier on top of tq 1.0: -0.029, t = -6.0, worse on only 11. Rows 0-1 reproduce D76 exactly.
- Honest reading: the gain is real (t -3.6 / -4.0) and hits the tail it was aimed at (HD95 -0.3 mm), but it is
  half the validation effect and the regression rate is higher (39-41% of test scans vs 18-23% on validation and
  sweep scans). The height-seed change trades many small losses for fewer large wins; it is not a dominance result
  like the classifier (5/300). Not re-tuned on test. Next: min-cut labelling, and per-scan regression sizes.

## D78. Graph-cut margin labelling (Boykov-Kolmogorov max-flow), and the fourth test evaluation
- Why: D77 traced the tail to first arrival: one gap in the crease lets a tooth front flood flat gingiva. A cut
  that pays for its boundary cannot do that cheaply. Energy over arch vertices (seeds as hard constraints):
  E = sum_v A_v D_v + mu * sum_cut l*_e / (1 + beta s_e); D from the two arrival distances, p = d_G / (d_T + d_G),
  D = -log p (tooth) or -log(1 - p) (gingiva); l*_e = cotan dual edge length (clamped at 0 on obtuse pairs);
  s = valley strength as in D71. After the geodesic graph cut of Price, Morse & Cohen (CVPR 2010), with seeds
  as hard constraints (Boykov & Jolly, ICCV 2001), moved from pixels to a mesh.
- Solver: Boykov & Kolmogorov (PAMI 2004), written from scratch (core/maxflow). Tests: CLRS Fig. 26.1 (23),
  brute-force min cut on 300 random graphs, Edmonds-Karp on 600 sparse random graphs, max-flow = cut value of
  the returned labels. Mutants: 5 tried; missing neighbour re-activation on freeing survived the small dense
  graphs and is killed by the sparse Edmonds-Karp test; the distance heuristic and one redundant re-activation
  are equivalent (speed only / already covered by the scan cursor; the redundant line was removed).
- Synthetic test: one false tooth seed on flat gingiva: Voronoi fake region 6.94 mm^2; cut 5.12, 1.70, 0.73,
  0.056 at mu 0.3, 1, 3, 10. mu = 0 reproduces Voronoi (the unary alone is the arrival order).
- Bug caught by the sweep: margin_eval computed valley strength only for Voronoi, so the first graph-cut sweep ran
  without the crease term (all beta identical, ASSD ~1.73). Fixed; the test-mode library path check would have
  caught it, validate mode has none.
- Tuning (60 sweep scans, paired vs Voronoi 0.558): sweep 1 mu {0.3..10} x beta {0, 100, 1000}: beta 0 worse at
  every mu (crease term essential); best at the mu edge. Sweep 2 mu {10..300}: 0.315 at (300, 100), edge again.
  Sweep 3 mu {300..1e6} x beta {30, 100, 300}: interior plateau along mu / beta ~ 3-10: (300,100) 0.315,
  (1000,100) 0.323, (1000,300) 0.323, (3000,300) 0.313, all within noise; mu 1e6 (unary negligible) 0.367, so the
  distance term helps. **Chosen (1000, 300), the plateau centre** (most stable neighbours); post-hoc criterion,
  stated as such.
- Validation (120 scans): Voronoi 0.518 -> cut **0.344 mm** (median 0.274), HD95 3.25 -> 2.38, IoU 0.872 -> 0.933;
  paired -0.174, t = -7.9, worse on 28 (2 by > 0.1 mm, largest +0.166). Classifier on top: 0.337.
- **Fourth test evaluation (disclosed, 2026-10-06; settings fixed before the run):** 300 scans. Voronoi 0.520 ->
  **GraphCut 0.338 mm** (median 0.409 -> 0.283), HD95 2.95 -> **2.34**, F1@0.25 0.764 -> 0.832, F1@0.5 0.819 ->
  **0.883**, IoU 0.847 -> **0.925**; paired -0.182 mm, t = -11.5, better on 247 / worse on 53, only 7 worse by
  > 0.1 mm (largest +0.191). Validation predicted it (-0.174 vs -0.182). With the classifier: 0.334 (adds -0.004,
  t -3.0): the cut makes the classifier nearly redundant, so the public build loses almost nothing without it.
  Rows 0 and 2 reproduce D77 exactly; library path agrees. Since the first test run: 0.551 -> 0.338 mm (-39%).
- Cost: labelling 36 -> 137 ms per scan native (two Dijkstras + max-flow on ~100k nodes), on top of ~1 s of
  shared inputs (curvature, cusps). Browser time not re-measured yet.
- Four test evaluations of the margin have now been run, each disclosed; the test set is no longer pristine
  in the strict sense (method choices were made after seeing earlier test numbers). Next margin claim should
  come from a fresh split or be flagged.

## D79. Where the margin result stands against published methods (2026-10-07)
- Published Teeth3DS results are supervised deep networks trained on ~1,200-1,440 labelled scans and scored per
  point over 17 classes (16 teeth + gingiva), mostly per-tooth IoU; none found reports tooth-gingiva boundary
  distance (ASSD/HD95 in mm). Closest numbers: CrossTooth (arXiv 2503.23702, 1,440/360 split) gingiva-class IoU
  96.41%, overall mIoU 95.86%, boundary IoU 82.06%; 3DTeethSAM (AAAI 2026, official 1,200/600) OA 95.48% (17-class),
  tooth-wise mIoU 91.90%; 3DTeethSeg'22 winner TSA 0.9859 (per-tooth F1). Prepared-tooth margin lines (a different
  task, private data, Alsheghri et al. 2024): median chamfer 0.137 mm, median Hausdorff 0.242 mm.
- To compare on a shared metric, margin_eval now also prints per-vertex (unweighted) binary metrics. Frozen
  methods, test split (300 scans), measurement only (no selection): GraphCut accuracy 0.9526, tooth IoU 0.9199,
  gingiva IoU **0.8967** (Voronoi 0.8656). Rows reproduce D78 exactly.
- Reading: on gingiva IoU the supervised state of the art is ~96.4% vs our 89.7%: about 3x our error (3.6% vs
  10.3%), on a different split, with ~1,400 labelled training scans vs our 3 tuned parameters and a 6-feature
  classifier that the cut no longer needs. The project's strength is the geometry, verification and evaluation
  discipline, not state-of-the-art accuracy. Boundary-distance metrics have no published counterpart to compare.

## D80. What limits gingiva IoU and F1 now: error decomposition, two cheap fixes tried and rejected
- Decomposition (validation, 120 scans, public graph cut; area-weighted): total error 1.6% of arch area. False
  gingiva in teeth missed entirely (< 50% covered) 18.9% of the error; in partially covered teeth 31.4%; false
  tooth 49.7%, of which only 15.8 points lie within 0.5 mm of a true tooth (boundary offset) and ~34 points
  farther (fake regions, mostly attached to real teeth). Missed teeth: 31 of 1,809 (incisors worst, 10/240);
  21 of the 31 carried a cusp seed, so the cut shrank around a point seed.
- Metric caveat: per-vertex gingiva IoU (0.909 here) is far below the area-weighted value (~0.98) because scans
  are much denser on teeth and near the margin: per-vertex scores magnify boundary error. Published per-point
  numbers use other samplings (D79), so cross-paper IoU comparisons carry this uncertainty.
- Tried, sweep scans then validation: (1) seed discs (each cusp seed grown to a geodesic disc, against the
  shrinking bias): worse at every radius, 0.5 mm +0.042 mm (t +2.4); discs also inflate false seeds. Removed.
  (2) Island removal (tooth components < A mm^2 -> gingiva, as in ToothGroupNetwork's post-processing): sweep
  A {2..20}, then {20, 40, 80}: best 20 mm^2, -0.027 mm (t -3.9, worse 3/60); 80 deletes real teeth. Validation:
  -0.005 mm (t -0.7), HD95 2.38 -> 2.43, one scan +0.69 mm. **Not adopted** (did not replicate; D75 rule).
  `min_tooth_region` stays, default 0 (off).
- Conclusion: the cheap geometric levers are spent. The remaining error (attached fake regions, missed and
  partially covered teeth) needs per-tooth evidence a learned model provides; the honest next comparisons are
  a state-of-the-art network scored by this tool on the same scans, and a fresh, never-used test set.

## D81. Fresh test data (Teeth3DS part 6), split audit, and a score mode for external predictions
- Split audit (official lists in the Teeth3DS archive): parts 5 and 6 are exactly the Teeth3DS test split (300 +
  300 jaws), part 1 is training. The split is per jaw, not per patient: 52 part-6 and 53 part-5 patients have
  their other jaw in part 1. Mild for this method (3 tuned parameters, a 6-feature classifier), stated anyway.
- Part 6 had never been used. **Frozen methods on part 6 (2026-10-07, first and only look):** 299 scans (1
  skipped). GraphCut **0.357 mm** (median 0.271), HD95 2.37, F1@0.25 0.829, F1@0.5 **0.880**, area IoU **0.927**;
  per-vertex accuracy 0.954, tooth IoU 0.920, gingiva IoU 0.901. Voronoi 0.580; paired -0.223 mm, t = -12.5,
  better on 258 / worse on 41 (3 by > 0.1 mm, largest +0.558). Classifier on top: 0.351 (-0.006, t -3.4).
  Part 5 gave 0.338 / 0.883 / 0.925: the result replicates on fresh data within 0.02 mm. This is the headline.
- ToothGroupNetwork comparison: its challenge checkpoints were trained on the 3DTeethSeg'22 training lists,
  which differ from the Teeth3DS split: 216 of part 6's scans are in the challenge TRAINING set. Only 84 part-6
  scans (90 in part 5) are in the challenge private test set; those 84 are unseen by both methods, so they are
  the comparison set. tools/tgn/prepare_input.py packages them; docs/TGN_COLAB.md runs the network on Colab.
- `margin_eval score <scan-dir> <pred-dir> [label]` scores any per-vertex prediction in the Teeth3DS label format
  (labels != 0 = tooth) with the same metrics, paired against the graph cut on the same scans; scans without a
  usable prediction (missing, or wrong vertex count) are skipped and counted. Checked: ground truth fed in as the
  prediction scores ASSD 0, F1 1, IoU 1; a truncated label file is skipped.

## D82. Why boundary F1@0.5 stays below 0.9: decomposition of the misses (validation, 120 scans)
- Boundary samples (edge midpoints) pooled over scans. Graph cut: precision 0.918, recall 0.855, pooled F1 0.885
  (the per-scan mean is 0.886). F1 is recall-limited.
- Recall misses, 14.5% of true samples: on missed teeth 1.5%; 0.5-1 mm off 3.6%; > 1 mm off on teeth that were
  found 9.4%. Of those 9.4%: the gingiva vertex beside the true edge is labelled tooth (over-extension) in 81.7%,
  it lies within 1.5 mm of two different teeth (interdental) in 61.0%, both in 60.1%. So **~5.6% of the true
  boundary is lost to interdental papillae labelled tooth (bridged between neighbouring crowns)**, the largest
  single F1 loss. It costs recall only: a bridged papilla removes true boundary without adding predicted boundary.
- Precision misses, 8.2% of predicted samples: fake regions on gingiva (tooth side > 0.5 mm from any tooth) 4.0%;
  0.5-1 mm off 1.1%; > 1 mm off but not fake 3.1%.
- Hypothesis tested and rejected: that the bridging comes from the cut's length term (going round a thin papilla
  costs more boundary than cutting across). Voronoi, which has no length term, bridges as much or more: recall
  0.801, far misses 15.1%, 44.9% of them bridged papillae = ~6.8% of the true boundary. The geometry around the
  papilla does not separate it for either method (plausibly weak or missing interproximal creases where the
  scanner cannot see between teeth; not verified).
- Fake regions are where the cut gained most over Voronoi: 10.8% -> 4.0% of predicted samples.
- Implication: the next F1 gain needs per-tooth knowledge (which papilla belongs between which two teeth), i.e.
  instance-level evidence from a learned model; tuning the cut will not recover papillae.

## D83. Interdental misses are not missing creases (validation, 120 scans)
- Measured kappa_min along every TRUE boundary edge (at the edge, and the strongest value within the one-rings of
  its endpoints), split interdental (gingiva side within 1.5 mm of two different teeth) vs cheek/tongue side,
  and matched vs missed by the graph cut; flat-gingiva baseline -0.24 /mm (21.7% of vertices below -1).
  Ring medians: cheek/tongue matched -4.49, missed -4.00; interdental matched -5.22, **missed -6.15** (q25 -8.52;
  98.9% below -1). The missed interdental boundary sits on the STRONGEST creases of the arch.
- Predicted boundary: only 2,621 predicted interdental edges lie > 1 mm from the true line, against 52,410 missed
  true interdental edges. The cut does not trace a competing crease elsewhere; it encloses the papilla in tooth and
  its boundary simply disappears (a narrow papilla keeps any replacement segment within 1 mm, so precision is
  unaffected). The competing-crease hypothesis (a shorter concave path between the crowns) is not supported.
- So the hypothesis "the scanner cannot see between the teeth, so there is no crease" (D82) is wrong too. The
  information is in the geometry. Working hypothesis for the mechanism (unverified): the tooth label enters the
  papilla where its two side creases converge at the tip, beneath the contact point, a crease gap like D77, and
  the arrival-distance unary then favours tooth over the whole wedge. Next: inspect the worst cases visually,
  then design (per-tooth labels, or a term that forbids entering a region enclosed by strong creases).
- Tooling: the interdental test first used per-tooth nearest-point queries, whose search expands across the
  arch for distant points; the run was stopped after 30 min. Replaced by a fixed-radius grid lookup (3 min).

## D84. Visual inspection of interdental misses: a sub-edge gingiva strip between touching crowns
- Tooling: the viewer gains a kappa_min overlay ("creases") and a dev-server-only URL loader
  (?scan=&labels=&focus=x,y,z&dist=&overlay=) that loads a local scan, detects the margin, overlays the labelled
  boundary and frames a point; margin_eval's errors experiment lists the 5 validation scans with the most missed
  interdental boundary, each with a viewer URL focused on its densest miss cluster.
- Seen on the two worst scans (QQKBAWLF_lower, JJEIMAPR_upper): between two adjacent crowns there is a narrow, deep
  valley (strong kappa_min across its whole width). The labels run a strip of gingiva up inside that valley,
  nearly to the contact point; the method labels the whole valley tooth and crosses at gum level. Handle loops
  (tunnels) pass beneath such contacts, consistent with crowns bridged by the scan at the contact.
- Measured (validation, 120 scans): at missed interdental samples, the distance to the neighbouring crown (another
  FDI label) is q25 0.11, **median 0.17**, q75 0.50 mm; 74.9% are below 0.5 mm. The labelled strip is about one
  vertex wide, below a typical edge length.
- Mechanism: in a BINARY tooth/gingiva labelling, two adjacent crowns are both "tooth", so the valley between them
  needs no boundary; keeping a one-vertex gingiva strip costs two long boundaries, merging costs none. Neither a
  crease term nor a better unary can express "tooth A | tooth B" in a two-label model. This supersedes the D83
  working hypothesis (entry at the papilla tip) and explains D82 (Voronoi bridges too).
- Design implication: per-tooth (multi-label) labelling makes the A|B boundary explicit; it runs along the valley,
  within 0.5 mm of both labelled strip edges for ~75% of these misses, which would recover most of the ~5.6% of
  true boundary lost here (upper bound for F1@0.5 roughly +0.03; not yet measured). Prerequisite: grouping cusp
  seeds into teeth.

## D85. Per-tooth labelling: cusp grouping + multi-label graph cut (alpha-expansion)
- From D84: two labels cannot express "tooth A | tooth B", so bridged interdental papillae cost no boundary.
  PerToothCut keeps the GraphCut result, groups the seed cusps into teeth and re-labels with gingiva + one label per
  group, minimizing the same energy (Potts boundary costs; unary p_l proportional to 1/d_l, which reduces to the
  binary p = d_G / (d_T + d_G) for two labels) by alpha-expansion (Boykov, Veksler, Zabih, PAMI 2001).
- Cusp grouping (group_cusps): tips connected inside the binary tooth region through vertices with kappa_min > tau.
  Checked against the FDI labels (60 sweep scans, 793 teeth): tau -1.25 pair precision 0.918, recall 0.884, clean
  teeth 86.1%; -1.0 82.0%; -1.5 85.6%; -2.0 74.0%.
- Alpha-expansion is its own module (core/multilabel), each move an exact s-t cut via the Kolmogorov-Zabih
  construction, using the D78 max-flow. Tests (brute force, n <= 7, up to 4 labels): energy never increases, the
  result is a local minimum over all expansion moves, within 2x the global optimum (the Potts bound); with two labels
  one move is the exact binary min cut. 4 mutants (sign of a linear term, dropped pair arc, swapped terminals, single
  sweep) all killed. Speed: nodes already labelled alpha are left out of each move (exact), and a candidate filter
  restricts a move to vertices within 3 mm grid cells of the label's current region (approximate): identical
  metrics on the sweep scans, 7.4 -> 2.5 s labelling per scan.
- Strip carving (the deeper vertex of each tooth|tooth edge becomes gingiva, reproducing the labels' one-vertex
  strip): **every variant worse** than the binary cut (+0.040..+0.080 mm, worse on 45-51 of 60). Off by default.
  The multi-label cut alone is what helps: separate tooth labels change both the unary and which boundaries are cheap.
- Sweeps (60 sweep scans, binary cut 0.323): tau -1.5 0.285; then -1.25 0.272 / -1.5 0.285 / -1.75 0.292 / -2.0
  0.311 (edge); then -0.75 0.276 / -1.0 0.276 / -1.25 0.272 / one label per tip 0.290 (worse on 21): **tau -1.25**,
  interior. Grouping matters.
- Validation (120 scans): 0.344 -> **0.283 mm** (median 0.237), HD95 2.38 -> 2.05, F1@0.5 0.886 -> **0.910**, IoU
  0.933 -> 0.940, per-vertex gingiva IoU 0.909 -> 0.915; paired -0.061, t = -4.7, better on 105, 3 worse by > 0.1 mm.
  F1 decomposition: recall 0.855 -> 0.888; recall misses > 1 mm off 9.4% -> 6.5%; fake regions 4.0% -> 2.3%;
  missed teeth 1.5% -> 1.9%.
- **Part 6 (second look; the first, D81, was the binary cut), 299 scans:** 0.357 -> **0.285 mm** (median 0.219),
  HD95 2.37 -> **1.98**, F1@0.25 0.829 -> 0.858, F1@0.5 0.880 -> **0.906**, IoU 0.927 -> **0.936**; per-vertex accuracy
  0.959, gingiva IoU 0.911; paired -0.072 mm, t = -7.3, better on 265, worse on 34 (1 by > 0.1 mm, largest +0.260).
  Validation predicted it (-0.061). Classifier on top: 0.282 (-0.004). Library path agrees.
- Cost: labelling 0.2 -> 2.2 s per scan native (300-scan mean); in the browser 6.8 s on a 164k-vertex scan, measured
  while another evaluation was running (an upper bound; to be re-measured). Now the public operating point.

## D86. Per-tooth cut speed: profile first, then two changes (one failed approach logged)
- Instrumented (MarginTimings, diagnostics only). 12 validation scans, labelling 1,715 ms/scan: alpha-expansion
  897, per-label Dijkstras 635 (17.7 labels per scan, each flooding the whole arch), binary stages 152, rest 26.
- Expansion: the candidate mask probed a hash map 27 times per node in every move (~50 moves per scan). Replaced
  by a dense cell grid dilated once per move: expansion 897 -> 450 ms (same run conditions), results identical.
- Per-label Dijkstras, attempt 1, FAILED: a cutoff on the crease-weighted distance (5-50 mm) destroyed the result
  (ASSD 0.26 -> 2.7-3.4 mm on all 12 scans): with alpha = 2560 one fissure adds thousands of weighted mm, so the
  cutoff truncated each tooth's field inside its own crown. The unary only uses ratios; absolute weighted distances
  have no scale. Attempt 2: each label's search only visits vertices within a straight-line radius of that label's
  seed tips. Same run: r = 0 1,131 ms, 20 714, **15 657**, 12 622, 9 587 ms; metrics identical (<= 0.001 mm per
  scan) down to 15, 12 within 0.005, 9 breaks (+0.059). **Default 15 mm.**
- Check on the 60 sweep scans: 0.274 mm, F1@0.5 0.916 (sweep 3 before the changes: 0.272 / 0.917); paired vs the
  binary cut -0.049, t = -3.1, better on 55 (unchanged). The radius changes results by ~0.002 mm on average.
- Net: labelling ~2.6x faster in like-for-like runs (~1.7 s -> ~0.65 s on the 12-scan slice); absolute timings in
  this session vary +-40% with machine load (a VM was using ~3 cores), so only same-run ratios are quoted.

## D87. Geodesic star-convexity prior: implemented, verified, no effect; not adopted
- Idea (D86 write-up): fake tooth regions attached to crowns (2.3% + 3.4% of predicted boundary) should violate a
  shape prior. Star convexity (Veksler, ECCV 2008), geodesic version (Gulshan et al., CVPR 2010): each tooth label
  gets a shortest-path tree from its seed tips; a vertex labelled l needs its parent labelled l.
- core/multilabel: optional star constraints in alpha_expansion. Both directions are infinite pairwise terms that
  are submodular in every move (taking alpha needs the alpha-parent to be alpha: arc parent -> node; keeping l
  forbids the l-parent from switching: arc node -> parent), so moves stay exact. star_repair makes an initial
  labelling feasible (parents first; violators fall back to an unconstrained label). Test vs brute force: feasible
  throughout, energy never rises, no feasible expansion move improves the result; 3 mutants (either arc dropped,
  forbidden nodes ignored) killed.
- PerToothCut: trees from the crease-weighted per-label Dijkstra (star_prior), or from plain edge lengths ("straight
  rays", star_plain). 12-scan slice: weighted -0.001 mm (t -1.0); plain +0.023 (2 scans worse by > 0.1 mm), too
  restrictive for real crowns. 60 sweep scans, weighted: 0.274 -> 0.275 mm, t = +0.4, better 32 / worse 27, at
  ~1.6x the labelling time. **Not adopted** (off by default).
- Why weighted trees change nothing: a region attached through a crease gap is reached by shortest paths through
  that same gap, so it is star-shaped in that metric; the remaining fake regions are mostly around false seeds,
  which are their own centres. The prior constrains shape, not seed correctness.

## D88. Learned per-vertex data term (logistic regression): mixed result, not adopted
- Model: P(tooth | 9 label-free features): log(d_G/d_T) of the crease-weighted distances, height quantile, normal .
  occlusal axis, kappa_min, kappa_max, log(1 + geodesic distance to the cut), log(1 + geodesic distance to the
  nearest tooth seed), squares of height quantile and normal term. Added to the cut as -w A log P (tooth labels) and
  -w A log(1 - P) (gingiva), in both the binary and per-tooth cuts. tools/vertex_train: folds {1,2} train, {0}
  tune w, {3,4} validation; per scan 1,500 random arch vertices + 1,500 within 1.5 mm of the true boundary.
  Unit test: a neutral model (P = 0.5) changes nothing; a confident tooth model only grows tooth regions.
- Validation ROC AUC 0.958 (training 0.956), BELOW the single feature log(d_G/d_T) (0.970): the linear model
  is miscalibrated on the clamped, non-linear distance ratio, and the strongest signal is the one the cut already
  uses. Next strongest: distance to the nearest seed (0.856), height quantile (0.840).
- Sweep on the 60 sweep scans (per-tooth cut 0.274 mm): w 0.25 0.271 (t -0.3, 2 scans worse by > 0.1 mm, up to
  +0.38), w 0.5 0.278, w 1 0.309, w 2 0.330, w 4 0.402. Region metrics improve slightly at w 0.25-0.5 (IoU 0.936 ->
  0.941, per-vertex gingiva IoU 0.922 -> 0.928), boundary metrics do not (F1@0.5 +0.002, HD95 1.94 -> 2.00), and new
  large regressions appear. **Not adopted** (D76/D77 reliability rule). Code stays, off by default; the model file
  is local (data/models, licence question as D76). A binned (additive) model could calibrate better; headroom is
  small because the dominant signal is already in the cut.

## D89. M10a: BVH ray casting (infrastructure for undercut and thickness tools)
- core/bvh: bounding-volume hierarchy built top-down with the binned surface-area heuristic (12 bins over all three
  axes, leaves of <= 4 triangles; MacDonald & Booth 1990, binning after Wald 2007), iterative near-child-first
  traversal with t_max shrinking, Moller-Trumbore ray-triangle tests (1997). Queries: nearest hit (t, face,
  barycentrics) and any-hit occlusion, each with a (t_min, t_max) range and one face to ignore (rays leaving a surface).
- Tests: analytic single triangle (t, barycentrics, miss, parallel, range, ignore); 7,000 random rays vs brute force
  on an icosphere, a torus and a 3,000-triangle random soup (same hit/miss, same t to 1e-12, occlusion consistent);
  rays from inside a closed sphere always hit within (0.99, 1]. 4 mutants (no t_max shrink, lone right child dropped,
  wrong barycentric bound, leaf skips its last triangle) killed.
- Process fix: the mutation harness could test a stale binary when a restore and the next edit fell in the same
  second as the previous build (make compares mtimes). Now every mutant and restore is followed by a 1 s pause and a
  touch; all 11 mutants of D85/D87/D89 were re-run that way and are all killed.

## D90. M10b: undercut map and best path of insertion (virtual surveyor)
- Definition: for a unit withdrawal direction d, a face is undercut if it faces away from d (n . d < -tolerance) or
  is hidden along d (a ray from its centroid towards +d hits the mesh, BVH any-hit; the origin is lifted 1e-7 of the
  bounding diagonal along the face normal, so faces parallel to d, like a vertical wall, are reachable instead of
  grazing their neighbours, which the first version did on 14 equatorial faces of an icosphere). Area-weighted over
  a face region; the whole mesh occludes. Best axis: Fibonacci spiral over the cap within max_tilt of a hint, then a
  pattern search on the sphere (step halving); ties go to the direction nearest the hint.
- Tests with known answers: sphere (undercut = exactly the back-facing faces, fraction 0.5), a plate floating over a
  base (exactly its 2 x 2 footprint, fraction 0.25), a tilted 6-degree-taper frustum (no undercut along its axis;
  the search lands inside the undercut-free cone), a 1-degree taper with only 12 spiral samples (only the local
  refinement reaches the 1-degree cone), and a tilt limit (the answer stays within it). 5 mutants killed (two only
  after adding the last two tests: no local refinement, tilt limit ignored).
- On scans (20 part-1 scans, 263 natural crowns defined by the FDI labels; median per tooth type):

  | tooth type | crowns | undercut along occlusal axis | best axis within 25 deg | tilt | ms per tooth |
  |---|---:|---:|---:|---:|---:|
  | incisors | 79 | 23.8% | 8.0% | 17.4 deg | 739 |
  | canines | 40 | 11.9% | 1.4% | 14.7 deg | 806 |
  | premolars | 78 | 8.1% | 5.0% | 8.2 deg | 2,051 |
  | molars | 66 | 5.5% | 1.1% | 12.0 deg | 4,414 |

  One common path for all crowns of an arch: 13.8% -> 9.3% (tilt 12.2 deg). Incisors carry the most undercut along
  the occlusal axis (labially inclined crowns), consistent with clinical practice of tilting the path for anteriors.
  These are unprepared natural teeth (height of contour, removable-appliance retention), not crown preparations.
- Speed: BVH build 111 ms per scan, ~0.9 M rays/s on one thread (measured while a VM used ~3 cores); a molar search
  takes ~4 s natively, the whole-arch search 13.7 s in the browser. Next: coarse-to-fine (face subsampling for the
  global stage) and threads. Viewer: "Undercut (occlusal axis)" / "Best path of insertion" buttons and an overlay.
