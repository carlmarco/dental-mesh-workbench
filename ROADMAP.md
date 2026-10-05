# Roadmap

One milestone at a time, in order. Each ends with a summary, check questions and a commit.

| # | Milestone | Paper | Status |
|---|---|---|---|
| 1 | Repo and toolchain: native build + tests, Emscripten module callable from TypeScript | | Done |
| 2 | Mesh I/O (OBJ, STL), vertex welding, half-edge construction | | Done |
| 3 | Topology diagnostics: boundary / non-manifold edges, components, Euler characteristic, genus and Betti numbers per component (D22) | | Done |
| 4 | Discrete curvature: mixed Voronoi areas, cotan mean-curvature normal, angle-defect Gaussian curvature; Gauss-Bonnet and sphere-convergence tests | Meyer, Desbrun, Schröder, Barr 2003 | Done |
| 5 | Viewer: WASM + Three.js, overlays for topology and curvature | | Done |
| 6 | Heat-method geodesics, operators written in DEC form (L = d0ᵀ ⋆1 d0); validated vs great-circle distance and Dijkstra | Crane, Weischedel, Wardetzky 2013 | Done |
| 6b | Intrinsic Delaunay cotan Laplacian (non-negative weights; addresses the obtuse-triangle caveat of 6) | Bobenko & Springborn 2007; flips: Fisher et al. 2007 (Sharp & Crane 2020 for the non-manifold extension, not implemented) | Done |
| 7 | Benchmarks and README | | Done |
| 8 | Stretch, only after 7 ships: signed heat method (implicit repair of open/noisy surfaces) | Feng & Crane 2024 | |
| 8b | Optional: homology handle loops via tree-cotree decomposition | Eppstein 2003 (to be read) | |

**Out of scope:** segmentation, restoration generation, remeshing, explicit mesh repair (P6), WebGPU compute,
machine learning.
