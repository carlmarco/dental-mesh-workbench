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
| 9 | Dental scans (Teeth3DS+, local only, D60): tooth-gingiva margin and cusp/ridge detection by feature-line extraction, measured against Teeth3DS labels and 3DTeethLand cusp landmarks; scan-size optimization; handle-loop scan QA (8b) | Meyer et al. 2003 + curvature-tensor / ridge literature (to be verified) | In progress: QA, handle loops, cusp tips, margin, failure analysis (D74-D77) done; min-cut margin labelling or signed heat (8) next |
| 8 | Stretch, only after 7 ships: signed heat method (implicit repair of open/noisy surfaces) | Feng & Crane 2024 | |
| 8b | Homology handle loops via tree-cotree decomposition (greedy shortest system) | Eppstein 2003; Erickson & Whittlesey 2005 | Done |

**Out of scope:** restoration generation, remeshing, explicit mesh repair (P6), WebGPU compute. Segmentation by
geometric feature lines is in scope (M9); ML segmentation is allowed if it becomes the obvious next step (D61).
