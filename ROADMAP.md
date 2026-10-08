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
| 9 | Dental scans (Teeth3DS+, local only, D60): tooth-gingiva margin and cusp/ridge detection by feature-line extraction, measured against Teeth3DS labels and 3DTeethLand cusp landmarks; scan-size optimization; handle-loop scan QA (8b) | Meyer et al. 2003 + curvature-tensor / ridge literature (to be verified) | In progress: QA, handle loops, cusp tips, margin, failure analysis (D74-D77), graph-cut labelling (D78), error analysis (D80-D84), per-tooth labelling (D85) done; publication pass, manufacturing-aware inspection or signed heat (8) next |
| 8 | Stretch, only after 7 ships: signed heat method (implicit repair of open/noisy surfaces) | Feng & Crane 2024 | |
| 8b | Homology handle loops via tree-cotree decomposition (greedy shortest system) | Eppstein 2003; Erickson & Whittlesey 2005 | Done |

## M10 plan: manufacturing-aware inspection and geometry generation

Ordered by value for effort; each item gets synthetic shapes with analytic answers as tests, then a measured
result on scans, like M9. Shared infrastructure first.

| Step | Tool | What it answers | Method (to cite and verify) | Test with a known answer |
|---|---|---|---|---|
| 10a | BVH ray casting | infrastructure for 10b-10c | SAH bounding-volume hierarchy, Möller-Trumbore ray-triangle | brute-force intersection on random rays |
| 10b | Insertion axis and undercut map | can the part be seated/removed or milled along a direction; which surface is undercut; the best path of draw | per-face visibility along d (n . d and occlusion by ray cast); minimize undercut area over directions on the hemisphere (sampling + local refinement) | sphere: undercut = everything below the equator for vertical d; a tilted cylinder: best axis = its own axis |
| 10c | Wall thickness map | minimum thickness for printing/milling (e.g. crowns, shells, printed models) | shape diameter function (Shapira, Shamir, Cohen-Or 2008): inward cone of rays, robust statistics | spherical shell: thickness = R_out - R_in; plate: its thickness |
| 10d | Signed distance on open scans (M8) | an implicit surface for offsets, without repairing holes first | signed heat method (Feng & Crane 2024) on a background grid | sphere with a hole: signed distance vs analytic |
| 10e | Offset surfaces: cement gap, hollowing for printing | GENERATES geometry: offset shells, hollowed models with a wall thickness | iso-surface of the signed distance (marching cubes; Lorensen & Cline 1987) | sphere offset by d: radius r + d within grid error, converging under refinement |
| 10f | Hole filling and remeshing | the mesh-editing work of many CAD roles | hole filling (Liepa 2003), isotropic remeshing (Botsch & Kobbelt 2004) | Euler characteristic after filling; edge-length statistics after remeshing |

10b and 10c are inspection (fast, directly dental CAD/CAM); 10d-10e are generation (the restoration and
manufacturing side); 10f is breadth. "Out of scope" below is lifted for 10e-10f by this plan.

**Out of scope:** restoration generation, remeshing, explicit mesh repair (P6), WebGPU compute. Segmentation by
geometric feature lines is in scope (M9); ML segmentation is allowed if it becomes the obvious next step (D61).
