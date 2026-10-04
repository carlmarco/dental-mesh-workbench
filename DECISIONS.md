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
  externalized for the browser. The build succeeds. Not yet verified at runtime in a browser.
  Candidate fix: build a web-only variant for Vite and keep the node variant for tests.
