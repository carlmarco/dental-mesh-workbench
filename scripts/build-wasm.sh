#!/usr/bin/env bash
# Build the WASM module. Pinned toolchain: emsdk 6.0.11 (see DECISIONS.md D5).
set -euo pipefail

EMSDK_DIR="${EMSDK_DIR:-$HOME/emsdk}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

# emsdk_env.sh only edits PATH/env vars of the *current shell*. Sourcing it here
# (rather than running it) is why emcc is available to the commands below.
# shellcheck disable=SC1091
source "$EMSDK_DIR/emsdk_env.sh" >/dev/null

# emcmake runs `cmake` with -DCMAKE_TOOLCHAIN_FILE=Emscripten.cmake, which swaps
# in em++/emar and sets the EMSCRIPTEN variable our CMakeLists checks.
emcmake cmake -S "$ROOT" -B "$ROOT/build-wasm" -DCMAKE_BUILD_TYPE=Release
cmake --build "$ROOT/build-wasm" -j
