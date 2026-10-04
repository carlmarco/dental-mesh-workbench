cmake_minimum_required(VERSION 3.24)
project(dental_mesh_workbench LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)  # -std=c++20, not gnu++20: stay portable

# Emit compile_commands.json so editors/clangd see the real flags.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_subdirectory(src/core)

# `EMSCRIPTEN` is set by emcmake's toolchain file. Native-only: tests (Catch2
# can't run in a browser). WASM-only: the bindings target.
if(EMSCRIPTEN)
    add_subdirectory(src/wasm)
else()
    include(FetchContent)
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.16.0  # pinned: see DECISIONS.md
    )
    FetchContent_MakeAvailable(Catch2)
    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)

    enable_testing()  # must be in the top-level file for `ctest` to find tests
    add_subdirectory(tests)
endif()
