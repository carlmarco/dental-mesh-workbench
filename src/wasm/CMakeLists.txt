# An executable target in Emscripten terms: produces dmw.js + dmw.wasm
# (not a native binary). `dmw` is the module name the JS side imports.
add_executable(dmw bindings.cpp)
target_link_libraries(dmw PRIVATE dmw_core embind)

# Link flags (passed to em++ at link time):
#  MODULARIZE=1       wrap the glue in a factory function: `const m = await createDmw()`
#                     instead of polluting a global `Module` object.
#  EXPORT_ES6=1       emit an ES module (`export default`) so TS/Vite can `import` it.
#  EXPORT_NAME        name of that factory function.
#  ENVIRONMENT        only include loader code for browser and Node (drops
#                     worker/shell variants, which shrinks the glue).
#  ALLOW_MEMORY_GROWTH  linear memory can grow when malloc runs out. Needed for
#                     meshes of unknown size. Cost: JS typed-array views over memory
#                     must be re-created after a grow (M5: the viewer copies at once).
#  EXPORTED_FUNCTIONS / EXPORTED_RUNTIME_METHODS
#                     _malloc/_free and the HEAPU8 byte view, so JS can copy file bytes
#                     into WASM memory and pass C++ an offset (D41).
#  --emit-tsd         generate TypeScript declarations from the embind registry,
#                     replacing the hand-written .d.ts (D8 revisited, D42).
target_link_options(dmw PRIVATE
    -sMODULARIZE=1
    -sEXPORT_ES6=1
    -sEXPORT_NAME=createDmw
    -sENVIRONMENT=web,node
    -sALLOW_MEMORY_GROWTH=1
    -sEXPORTED_FUNCTIONS=_malloc,_free
    -sEXPORTED_RUNTIME_METHODS=HEAPU8
    "SHELL:--emit-tsd ${CMAKE_BINARY_DIR}/dist/dmw.d.ts")

# Put the artifacts at a stable, predictable path for the web package to import.
set_target_properties(dmw PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/dist
    SUFFIX ".mjs")
