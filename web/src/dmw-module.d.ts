// The Emscripten glue (src/wasm/dmw.mjs) is generated JS with no types. This
// declares its shape for the type checker. It must be kept in sync by hand with
// src/wasm/bindings.cpp: anything registered in EMSCRIPTEN_BINDINGS appears here.
// (Relative paths aren't allowed in `declare module`, hence the wildcard.)
declare module "*/dmw.mjs" {
  export interface DmwModule {
    add(a: number, b: number): number;
  }
  // MODULARIZE=1 + EXPORT_ES6=1 + EXPORT_NAME=createDmw: default export is an
  // async factory that fetches/compiles the .wasm and resolves to the module.
  const createDmw: () => Promise<DmwModule>;
  export default createDmw;
}
