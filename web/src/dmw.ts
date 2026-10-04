import createDmw from "./wasm/dmw.mjs";
import type { DmwModule } from "./wasm/dmw.mjs";

export type Dmw = DmwModule;

// Instantiation is async: the .wasm must be fetched and compiled before any
// C++ can run. Callers await this once and then call synchronously.
export function loadDmw(): Promise<Dmw> {
  return createDmw();
}
