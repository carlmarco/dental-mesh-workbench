import { loadDmw } from "./dmw.ts";

const out = document.getElementById("out")!;
const dmw = await loadDmw();
out.textContent = `add(2, 3) = ${dmw.add(2, 3)}  (computed in C++ via WASM)`;
