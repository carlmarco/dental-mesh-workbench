// Headless check that TS -> WASM -> C++ works. Run: npm run smoke
// Mirrors tests/test_smoke.cpp so native and WASM agree.
import assert from "node:assert/strict";
import { loadDmw } from "./dmw.ts";

const dmw = await loadDmw();
assert.equal(dmw.add(2, 3), 5);
assert.equal(dmw.add(-4, 4), 0);
assert.equal(dmw.add(0, 0), 0);
console.log("smoke ok: add() callable from TypeScript via WASM");
