import { defineConfig } from "vite";

// Relative asset paths: GitHub Pages serves a project site from /<repo>/, not from the root.
export default defineConfig({ base: "./" });
