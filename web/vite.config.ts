import { defineConfig } from "vite";

export default defineConfig({
  // Relative asset paths: GitHub Pages serves a project site from /<repo>/, not from the root.
  base: "./",
  // Dev server only (localhost): also serve the git-ignored ../data folder, so local,
  // license-restricted scans can be opened for evaluation (D60). Not part of the built site.
  server: { fs: { allow: [".", "../data"] } },
});
