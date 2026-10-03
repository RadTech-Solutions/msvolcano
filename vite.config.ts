import { readFileSync } from "node:fs";
import { defineConfig } from "vite";

const pkg = JSON.parse(readFileSync(new URL("./package.json", import.meta.url), "utf8"));

// base "./" lets the build run from a GitHub Pages project path or a local file server.
export default defineConfig({
  base: "./",
  define: { __APP_VERSION__: JSON.stringify(pkg.version) },
  test: { include: ["tests/**/*.test.ts"] },
});
