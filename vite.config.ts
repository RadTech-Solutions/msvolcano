import { defineConfig } from "vite";

// base "./" lets the build run from a GitHub Pages project path or a local file server.
export default defineConfig({ base: "./", test: { include: ["tests/**/*.test.ts"] } });
