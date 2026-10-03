import { describe, expect, it } from "vitest";
import { parseDelimited } from "../src/engine/parse";
import text from "../public/example/proteinGroups_example.txt?raw";
import { computeVolcano, defaultParams, lfqColumns } from "../src/engine/volcano";

const table = parseDelimited(text, "\t");
const bait = lfqColumns(table).filter((c) => /BAIT_/.test(c));
const control = lfqColumns(table).filter((c) => /CTRL_/.test(c));
const result = computeVolcano(table, { ...defaultParams, bait, control });

const trueInteractors = ["RLN1A", "RLN1B", "TPC2", "TPC4", "ASM3", "ASM7", "KDR5", "CHP9", "NUF2L", "SRB4", "PMC1", "HELX3"];

describe("synthetic example", () => {
  it("has the expected shape", () => {
    expect(table.rows.length).toBeGreaterThan(500);
    expect(bait).toHaveLength(4);
    expect(control).toHaveLength(4);
    expect(text.length).toBeLessThan(400_000);
    expect(table.rows.every((r) => r[0].startsWith("SYN") || /^(REV|CON)__SYN/.test(r[0]))).toBe(true);
  });

  it("drops reverse and contaminant rows", () => {
    expect(result.dropped.contaminantOrReverse).toBe(16);
    expect(result.proteins.some((p) => p.id.startsWith("REV__") || p.id.startsWith("CON__"))).toBe(false);
  });

  it("flags all 12 true interactors as significant", () => {
    const sig = new Set(result.proteins.filter((p) => p.significant).map((p) => p.gene));
    for (const g of trueInteractors) expect(sig.has(g), g).toBe(true);
  });

  it("keeps background false positives low and borderline binders mixed", () => {
    const known = new Set(["BAITX", ...trueInteractors]);
    const background = result.proteins.filter((p) => !known.has(p.gene) && !p.gene.startsWith("WKB"));
    expect(background.filter((p) => p.significant).length).toBeLessThanOrEqual(5);
    const weak = result.proteins.filter((p) => p.gene.startsWith("WKB"));
    const nSig = weak.filter((p) => p.significant).length;
    expect(weak).toHaveLength(10);
    expect(nSig).toBeGreaterThanOrEqual(1);
    expect(nSig).toBeLessThan(10);
  });
});
