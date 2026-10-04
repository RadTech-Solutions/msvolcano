import { describe, expect, it } from "vitest";
import { computeStoichiometry, countTrypticPeptides, parsePeptideTable } from "../src/engine/stoich";
import type { Protein } from "../src/engine/volcano";

const a = (n: number) => "A".repeat(n);

describe("countTrypticPeptides", () => {
  it("counts simple cleavage after K and R", () => {
    expect(countTrypticPeptides(a(6) + "K" + a(6) + "R")).toBe(2);
  });
  it("does not cleave before proline", () => {
    // K followed by P stays joined: one peptide of 15 residues
    expect(countTrypticPeptides(a(6) + "KP" + a(6) + "K")).toBe(1);
    // R before P too, but a later cleavage site still works
    expect(countTrypticPeptides(a(6) + "RP" + a(5) + "K")).toBe(1);
  });
  it("applies the length limits at the boundaries", () => {
    expect(countTrypticPeptides(a(5) + "K")).toBe(0); // 6
    expect(countTrypticPeptides(a(6) + "K")).toBe(1); // 7
    expect(countTrypticPeptides(a(29) + "K")).toBe(1); // 30
    expect(countTrypticPeptides(a(30) + "K")).toBe(0); // 31
  });
  it("handles sequences without a cleavage site", () => {
    expect(countTrypticPeptides(a(10))).toBe(1);
    expect(countTrypticPeptides(a(4))).toBe(0);
    expect(countTrypticPeptides("")).toBe(0);
  });
  it("counts the C terminal peptide and drops short tails", () => {
    expect(countTrypticPeptides(a(7) + "K" + a(3))).toBe(1);
    expect(countTrypticPeptides(a(7) + "K" + a(7))).toBe(2);
  });
  it("treats a terminal K or R followed by nothing as a normal end", () => {
    expect(countTrypticPeptides(a(7) + "KK")).toBe(1); // 8 and 1
  });
});

describe("parsePeptideTable", () => {
  it("parses rows and skips the header and blanks", () => {
    const m = parsePeptideTable("accession\tgene\ttryptic_peptides\nP1\tAbc\t12\nP2\t\t3\n\n");
    expect(m.size).toBe(2);
    expect(m.get("P1")).toEqual({ gene: "Abc", tryptic: 12 });
    expect(m.get("P2")).toEqual({ gene: "", tryptic: 3 });
  });
});

const prot = (row: number, gene: string, id: string, x: number, baitMean: number): Protein => ({
  row, gene, id, x, y: 0, baitMean, se: 0, values: [], raw: [], imputed: [], significant: false,
});

describe("computeStoichiometry", () => {
  // Ip = 2^baitMean - 2^(baitMean - x). With baitMean 10, x 1: 1024 - 512 = 512.
  // baitMean 8, x 1: 256 - 128 = 128.
  const peptides = new Map([
    ["BAIT", { gene: "Bait", tryptic: 16 }],
    ["PREY", { gene: "Prey", tryptic: 4 }],
    ["ISO", { gene: "Iso", tryptic: 8 }],
  ]);
  it("computes exact ratios", () => {
    const proteins = [prot(0, "Bait", "BAIT", 1, 10), prot(1, "Prey", "PREY", 1, 8), prot(2, "Iso", "ISO-2", 1, 8)];
    const st = computeStoichiometry(proteins, peptides, "Bait");
    expect(st.get(0)).toBe(1);
    // (128/4) / (512/16) = 1
    expect(st.get(1)).toBe(1);
    // (128/8) / 32 = 0.5, isoform suffix stripped
    expect(st.get(2)).toBe(0.5);
  });
  it("finds the bait by accession and uses 0.01 for unknown proteins", () => {
    const proteins = [prot(0, "Bait", "BAIT", 1, 10), prot(1, "X", "NOPE", 1, 8)];
    const st = computeStoichiometry(proteins, peptides, "BAIT");
    expect(st.get(1)).toBeCloseTo(128 / 0.01 / 32, 6);
  });
  it("returns NaN for proteins not enriched over control", () => {
    const proteins = [prot(0, "Bait", "BAIT", 1, 10), prot(1, "Prey", "PREY", -1, 8)];
    expect(computeStoichiometry(proteins, peptides, "Bait").get(1)).toBeNaN();
  });
  it("throws for a missing or non enriched bait", () => {
    const proteins = [prot(0, "Bait", "BAIT", -1, 10), prot(1, "Prey", "PREY", 1, 8)];
    expect(() => computeStoichiometry(proteins, peptides, "Nothing")).toThrow(/not found/);
    expect(() => computeStoichiometry(proteins, peptides, "Bait")).toThrow(/non positive/);
  });
});
