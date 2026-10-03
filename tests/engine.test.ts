import { describe, expect, it } from "vitest";
import { parseDelimited } from "../src/engine/parse";
import { tTest, tTwoSidedP } from "../src/engine/stats";
import { computeVolcano, cutoffCurve, isSignificant } from "../src/engine/volcano";

describe("t distribution", () => {
  // Reference values from scipy: 2 * t.sf(2, 10)
  it("matches scipy two sided p values", () => {
    expect(tTwoSidedP(2, 10)).toBeCloseTo(0.07338803, 7);
    expect(tTwoSidedP(1, 1)).toBeCloseTo(0.5, 10);
    expect(tTwoSidedP(0, 7)).toBeCloseTo(1, 10);
  });
});

describe("tTest", () => {
  const a = [5.1, 4.9, 5.6, 5.8, 6.0];
  const b = [4.0, 4.2, 3.9, 4.4, 4.1];
  // Reference values from scipy.stats.ttest_ind
  it("Student matches scipy", () => {
    const r = tTest(a, b, false)!;
    expect(r.df).toBe(8);
    expect(r.t).toBeCloseTo(6.03402426, 7);
    expect(r.p).toBeCloseTo(0.000311383454, 9);
    expect(r.diff).toBeCloseTo(1.36, 10);
  });
  it("Welch matches scipy", () => {
    const r = tTest(a, b, true)!;
    expect(r.df).toBeCloseTo(5.325519006, 6);
    expect(r.p).toBeCloseTo(0.00144174006, 8);
  });
  it("returns null for zero variance or n < 2", () => {
    expect(tTest([1, 1], [1, 1], false)).toBeNull();
    expect(tTest([1], [1, 2], false)).toBeNull();
  });
});

describe("cutoff", () => {
  it("marks points beyond the hyperbola", () => {
    expect(isSignificant(5, 10, 3, 3)).toBe(true); // 3 / 2 = 1.5 < 10
    expect(isSignificant(3.1, 2, 3, 3)).toBe(false); // 3 / 0.1 = 30
    expect(isSignificant(-6, 10, 3, 3)).toBe(false); // right side only
  });
  it("curve starts where it meets the top of the plot", () => {
    const { right } = cutoffCurve(3, 3, 15, 5);
    expect(right[0].x).toBeCloseTo((3 + 5 * 3) / 5, 10);
    expect(right[0].y).toBeCloseTo(5, 6);
  });
});

describe("computeVolcano", () => {
  const tsv = [
    "Majority protein IDs\tGene names\tReverse\tPotential contaminant\tLFQ intensity B1\tLFQ intensity B2\tLFQ intensity B3\tLFQ intensity C1\tLFQ intensity C2\tLFQ intensity C3",
    "P1\tBAIT\t\t\t30\t31\t30.5\t20\t20.2\t19.8",
    "P2\tPREY\t\t\t28\t29\t28.5\t18\t17.5\t18.2",
    "P3\tFLAT\t\t\t24\t24.5\t23.8\t24.1\t23.9\t24.2",
    "P4\tREV\t+\t\t30\t30\t30\t10\t10\t10",
    "P5\tCONT\t\t+\t30\t30\t30\t10\t10\t10",
    "P6\tMISSING\t\t\t\t\t\t\t\t",
  ].join("\n");
  const table = parseDelimited(tsv, "\t");
  const bait = table.columns.filter((c) => /B\d$/.test(c));
  const control = table.columns.filter((c) => /C\d$/.test(c));

  it("filters reverse and contaminants and calls enriched proteins", () => {
    const r = computeVolcano(table, { bait, control, welch: false, shift: 1.8, shrink: 0.3, minFoldChange: 3, curvature: 3, seed: 1 });
    expect(r.dropped.contaminantOrReverse).toBe(2);
    expect(r.logTransformed).toBe(false);
    const by = Object.fromEntries(r.proteins.map((p) => [p.gene, p]));
    expect(by.BAIT.significant).toBe(true);
    expect(by.PREY.significant).toBe(true);
    expect(by.FLAT.significant).toBe(false);
    expect(by.BAIT.x).toBeGreaterThan(9);
  });
  it("is reproducible with a fixed seed", () => {
    const run = () => computeVolcano(table, { bait, control, welch: true, shift: 1.8, shrink: 0.3, minFoldChange: 3, curvature: 3, seed: 7 });
    expect(run().proteins).toEqual(run().proteins);
  });
  it("requires two columns per group", () => {
    expect(() => computeVolcano(table, { bait: [bait[0]], control, welch: false, shift: 1.8, shrink: 0.3, minFoldChange: 3, curvature: 3, seed: 1 })).toThrow();
  });
});
