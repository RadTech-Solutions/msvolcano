import { describe, expect, it } from "vitest";
import { parseDelimited } from "../src/engine/parse";
import { parseNumber, usesDecimalComma } from "../src/engine/numbers";
import { permutationFdr } from "../src/engine/experimental";
import { sniffSeparator } from "../src/engine/parse";
import { baitRecovery, coverageCheck, loadingCheck, observedCorrelation, samplePca, symmetricEigen } from "../src/engine/qc";
import { tTest } from "../src/engine/stats";
import { computeVolcano, presenceCsv, shiftSensitivity, type Params, type Protein } from "../src/engine/volcano";

const HEAD = "Majority protein IDs\tGene names\tReverse\tPotential contaminant\tLFQ intensity B_1\tLFQ intensity B_2\tLFQ intensity B_3\tLFQ intensity B_4\tLFQ intensity C_1\tLFQ intensity C_2\tLFQ intensity C_3\tLFQ intensity C_4";
const row = (id: string, gene: string, b: number[], c: number[]) => [id, gene, "", "", ...b, ...c].join("\t");
const base: Omit<Params, "bait" | "control"> = { welch: false, shift: 1.8, shrink: 0.3, minFoldChange: 3, curvature: 3, seed: 1 };

// 40 flat proteins plus the cases under test
const flat = Array.from({ length: 40 }, (_, i) => row(`F${i}`, `FLAT${i}`, [2e7 * (1 + (i % 5) / 10), 2.1e7, 1.9e7, 2.05e7], [2e7, 2.1e7, 1.95e7, 2.02e7]));
const cases = [
  row("P_ONLY", "ONLYBAIT", [5e7, 6e7, 5.5e7, 0], [0, 0, 0, 0]), // bait 3 of 4, no control
  row("P_SPARSE", "SPARSE", [5e7, 0, 0, 0], [0, 0, 0, 0]), // bait 1 of 4, no control
  row("P_MIXED", "MIXED", [5e7, 6e7, 5.5e7, 5.2e7], [1e6, 0, 1.2e6, 1.1e6]), // one control missing
  row("P_NOBAIT", "NOBAIT", [0, 0, 0, 0], [3e7, 3e7, 3e7, 3e7]),
];
const table = parseDelimited([HEAD, ...flat, ...cases].join("\n"), "\t");
const bait = table.columns.filter((c) => /_B_|B_\d/.test(c) && c.includes("B_"));
const control = table.columns.filter((c) => c.includes("C_"));
const params: Params = { ...base, bait, control };

describe("presence only handling", () => {
  it("as published keeps absent in control proteins in the test with imputed values", () => {
    const r = computeVolcano(table, { ...params, presence: "include" });
    expect(r.presenceOnly).toHaveLength(0);
    const q = r.proteins.find((p) => p.gene === "ONLYBAIT")!;
    expect(q.imputed.slice(4).every(Boolean)).toBe(true); // all control values were filled in
    expect(q.imputed[3]).toBe(true);
  });
  it("separate mode moves them out of the test and lists them with counts", () => {
    const r = computeVolcano(table, { ...params, presence: "separate" });
    expect(r.proteins.find((p) => p.gene === "ONLYBAIT")).toBeUndefined();
    expect(r.presenceOnly.map((p) => p.gene)).toEqual(["ONLYBAIT"]);
    expect(r.presenceOnly[0].nObs).toBe(3);
    expect(r.presenceOnly[0].nBait).toBe(4);
    expect(r.dropped.presenceOnly).toBe(1);
  });
  it("leaves out a protein measured only once in total, which would be tested on filled-in values alone", () => {
    const r = computeVolcano(table, { ...params, presence: "separate" });
    expect(r.proteins.find((p) => p.gene === "SPARSE")).toBeUndefined();
    expect(r.presenceOnly.find((p) => p.gene === "SPARSE")).toBeUndefined();
    expect(r.dropped.tooFewValues).toBe(1);
  });
  it("as published, the same protein is still tested", () => {
    const r = computeVolcano(table, { ...params, presence: "include" });
    expect(r.proteins.find((p) => p.gene === "SPARSE")).toBeDefined();
  });
  it("still tests a protein with some control values", () => {
    const r = computeVolcano(table, { ...params, presence: "separate" });
    const q = r.proteins.find((p) => p.gene === "MIXED")!;
    expect(q.imputed.filter(Boolean)).toHaveLength(1);
  });
  it("drops proteins that are absent in bait", () => {
    const r = computeVolcano(table, { ...params, presence: "separate" });
    expect(r.dropped.absentInBait).toBe(1);
    expect(r.proteins.find((p) => p.gene === "NOBAIT")).toBeUndefined();
  });
  it("writes a CSV of the presence only list", () => {
    const r = computeVolcano(table, { ...params, presence: "separate" });
    const csv = presenceCsv(table, r.presenceOnly, params).split("\n");
    expect(csv[0].startsWith("gene,protein_id,bait_replicates_with_value")).toBe(true);
    expect(csv[1].startsWith("ONLYBAIT,P_ONLY,3,4")).toBe(true);
  });
});

describe("no imputation mode", () => {
  it("tests only the measured values and matches a direct t-test", () => {
    const r = computeVolcano(table, { ...params, imputation: "none", presence: "include" });
    const q = r.proteins.find((p) => p.gene === "MIXED")!;
    const rawB = q.raw.slice(0, 4).filter((v) => !Number.isNaN(v));
    const rawC = q.raw.slice(4).filter((v) => !Number.isNaN(v));
    const t = tTest(rawB, rawC, false)!;
    expect(q.x).toBeCloseTo(t.diff, 10);
    expect(q.y).toBeCloseTo(-Math.log10(t.p), 10);
    expect(q.imputed.some(Boolean)).toBe(true); // flags still show what was missing
  });
  it("excludes proteins with fewer than two measured values in a group", () => {
    const r = computeVolcano(table, { ...params, imputation: "none", presence: "include" });
    expect(r.proteins.find((p) => p.gene === "ONLYBAIT")).toBeUndefined();
    expect(r.dropped.untestable).toBeGreaterThan(0);
  });
  it("has no sensitivity to report", () => {
    expect(shiftSensitivity(table, { ...params, imputation: "none" }, [1.5, 1.8])).toBeNull();
  });
});

describe("shift sensitivity", () => {
  it("reports hit counts and overlap with the current setting", () => {
    const rows = shiftSensitivity(table, { ...params, presence: "include" }, [1.2, 1.8, 2.4])!;
    expect(rows).toHaveLength(3);
    const b = rows.find((x) => x.base)!;
    expect(b.shared).toBe(b.hits);
    expect(b.lost).toBe(0);
    for (const x of rows) expect(x.shared + x.gained).toBe(x.hits);
  });
});

describe("quality checks", () => {
  const r = computeVolcano(table, { ...params, presence: "separate" });
  it("computes correlation from measured values only", () => {
    const c = observedCorrelation(r.proteins);
    expect(c).toHaveLength(8);
    expect(c[0][0]).toBe(1);
    expect(c[0][1]).toBeCloseTo(c[1][0], 12);
  });
  it("finds the bait among enriched proteins", () => {
    expect(baitRecovery(r.proteins, r.presenceOnly, "ONLYBAIT").status).toBe("ok");
    expect(baitRecovery(r.proteins, r.presenceOnly, "MIXED").status).toBe("ok");
    // 5 of 40 or so flat proteins are within the top group, so pick the weakest one to see a warning
    const weakest = [...r.proteins].sort((a, b) => a.x - b.x)[0];
    expect(baitRecovery(r.proteins, r.presenceOnly, weakest.gene).status).toBe("warn");
    // rank 4 of many is fine
    const fourth = [...r.proteins].sort((a, b) => b.x - a.x)[3];
    expect(baitRecovery(r.proteins, r.presenceOnly, fourth.gene).status).toBe("ok");
    expect(baitRecovery(r.proteins, r.presenceOnly, "NOPE").found).toBe(false);
    expect(baitRecovery(r.proteins, r.presenceOnly, "").status).toBe("info");
  });
  it("warns about unequal loading", () => {
    const shifted: Protein[] = r.proteins.map((q) => ({ ...q, raw: q.raw.map((v, j) => (j < 4 ? v + 2 : v)) }));
    expect(loadingCheck(shifted, 4).status).toBe("warn");
    expect(loadingCheck(r.proteins, 4).status).toBe("ok");
  });
  it("warns when one sample lags the others of its group, not when a whole group is lower", () => {
    const lag = [
      { column: "b1", observed: 90, total: 100, group: "bait" as const }, { column: "b2", observed: 92, total: 100, group: "bait" as const }, { column: "b3", observed: 20, total: 100, group: "bait" as const },
      { column: "c1", observed: 90, total: 100, group: "control" as const },
    ];
    expect(coverageCheck(lag).status).toBe("warn");
    expect(coverageCheck(lag).detail).toContain("b3");
    const controlsLower = [
      { column: "b1", observed: 95, total: 100, group: "bait" as const }, { column: "b2", observed: 96, total: 100, group: "bait" as const },
      { column: "c1", observed: 45, total: 100, group: "control" as const }, { column: "c2", observed: 47, total: 100, group: "control" as const },
    ];
    expect(coverageCheck(controlsLower).status).toBe("ok");
  });
});

describe("pca", () => {
  it("diagonalises a small symmetric matrix", () => {
    const { values } = symmetricEigen([[2, 1], [1, 2]]);
    expect(values[0]).toBeCloseTo(3, 10);
    expect(values[1]).toBeCloseTo(1, 10);
  });
  it("separates two clear groups on PC1", () => {
    const mk = (g: number, i: number): Protein["raw"] => [0, 1, 2, 3, 4, 5].map((s) => (s < 3 ? 10 : 10 + g) + i * 0.3 + (s % 3) * 0.05);
    const prots = Array.from({ length: 30 }, (_, i) => ({ row: i, gene: "g", id: "i", x: 0, y: 0, baitMean: 0, se: 0, values: mk(2 + (i % 3), i), raw: mk(2 + (i % 3), i), imputed: new Array(6).fill(false), significant: false })) as Protein[];
    const res = samplePca(prots)!;
    expect(res.proteinsUsed).toBe(30);
    expect(res.explained[0]).toBeGreaterThan(0.8);
    const left = res.scores.slice(0, 3).map((s) => s[0]), right = res.scores.slice(3).map((s) => s[0]);
    expect(Math.max(...left) * Math.min(...right) < 0 || Math.min(...left) * Math.max(...right) < 0).toBe(true);
  });
  it("returns null with too few complete proteins", () => {
    const r = computeVolcano(table, { ...params, presence: "separate" });
    expect(samplePca(r.proteins.slice(0, 3))).toBeNull();
  });
});


describe("audit fixes", () => {
  it("refuses a selected column with no values instead of mis-detecting the scale", () => {
    const t = parseDelimited([HEAD, row("A", "GA", [1e7, 1e7, 1e7, 0], [1e7, 1e7, 1e7, 1e7]), row("B", "GB", [2e7, 2e7, 2e7, 0], [1e7, 1e7, 1e7, 1e7])].join("\n"), "\t");
    expect(() => computeVolcano(t, { ...params })).toThrow(/has no values/);
  });
  it("detects the scale from the typical value and honours the user's choice", () => {
    expect(computeVolcano(table, { ...params }).scale.mode).toBe("raw");
    const logged = parseDelimited([HEAD, ...Array.from({ length: 30 }, (_, i) => row(`L${i}`, `L${i}`, [25 + i / 10, 25.2, 24.9, 25.1], [24, 24.1, 23.9, 24.2]))].join("\n"), "\t");
    const auto = computeVolcano(logged, { ...params, bait, control });
    expect(auto.scale.mode).toBe("log2");
    expect(auto.logTransformed).toBe(false);
    expect(computeVolcano(logged, { ...params, bait, control, logMode: "raw" }).scale.mode).toBe("raw");
  });
  it("warns when log2 is forced on raw looking values", () => {
    expect(computeVolcano(table, { ...params, logMode: "log2" }).scale.warning).toContain("raw intensities");
  });
  it("warns when values look like log10 or natural log", () => {
    const small = parseDelimited([HEAD, ...Array.from({ length: 30 }, (_, i) => row(`L${i}`, `L${i}`, [6 + i / 100, 6.1, 5.9, 6.05], [5, 5.1, 4.9, 5.2]))].join("\n"), "\t");
    expect(computeVolcano(small, { ...params, bait, control }).scale.warning).toContain("log10");
  });
  it("gives every filled-in cell its own draw, so changing the presence setting does not reshuffle other proteins", () => {
    const a = computeVolcano(table, { ...params, presence: "include" }).proteins.find((p) => p.gene === "MIXED")!;
    const b = computeVolcano(table, { ...params, presence: "separate" }).proteins.find((p) => p.gene === "MIXED")!;
    expect(a.values).toEqual(b.values);
    expect(a.x).toBe(b.x);
  });
  it("is reproducible for the same seed and differs for another", () => {
    const f = (seed: number) => computeVolcano(table, { ...params, seed }).proteins.find((p) => p.gene === "MIXED")!.values;
    expect(f(1)).toEqual(f(1));
    expect(f(1)).not.toEqual(f(2));
  });
  it("returns no test when the variance is only rounding noise", () => {
    expect(tTest([0.1, 0.1, 0.1], [0.3, 0.3, 0.3], false)).toBeNull();
    const t = tTest([1, 2, 3], [2, 3, 5], false)!;
    expect(t.se).toBeCloseTo(t.diff / t.t, 12);
  });
  it("cannot run permutation FDR when values are missing", () => {
    const r = computeVolcano(table, { ...params, imputation: "none", presence: "include" });
    expect(permutationFdr(r.proteins, 4, 0.1, 0.05, 20, 1)).toBeNull();
  });
  it("needs 20 shared proteins before reporting a correlation", () => {
    const few = computeVolcano(table, { ...params }).proteins.slice(0, 10);
    expect(Number.isNaN(observedCorrelation(few)[0][1])).toBe(true);
  });
});

describe("numbers and separators", () => {
  it("reads thousands separators and decimal commas correctly", () => {
    expect(parseNumber("1,234.5")).toBe(1234.5);
    expect(parseNumber("1.234,5", true)).toBe(1234.5);
    expect(parseNumber("12,5", true)).toBe(12.5);
    expect(parseNumber("1,234,567")).toBe(1234567);
    expect(parseNumber("n.d.")).toBeNaN();
    expect(parseNumber("N/A")).toBeNaN();
    expect(parseNumber("abc")).toBeNaN();
  });
  it("decides the style from the values", () => {
    expect(usesDecimalComma(["1,5", "2,25", "3"])).toBe(true);
    expect(usesDecimalComma(["1.5", "2,25"])).toBe(false);
    expect(usesDecimalComma(["1,234", "2,345,678"])).toBe(false);
    expect(usesDecimalComma(["24,512", "30,118", ""])).toBe(true);
    expect(usesDecimalComma(["1,234.5", "2,345.6"])).toBe(false);
  });
  it("sniffs tab, semicolon and comma", () => {
    expect(sniffSeparator("a\tb\tc\n1\t2\t3")).toBe("\t");
    expect(sniffSeparator("a;b;c\n1,5;2,5;3,5")).toBe(";");
    expect(sniffSeparator("a,b,c\n1,2,3")).toBe(",");
  });
});
