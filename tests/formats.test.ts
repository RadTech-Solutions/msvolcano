import { describe, expect, it } from "vitest";
import { detectFormat, explainUnsupported } from "../src/engine/formats";
import { parseDelimited } from "../src/engine/parse";
import { computeVolcano } from "../src/engine/volcano";
import { groupColumns, suggestRole } from "../src/groups";

const tsv = (...lines: string[]) => parseDelimited(lines.join("\n"), "\t");

describe("MaxQuant", () => {
  const t = tsv(
    "Protein IDs\tMajority protein IDs\tGene names\tOnly identified by site\tReverse\tPotential contaminant\tIntensity\tLFQ intensity A_1\tLFQ intensity A_2\tLFQ intensity B_1\tLFQ intensity B_2\tiBAQ A_1",
    "P1;P2\tP1;P2\tG1\t\t\t\t9e9\t1e7\t2e7\t3e7\t4e7\t5",
    "P3\tP3\tG3\t\t+\t\t1e9\t1e7\t2e7\t3e7\t4e7\t5",
  );
  it("takes only the LFQ columns, not iBAQ or the total intensity", () => {
    const f = detectFormat(t)!;
    expect(f.id).toBe("maxquant");
    expect(f.sampleColumns).toEqual(["LFQ intensity A_1", "LFQ intensity A_2", "LFQ intensity B_1", "LFQ intensity B_2"]);
    expect(f.isFlagged(t.rows[1])).toBe(true);
    expect(f.isFlagged(t.rows[0])).toBe(false);
  });
});

describe("FragPipe", () => {
  const t = tsv(
    "Protein\tProtein ID\tEntry Name\tGene\tProtein Length\tCoverage\tOrganism\tProtein Existence\tDescription\tProtein Probability\tTop Peptide Probability\tCombined Total Peptides\tCombined Spectral Count\tCombined Unique Spectral Count\tCombined Total Spectral Count\tbait_1 Spectral Count\tbait_1 Intensity\tbait_1 MaxLFQ Intensity\tbait_2 MaxLFQ Intensity\tctrl_1 MaxLFQ Intensity\tctrl_2 MaxLFQ Intensity\tIndistinguishable Proteins",
    "sp|P12345|AAA_HUMAN\tP12345\tAAA_HUMAN\tAAA\t100\t10\tHuman\t1\tdesc\t1\t1\t5\t5\t5\t5\t3\t1e6\t2e6\t3e6\t1e5\t2e5\t",
    "contam_sp|P99999|KRT1\tP99999\tKRT1\tKRT1\t100\t10\tHuman\t1\tdesc\t1\t1\t5\t5\t5\t5\t3\t1e6\t2e6\t3e6\t1e5\t2e5\t",
  );
  it("uses MaxLFQ columns only and flags contaminants", () => {
    const f = detectFormat(t)!;
    expect(f.id).toBe("fragpipe");
    expect(f.sampleColumns).toEqual(["bait_1 MaxLFQ Intensity", "bait_2 MaxLFQ Intensity", "ctrl_1 MaxLFQ Intensity", "ctrl_2 MaxLFQ Intensity"]);
    expect(f.clean("bait_1 MaxLFQ Intensity")).toBe("bait_1");
    expect(f.isFlagged(t.rows[1])).toBe(true);
    expect(f.isFlagged(t.rows[0])).toBe(false);
    expect(groupColumns(f.sampleColumns, f.clean).map((g) => g.name)).toEqual(["bait", "ctrl"]);
  });
});

describe("DIA-NN", () => {
  const t = tsv(
    "Protein.Group\tProtein.Names\tGenes\tFirst.Protein.Description\tD:\\data\\Bait_1.raw\tD:\\data\\Bait_2.raw\t/mnt/Ctrl_1.d\t/mnt/Ctrl_2.d",
    "Q9Y6Y8;A0A994\tS23IP_HUMAN;X_HUMAN\tSEC23IP;X\tSEC23-interacting protein\t1.5e6\t1.6e6\t\t2e5",
    "Cont_P00761\tTRYP_PIG\t\ttrypsin\t5e5\t6e5\t5e5\t5e5",
  );
  it("takes the run columns, ignores metadata, cleans file paths", () => {
    const f = detectFormat(t)!;
    expect(f.id).toBe("diann");
    expect(f.sampleColumns).toHaveLength(4);
    expect(f.sampleColumns.map(f.clean)).toEqual(["Bait_1", "Bait_2", "Ctrl_1", "Ctrl_2"]);
    expect(f.isFlagged(t.rows[1])).toBe(true);
  });
  it("treats empty cells as missing and keeps the first id", () => {
    const f = detectFormat(t)!;
    const many = tsv(t.columns.join("\t"), ...Array.from({ length: 30 }, (_, i) => `P${i}\tN${i}\tG${i}\tdesc\t${1e6 + i}\t${1.1e6}\t${1e6}\t${1.05e6}`));
    const r = computeVolcano(many, { bait: f.sampleColumns.slice(0, 2), control: f.sampleColumns.slice(2), welch: false, shift: 1.8, shrink: 0.3, minFoldChange: 3, curvature: 3, seed: 1 });
    expect(r.format).toContain("DIA-NN");
    expect(r.proteins[0].id).toBe("P0");
  });
});

describe("Spectronaut", () => {
  const t = tsv(
    "PG.ProteinGroups\tPG.Genes\tPG.ProteinNames\t[1] bait_1.raw.PG.Quantity\t[2] bait_2.raw.PG.Quantity\t[3] ctrl_1.raw.PG.Quantity\t[4] ctrl_2.raw.PG.Quantity",
    "P1;P2\tG1\tN1\t1000\tFiltered\t200\t300",
  );
  it("reads PG.Quantity columns and treats Filtered as missing", () => {
    const f = detectFormat(t)!;
    expect(f.id).toBe("spectronaut");
    expect(f.sampleColumns).toHaveLength(4);
    expect(f.sampleColumns.map(f.clean)).toEqual(["bait_1", "bait_2", "ctrl_1", "ctrl_2"]);
  });
  it("explains a long format report", () => {
    expect(explainUnsupported(tsv("R.FileName\tPG.ProteinGroups\tPG.Quantity", "a\tb\t1"))).toContain("long format");
  });
});

describe("Proteome Discoverer", () => {
  const t = tsv(
    "Accession\tDescription\tGene Symbol\tContaminant\tAbundances (Normalized): F1: Bait 1\tAbundances (Normalized): F2: Bait 2\tAbundances (Normalized): F3: Ctrl 1\tAbundances (Normalized): F4: Ctrl 2\tAbundance: F1: Bait 1",
    "P1\td\tG1\tFalse\t100\t110\t10\t12\t9",
    "P2\td\tG2\tTrue\t100\t110\t10\t12\t9",
  );
  it("prefers normalized abundances and flags contaminants", () => {
    const f = detectFormat(t)!;
    expect(f.id).toBe("proteome-discoverer");
    expect(f.sampleColumns).toHaveLength(4);
    expect(f.sampleColumns.map(f.clean)).toEqual(["Bait 1", "Bait 2", "Ctrl 1", "Ctrl 2"]);
    expect(f.isFlagged(t.rows[1])).toBe(true);
    expect(f.isFlagged(t.rows[0])).toBe(false);
  });
});

describe("Perseus matrix", () => {
  const t = tsv(
    "Protein IDs\tGene names\tbait_1\tbait_2\tctrl_1\tctrl_2\tNote",
    "#!{Type}T\tT\tE\tE\tE\tE\tT",
    "P1\tG1\t25.1\tNaN\t20.2\t20.0\tx",
  );
  it("drops annotation rows and uses Main columns", () => {
    expect(t.rows).toHaveLength(1);
    const f = detectFormat(t)!;
    expect(f.id).toBe("perseus");
    expect(f.sampleColumns).toEqual(["bait_1", "bait_2", "ctrl_1", "ctrl_2"]);
  });
});

describe("generic table and decimal commas", () => {
  const rows = Array.from({ length: 25 }, (_, i) => `P${i}\t${(10 + i / 10).toFixed(2).replace(".", ",")}\t10,5\t9,8\t10,1`);
  const t = tsv("Protein\ta_1\ta_2\tb_1\tb_2", ...rows);
  it("falls back to numeric columns and reads decimal commas", () => {
    const f = detectFormat(t)!;
    expect(f.id).toBe("generic");
    expect(f.sampleColumns).toHaveLength(4);
    const r = computeVolcano(t, { bait: ["a_1", "a_2"], control: ["b_1", "b_2"], welch: false, shift: 1.8, shrink: 0.3, minFoldChange: 3, curvature: 3, seed: 1 });
    expect(r.proteins.length).toBeGreaterThan(0);
    expect(Number.isFinite(r.proteins[0].x)).toBe(true);
  });
  it("returns null when nothing looks like samples", () => {
    expect(detectFormat(tsv("a\tb", "x\ty"))).toBeNull();
  });
});


describe("parser anchoring", () => {
  it("keeps DIA-NN runs whose path contains Sequences or starts with N.", () => {
    const t = tsv("Protein.Group\tGenes\tD:\\Sequences\\run_1.raw\tD:\\Sequences\\run_2.raw\tN.run3.raw\tN.run4.raw\tN.All.Sequences", ...Array.from({ length: 30 }, (_, i) => `P${i}\tG${i}\t1e6\t1e6\t1e6\t1e6\t12`));
    const f = detectFormat(t)!;
    expect(f.sampleColumns).toHaveLength(4);
    expect(f.sampleColumns).not.toContain("N.All.Sequences");
  });
  it("keeps a FragPipe sample whose name starts with Combined or contains Total", () => {
    const t = tsv("Protein\tProtein ID\tGene\tCombined_1 MaxLFQ Intensity\tCombined_2 MaxLFQ Intensity\tTotal_lysate MaxLFQ Intensity\tbait MaxLFQ Intensity\tCombined Total Peptides", "p\tP1\tG\t1\t2\t3\t4\t5");
    expect(detectFormat(t)!.sampleColumns).toHaveLength(4);
  });
  it("cleans raw Proteome Discoverer abundance names", () => {
    const t = tsv("Accession\tGene Symbol\tAbundance: F1: Bait 1\tAbundance: F2: Bait 2\tAbundance: F3: Ctrl 1\tAbundance: F4: Ctrl 2", "P1\tG\t1\t2\t3\t4");
    const f = detectFormat(t)!;
    expect(f.sampleColumns.map(f.clean)).toEqual(["Bait 1", "Bait 2", "Ctrl 1", "Ctrl 2"]);
  });
  it("treats n.d. and N/A as missing when deciding if a column is numeric", () => {
    const t = tsv("Protein\ta_1\ta_2\tb_1\tb_2", ...Array.from({ length: 25 }, (_, i) => `P${i}\t${i % 5 === 0 ? "n.d." : 10 + i}\t11\t${i % 7 === 0 ? "N/A" : 12}\t13`));
    expect(detectFormat(t)!.sampleColumns).toHaveLength(4);
  });
  it("suggests bait before control when both words appear", () => {
    expect(suggestRole("GFP_bait")).toBe("bait");
    expect(suggestRole("GFP_ctrl")).toBe("control");
  });
});
