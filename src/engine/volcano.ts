import { detectFormat } from "./formats";
import type { Table } from "./parse";
import { parseNumber, usesDecimalComma } from "./numbers";
import { mean, quantile, randNormal, rng, sd, tTest } from "./stats";

export type ImputationMethod = "normal" | "mindet" | "minprob" | "none";
/** "separate" keeps proteins seen in bait but never in a control out of the t-test; "include" imputes them as the 2016 tool did. */
export type PresenceMode = "separate" | "include";
/** "auto" decides from the intensities. "raw" means not logged yet. "log2" means already log2. */
export type LogMode = "auto" | "raw" | "log2";

export interface Params {
  bait: string[];
  control: string[];
  welch: boolean;
  /** Imputation cloud: mean = group mean - shift * sd, sd = shrink * group sd. */
  shift: number;
  shrink: number;
  minFoldChange: number;
  curvature: number;
  seed: number;
  /** "normal" is the method from the 2016 paper. The others are new in version 2. */
  imputation?: ImputationMethod;
  /** Defaults to "include" in the engine (the published behaviour). The interface passes "separate". */
  presence?: PresenceMode;
  /** Minimum number of bait replicates with a value for a protein to count as present only in bait. */
  presenceMin?: number;
  /** In "separate" mode a protein needs this many measured values in at least one group to be tested. Default 2. */
  minValid?: number;
  logMode?: LogMode;
}

export interface Protein {
  row: number;
  gene: string;
  id: string;
  x: number;
  y: number;
  /** Mean log2 intensity of the bait columns after imputation. */
  baitMean: number;
  /** Standard error of the difference, used by the experimental s0 statistic. */
  se: number;
  /** Log2 values used for the test, bait columns first then control columns. Imputed where a value was missing. */
  values: number[];
  /** The measured log2 values with NaN where nothing was quantified. Same order as values. */
  raw: number[];
  /** True where the value in `values` was filled in. */
  imputed: boolean[];
  significant: boolean;
}

/** A protein quantified in the bait samples and in none of the controls. Not given a p value. */
export interface PresenceOnly {
  row: number;
  gene: string;
  id: string;
  /** Replicates with a value, out of nBait. */
  nObs: number;
  nBait: number;
  /** Mean log2 intensity over the bait replicates that have a value. */
  meanLog2: number;
}

export interface VolcanoResult {
  proteins: Protein[];
  presenceOnly: PresenceOnly[];
  /** True when intensities looked unlogged and were log2 transformed. */
  logTransformed: boolean;
  dropped: { contaminantOrReverse: number; absentInBait: number; untestable: number; presenceOnly: number; tooFewValues: number };
  /** How the intensities were read, with the evidence, so the interface can show it. */
  scale: { mode: "raw" | "log2"; decided: "auto" | "user"; median: number; max: number; warning: string; decimalComma: boolean };
  /** Name of the detected input format. */
  format: string;
  /** Sample columns in the order bait then control, with how many of the tested proteins had a value in each. */
  coverage: { column: string; observed: number; total: number; group: "bait" | "control" }[];
}

export const defaultParams: Omit<Params, "bait" | "control"> = {
  welch: false,
  shift: 1.8,
  shrink: 0.3,
  minFoldChange: 3,
  curvature: 3,
  seed: 1,
};

/** Sample columns of a protein table, whichever software wrote it. */
export function lfqColumns(t: Table): string[] {
  return detectFormat(t)?.sampleColumns ?? [];
}

/** Hyperbolic cutoff: significant when x > minFoldChange and y > curvature / (x - minFoldChange). */
export function isSignificant(x: number, y: number, minFoldChange: number, curvature: number): boolean {
  return x > minFoldChange && curvature / (x - minFoldChange) < y;
}

export function cutoffCurve(minFoldChange: number, curvature: number, xMax: number, yMax: number, step = 0.05) {
  const x0 = (curvature + yMax * minFoldChange) / yMax;
  const right: { x: number; y: number }[] = [];
  for (let x = x0; x <= xMax; x += step) right.push({ x, y: curvature / (x - minFoldChange) });
  return {
    right,
    left: right.map((p) => ({ x: -p.x, y: p.y })),
  };
}

const observed = (v: number[]) => v.filter((x) => !Number.isNaN(x));

function median(v: number[]): number {
  if (!v.length) return NaN;
  const a = [...v].sort((x, y) => x - y);
  const m = a.length >> 1;
  return a.length % 2 ? a[m] : (a[m - 1] + a[m]) / 2;
}

/** A small hash so each filled-in cell gets its own random draw, independent of which other proteins are present. */
function mix(a: number, b: number, c: number): number {
  let h = (Math.imul(a | 0, 374761393) + Math.imul(b | 0, 668265263) + Math.imul(c | 0, 2147483629)) | 0;
  h = Math.imul(h ^ (h >>> 13), 1274126177);
  h ^= h >>> 16;
  return h >>> 0;
}

export function computeVolcano(table: Table, p: Params): VolcanoResult {
  if (p.bait.length < 2 || p.control.length < 2) {
    throw new Error("Select at least two sample columns for bait and for control.");
  }
  const fmt = detectFormat(table);
  if (!fmt) throw new Error("No sample columns found in this table.");
  const idx = (names: string[]) =>
    names.map((n) => {
      const i = table.columns.indexOf(n);
      if (i < 0) throw new Error(`Column not found: ${n}`);
      return i;
    });
  const bi = idx(p.bait), ci = idx(p.control);
  const method = p.imputation ?? "normal";
  const presenceMode = p.presence ?? "include";
  const presenceMin = p.presenceMin ?? 2;
  const minValid = p.minValid ?? 2;

  const dropped = { contaminantOrReverse: 0, absentInBait: 0, untestable: 0, presenceOnly: 0, tooFewValues: 0 };
  let rows = table.rows.map((r, row) => ({ row, r }));
  rows = rows.filter(({ r }) => {
    const bad = fmt.isFlagged(r);
    if (bad) dropped.contaminantOrReverse++;
    return !bad;
  });

  // One number style per file: a decimal comma in one sample column means all of them.
  const sampleCells = rows.slice(0, 500).flatMap(({ r }) => [...bi, ...ci].map((i) => r[i]));
  const decimalComma = usesDecimalComma(sampleCells);
  const num = (s: string | undefined) => parseNumber(s, decimalComma);

  let bait = rows.map(({ r }) => bi.map((i) => num(r[i])));
  let ctrl = rows.map(({ r }) => ci.map((i) => num(r[i])));

  // Are the intensities already log2? Decided from the typical positive value over every selected column.
  const positives = [...bait.flat(), ...ctrl.flat()].filter((x) => x > 0);
  const med = median(positives), max = positives.length ? Math.max(...positives) : NaN;
  const auto = p.logMode ?? "auto";
  const unlogged = auto === "raw" ? true : auto === "log2" ? false : med >= 100;
  let warning = "";
  if (!unlogged && med > 100) warning = "These values are very large for log2 intensities (median " + med.toExponential(1) + "). They look like raw intensities. Choose the intensity scale in step 3 if this is wrong.";
  if (!unlogged && max < 12) warning = "These values are small for log2 intensities (the largest is " + max.toFixed(1) + "). They may be log10 or natural log, which would make every difference and the cutoff meaningless. Choose the intensity scale in step 3 if this is wrong.";
  if (unlogged && med < 1000) warning = "These look like raw values but are small for intensities (median " + med.toFixed(0) + "). Spectral counts or ratios would also look like this. Choose the intensity scale in step 3 if this is wrong.";
  if (unlogged) {
    // Zero means "not quantified" in the software output, so it becomes missing.
    const lg = (m: number[][]) => m.map((v) => v.map((x) => (x > 0 ? Math.log2(x) : NaN)));
    bait = lg(bait);
    ctrl = lg(ctrl);
  }

  // A selected column with no measured values at all (after zeros became missing) is an input problem, not something to compute around.
  const empty = [...p.bait, ...p.control].filter((_, j) => !rows.some((_, k) => !Number.isNaN((j < p.bait.length ? bait[k][j] : ctrl[k][j - p.bait.length]))));
  if (empty.length) throw new Error(`The column "${empty[0]}" has no values. Untick it in step 2 or check the file.`);

  // Proteins with no value in any bait replicate cannot be enriched in bait.
  const keep = bait.map((v) => observed(v).length > 0);
  dropped.absentInBait = keep.filter((k) => !k).length;
  rows = rows.filter((_, i) => keep[i]);
  bait = bait.filter((_, i) => keep[i]);
  ctrl = ctrl.filter((_, i) => keep[i]);

  const rawBait = bait, rawCtrl = ctrl;
  const nb = p.bait.length;
  // Every filled-in cell has its own random draw, keyed by its row and column, so changing a filter
  // or a setting elsewhere never reshuffles the values of other proteins.
  const cellNormal = (rowId: number, col: number) => randNormal(rng(mix(p.seed, rowId, col)));
  const imputeLow = (m: number[][], probabilistic: boolean, offset: number) => {
    // Left censored imputation from a low quantile of each column (Lazar et al. 2016).
    const nCols = m[0]?.length ?? 0;
    const q = Array.from({ length: nCols }, (_, j) => quantile(observed(m.map((v) => v[j])), 0.01));
    const rowSds = m.map((v) => observed(v)).filter((v) => v.length > 1).map(sd);
    const spread = probabilistic ? quantile(rowSds, 0.5) : 0;
    return m.map((v, i) => v.map((x, j) => (Number.isNaN(x) ? q[j] + spread * (probabilistic ? cellNormal(rows[i].row, offset + j) : 0) : x)));
  };
  const impute = (m: number[][], offset: number) => {
    if (method === "none") return m;
    if (method === "mindet") return imputeLow(m, false, offset);
    if (method === "minprob") return imputeLow(m, true, offset);
    // As published, each group gets its own cloud: the group mean minus shift standard deviations.
    const all = observed(m.flat());
    const mu = mean(all) - p.shift * sd(all);
    const sigma = p.shrink * sd(all);
    return m.map((v, i) => v.map((x, j) => (Number.isNaN(x) ? mu + sigma * cellNormal(rows[i].row, offset + j) : x)));
  };
  // Imputation statistics come from every kept protein, including the ones set aside below.
  const filledBait = impute(rawBait, 0);
  const filledCtrl = impute(rawCtrl, nb);

  const proteins: Protein[] = [];
  const presenceOnly: PresenceOnly[] = [];
  const geneOf = (r: string[]) => (fmt.geneCol >= 0 ? (r[fmt.geneCol] ?? "").split(";")[0] : "");
  const idOf = (r: string[]) => (fmt.idCol >= 0 ? (r[fmt.idCol] ?? "").split(";")[0] : "");

  rows.forEach(({ row, r }, i) => {
    const nObsB = observed(rawBait[i]).length, nObsC = observed(rawCtrl[i]).length;
    const gene = geneOf(r), id = idOf(r);
    if (presenceMode === "separate" && nObsC === 0 && nObsB >= presenceMin) {
      dropped.presenceOnly++;
      presenceOnly.push({ row, gene: gene || id, id, nObs: nObsB, nBait: rawBait[i].length, meanLog2: mean(observed(rawBait[i])) });
      return;
    }
    // A protein measured once in bait and never in a control would be tested on filled-in values alone.
    if (presenceMode === "separate" && nObsB < minValid && nObsC < minValid) { dropped.tooFewValues++; return; }
    // Without imputation the test uses only the measured values and needs two per group.
    const a = method === "none" ? observed(rawBait[i]) : filledBait[i];
    const b = method === "none" ? observed(rawCtrl[i]) : filledCtrl[i];
    const t = tTest(a, b, p.welch);
    if (!t) { dropped.untestable++; return; }
    const x = t.diff;
    const y = -Math.log10(t.p);
    const values = [...filledBait[i], ...filledCtrl[i]];
    const raw = [...rawBait[i], ...rawCtrl[i]];
    proteins.push({
      row, gene: gene || id, id, x, y, baitMean: mean(observed(filledBait[i])), se: t.se,
      values, raw, imputed: raw.map((v) => Number.isNaN(v)), significant: isSignificant(x, y, p.minFoldChange, p.curvature),
    });
  });
  presenceOnly.sort((a, b) => b.nObs - a.nObs || b.meanLog2 - a.meanLog2);

  const columns = [...p.bait, ...p.control];
  const coverage = columns.map((column, j) => ({
    column, group: (j < p.bait.length ? "bait" : "control") as "bait" | "control",
    observed: proteins.filter((q) => !Number.isNaN(q.raw[j])).length, total: proteins.length,
  }));
  return {
    proteins, presenceOnly, logTransformed: unlogged, dropped, format: fmt.label, coverage,
    scale: { mode: unlogged ? "raw" : "log2", decided: auto === "auto" ? "auto" : "user", median: med, max, warning, decimalComma },
  };
}

export function toCsv(table: Table, proteins: Protein[], p: Params, stoich?: Map<number, number>): string {
  const header = ["gene", "protein_id", "log2_difference", "neg_log10_p", ...(stoich?.size ? ["stoichiometry"] : []), ...p.bait, ...p.control];
  const bi = p.bait.map((n) => table.columns.indexOf(n));
  const ci = p.control.map((n) => table.columns.indexOf(n));
  const q = (s: string) => (/[",\n\t]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s);
  const lines = [header.join(",")];
  for (const pr of proteins) {
    const r = table.rows[pr.row];
    lines.push([pr.gene, pr.id, pr.x.toFixed(4), pr.y.toFixed(4), ...(stoich?.size ? [String(stoich.get(pr.row) ?? "")] : []), ...bi.map((i) => r[i]), ...ci.map((i) => r[i])].map(q).join(","));
  }
  return lines.join("\n");
}

/** CSV of the proteins seen in bait but in none of the controls. */
export function presenceCsv(table: Table, rows: PresenceOnly[], p: Params): string {
  const header = ["gene", "protein_id", "bait_replicates_with_value", "bait_replicates", "mean_log2_in_bait", ...p.bait, ...p.control];
  const bi = p.bait.map((n) => table.columns.indexOf(n));
  const ci = p.control.map((n) => table.columns.indexOf(n));
  const q = (s: string) => (/[",\n\t]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s);
  const lines = [header.join(",")];
  for (const pr of rows) {
    const r = table.rows[pr.row];
    lines.push([pr.gene, pr.id, String(pr.nObs), String(pr.nBait), pr.meanLog2.toFixed(3), ...bi.map((i) => r[i]), ...ci.map((i) => r[i])].map(q).join(","));
  }
  return lines.join("\n");
}

export interface SensitivityRow { label: string; hits: number; shared: number; lost: number; gained: number; base: boolean }

/**
 * Re-run the analysis with other fill-in shifts and other random draws, and compare the hyperbolic-curve
 * hit list with the current one. Shows how much of the result depends on the arbitrary settings of the
 * fill-in. Returns null when nothing is filled in.
 */
export function shiftSensitivity(table: Table, p: Params, shifts: number[], seeds: number[] = []): SensitivityRow[] | null {
  if ((p.imputation ?? "normal") === "none") return null;
  const hitSet = (q: Params) => {
    const r = computeVolcano(table, q);
    return new Set(r.proteins.filter((x) => isSignificant(x.x, x.y, p.minFoldChange, p.curvature)).map((x) => x.row));
  };
  const base = hitSet(p);
  const row = (label: string, set: Set<number>, isBase: boolean): SensitivityRow => {
    const shared = [...set].filter((r) => base.has(r)).length;
    return { label, hits: set.size, shared, lost: base.size - shared, gained: set.size - shared, base: isBase };
  };
  const out: SensitivityRow[] = [];
  for (const shift of shifts) out.push(shift === p.shift ? row(`shift ${shift} (current)`, base, true) : row(`shift ${shift}`, hitSet({ ...p, shift }), false));
  for (const seed of seeds) out.push(row(`other random draw ${seed}`, hitSet({ ...p, seed }), false));
  return out;
}
