import { detectFormat } from "./formats";
import type { Table } from "./parse";
import { mean, quantile, randNormal, rng, sd, tTest } from "./stats";

export type ImputationMethod = "normal" | "mindet" | "minprob" | "none";
/** "separate" keeps proteins seen in bait but never in a control out of the t-test; "include" imputes them as the 2016 tool did. */
export type PresenceMode = "separate" | "include";

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
  dropped: { contaminantOrReverse: number; absentInBait: number; untestable: number; presenceOnly: number };
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

const MISSING = new Set(["", "NaN", "NA", "N/A", "#N/A", "Filtered", "nan", "NULL"]);
function toNum(s: string | undefined): number {
  if (s === undefined || MISSING.has(s.trim())) return NaN;
  // A decimal comma (European export) is accepted when there is no decimal point.
  const v = Number(s.includes(".") ? s : s.replace(",", "."));
  return Number.isFinite(v) ? v : NaN;
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

  const dropped = { contaminantOrReverse: 0, absentInBait: 0, untestable: 0, presenceOnly: 0 };
  let rows = table.rows.map((r, row) => ({ row, r }));
  rows = rows.filter(({ r }) => {
    const bad = fmt.isFlagged(r);
    if (bad) dropped.contaminantOrReverse++;
    return !bad;
  });

  let bait = rows.map(({ r }) => bi.map((i) => toNum(r[i])));
  let ctrl = rows.map(({ r }) => ci.map((i) => toNum(r[i])));

  // Intensities are treated as already log2 when the mean of the bait column means is below 100.
  const colMeans = bi.map((_, j) => mean(observed(bait.map((v) => v[j]))));
  const logTransformed = mean(colMeans) >= 100;
  if (logTransformed) {
    // Zero means "not quantified" in the software output, so it becomes missing.
    const lg = (m: number[][]) => m.map((v) => v.map((x) => (x > 0 ? Math.log2(x) : NaN)));
    bait = lg(bait);
    ctrl = lg(ctrl);
  }

  // Proteins with no value in any bait replicate cannot be enriched in bait.
  const keep = bait.map((v) => observed(v).length > 0);
  dropped.absentInBait = keep.filter((k) => !k).length;
  rows = rows.filter((_, i) => keep[i]);
  bait = bait.filter((_, i) => keep[i]);
  ctrl = ctrl.filter((_, i) => keep[i]);

  const rawBait = bait, rawCtrl = ctrl;
  const rand = rng(p.seed);
  const imputeLow = (m: number[][], probabilistic: boolean) => {
    // Left censored imputation from a low quantile of each column (Lazar et al. 2016).
    const nCols = m[0]?.length ?? 0;
    const q = Array.from({ length: nCols }, (_, j) => quantile(observed(m.map((v) => v[j])), 0.01));
    const rowSds = m.map((v) => observed(v)).filter((v) => v.length > 1).map(sd);
    const spread = probabilistic ? quantile(rowSds, 0.5) : 0;
    return m.map((v) => v.map((x, j) => (Number.isNaN(x) ? q[j] + spread * randNormal(rand) : x)));
  };
  const impute = (m: number[][]) => {
    if (method === "none") return m;
    if (method === "mindet") return imputeLow(m, false);
    if (method === "minprob") return imputeLow(m, true);
    const all = observed(m.flat());
    const mu = mean(all) - p.shift * sd(all);
    const sigma = p.shrink * sd(all);
    return m.map((v) => v.map((x) => (Number.isNaN(x) ? mu + sigma * randNormal(rand) : x)));
  };
  // Imputation statistics come from every kept protein, including the ones set aside below.
  const filledBait = impute(rawBait);
  const filledCtrl = impute(rawCtrl);

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
      row, gene: gene || id, id, x, y, baitMean: mean(observed(filledBait[i])), se: t.diff / t.t,
      values, raw, imputed: raw.map((v) => Number.isNaN(v)), significant: isSignificant(x, y, p.minFoldChange, p.curvature),
    });
  });
  presenceOnly.sort((a, b) => b.nObs - a.nObs || b.meanLog2 - a.meanLog2);

  const columns = [...p.bait, ...p.control];
  const coverage = columns.map((column, j) => ({
    column, group: (j < p.bait.length ? "bait" : "control") as "bait" | "control",
    observed: proteins.filter((q) => !Number.isNaN(q.raw[j])).length, total: proteins.length,
  }));
  return { proteins, presenceOnly, logTransformed, dropped, format: fmt.label, coverage };
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

export interface SensitivityRow { shift: number; hits: number; shared: number; lost: number; gained: number; base: boolean }

/**
 * Re-run the analysis at other imputation shifts and compare the hyperbolic-curve hit list with the current one.
 * Shows how much of the result depends on the arbitrary shift. Returns null when nothing is imputed.
 */
export function shiftSensitivity(table: Table, p: Params, shifts: number[]): SensitivityRow[] | null {
  if ((p.imputation ?? "normal") === "none") return null;
  const hitSet = (shift: number) => {
    const r = computeVolcano(table, { ...p, shift });
    return new Set(r.proteins.filter((q) => isSignificant(q.x, q.y, p.minFoldChange, p.curvature)).map((q) => q.row));
  };
  const base = hitSet(p.shift);
  return shifts.map((shift) => {
    const s = shift === p.shift ? base : hitSet(shift);
    const shared = [...s].filter((r) => base.has(r)).length;
    return { shift, hits: s.size, shared, lost: base.size - shared, gained: s.size - shared, base: shift === p.shift };
  });
}
