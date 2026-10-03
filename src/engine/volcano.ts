import type { Table } from "./parse";
import { mean, randNormal, rng, sd, tTest } from "./stats";

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
}

export interface Protein {
  row: number;
  gene: string;
  id: string;
  x: number;
  y: number;
  /** Mean log2 LFQ of the bait columns after imputation. */
  baitMean: number;
  significant: boolean;
}

export interface VolcanoResult {
  proteins: Protein[];
  /** True when LFQ values looked unlogged and were log2 transformed. */
  logTransformed: boolean;
  dropped: { contaminantOrReverse: number; absentInBait: number; untestable: number };
}

export const defaultParams: Omit<Params, "bait" | "control"> = {
  welch: false,
  shift: 1.8,
  shrink: 0.3,
  minFoldChange: 3,
  curvature: 3,
  seed: 1,
};

const col = (t: Table, re: RegExp) => t.columns.findIndex((c) => re.test(c));

export function lfqColumns(t: Table): string[] {
  return t.columns.filter((c) => /lfq/i.test(c));
}

function toNum(s: string | undefined): number {
  if (s === undefined || s === "" || s === "NaN") return NaN;
  const v = Number(s);
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

export function computeVolcano(table: Table, p: Params): VolcanoResult {
  if (p.bait.length < 2 || p.control.length < 2) {
    throw new Error("Select at least two LFQ columns for bait and for control.");
  }
  const idx = (names: string[]) =>
    names.map((n) => {
      const i = table.columns.indexOf(n);
      if (i < 0) throw new Error(`Column not found: ${n}`);
      return i;
    });
  const bi = idx(p.bait), ci = idx(p.control);
  const geneCol = col(table, /^gene.?names?$/i);
  const idCol = col(table, /^majority.protein.ids?$/i) >= 0 ? col(table, /^majority.protein.ids?$/i) : col(table, /^protein.ids?$/i);
  const flagCols = [col(table, /contaminant/i), col(table, /^reverse/i)].filter((i) => i >= 0);

  const dropped = { contaminantOrReverse: 0, absentInBait: 0, untestable: 0 };
  let rows = table.rows.map((r, row) => ({ row, r }));
  rows = rows.filter(({ r }) => {
    const bad = flagCols.some((i) => r[i] === "+");
    if (bad) dropped.contaminantOrReverse++;
    return !bad;
  });

  let bait = rows.map(({ r }) => bi.map((i) => toNum(r[i])));
  let ctrl = rows.map(({ r }) => ci.map((i) => toNum(r[i])));

  // Drop proteins whose bait values are all at the group minimum (typically all zero).
  const baitMin = Math.min(...bait.flat().filter((v) => !Number.isNaN(v)));
  const keep = bait.map((vals) => {
    const m = mean(vals.filter((v) => !Number.isNaN(v)));
    return !(vals.every((v) => Number.isNaN(v)) || m <= baitMin);
  });
  dropped.absentInBait = keep.filter((k) => !k).length;
  rows = rows.filter((_, i) => keep[i]);
  bait = bait.filter((_, i) => keep[i]);
  ctrl = ctrl.filter((_, i) => keep[i]);

  // Values are treated as already log2 when the mean of the bait column means is below 100.
  const colMeans = bi.map((_, j) => mean(bait.map((v) => v[j]).filter((v) => !Number.isNaN(v))));
  const logTransformed = mean(colMeans) >= 100;
  if (logTransformed) {
    const lg = (m: number[][]) => m.map((v) => v.map((x) => (x > 0 ? Math.log2(x) : NaN)));
    bait = lg(bait);
    ctrl = lg(ctrl);
  }

  const rand = rng(p.seed);
  const impute = (m: number[][]) => {
    const all = m.flat().filter((v) => !Number.isNaN(v));
    const mu = mean(all) - p.shift * sd(all);
    const sigma = p.shrink * sd(all);
    return m.map((v) => v.map((x) => (Number.isNaN(x) ? mu + sigma * randNormal(rand) : x)));
  };
  bait = impute(bait);
  ctrl = impute(ctrl);

  const proteins: Protein[] = [];
  rows.forEach(({ row, r }, i) => {
    const t = tTest(bait[i], ctrl[i], p.welch);
    if (!t) { dropped.untestable++; return; }
    const x = t.diff;
    const y = -Math.log10(t.p);
    const gene = geneCol >= 0 ? (r[geneCol] ?? "").split(";")[0] : "";
    const id = idCol >= 0 ? (r[idCol] ?? "").split(";")[0] : "";
    proteins.push({ row, gene: gene || id, id, x, y, baitMean: mean(bait[i]), significant: isSignificant(x, y, p.minFoldChange, p.curvature) });
  });
  return { proteins, logTransformed, dropped };
}

export function toCsv(table: Table, proteins: Protein[], p: Params): string {
  const header = ["gene", "protein_id", "log2_difference", "neg_log10_p", ...p.bait, ...p.control];
  const bi = p.bait.map((n) => table.columns.indexOf(n));
  const ci = p.control.map((n) => table.columns.indexOf(n));
  const q = (s: string) => (/[",\n\t]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s);
  const lines = [header.join(",")];
  for (const pr of proteins) {
    const r = table.rows[pr.row];
    lines.push([pr.gene, pr.id, pr.x.toFixed(4), pr.y.toFixed(4), ...bi.map((i) => r[i]), ...ci.map((i) => r[i])].map(q).join(","));
  }
  return lines.join("\n");
}
