import { mean } from "./stats";
import type { PresenceOnly, Protein } from "./volcano";

const observed = (v: number[]) => v.filter((x) => !Number.isNaN(x));

function median(v: number[]): number {
  if (!v.length) return NaN;
  const a = [...v].sort((x, y) => x - y);
  const m = a.length >> 1;
  return a.length % 2 ? a[m] : (a[m - 1] + a[m]) / 2;
}

/** Pearson correlation between sample columns using only proteins measured in both. Imputed values are ignored. */
export function observedCorrelation(proteins: Protein[]): number[][] {
  const n = proteins[0]?.raw.length ?? 0;
  const out: number[][] = Array.from({ length: n }, () => new Array<number>(n).fill(NaN));
  for (let a = 0; a < n; a++) {
    out[a][a] = 1;
    for (let b = a + 1; b < n; b++) {
      const xs: number[] = [], ys: number[] = [];
      for (const q of proteins) if (!Number.isNaN(q.raw[a]) && !Number.isNaN(q.raw[b])) { xs.push(q.raw[a]); ys.push(q.raw[b]); }
      let r = NaN;
      // Too few shared proteins give an unstable correlation, so leave it blank.
      if (xs.length >= 20) {
        const mx = mean(xs), my = mean(ys);
        let sxy = 0, sxx = 0, syy = 0;
        for (let i = 0; i < xs.length; i++) { const dx = xs[i] - mx, dy = ys[i] - my; sxy += dx * dy; sxx += dx * dx; syy += dy * dy; }
        r = sxy / Math.sqrt(sxx * syy);
      }
      out[a][b] = out[b][a] = r;
    }
  }
  return out;
}

/** Eigen decomposition of a small symmetric matrix by Jacobi rotations. Eigenvalues sorted descending. */
export function symmetricEigen(A: number[][]): { values: number[]; vectors: number[][] } {
  const n = A.length;
  const a = A.map((r) => [...r]);
  const v: number[][] = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));
  let norm = 0;
  for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) norm += A[i][j] ** 2;
  for (let sweep = 0; sweep < 100; sweep++) {
    let off = 0;
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) off += a[i][j] ** 2;
    if (off <= 1e-22 * Math.max(norm, 1e-300)) break;
    for (let p = 0; p < n - 1; p++) {
      for (let q = p + 1; q < n; q++) {
        if (Math.abs(a[p][q]) < 1e-300) continue;
        const theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
        const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
        const c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (let k = 0; k < n; k++) { const akp = a[k][p], akq = a[k][q]; a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq; }
        for (let k = 0; k < n; k++) { const apk = a[p][k], aqk = a[q][k]; a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk; }
        for (let k = 0; k < n; k++) { const vkp = v[k][p], vkq = v[k][q]; v[k][p] = c * vkp - s * vkq; v[k][q] = s * vkp + c * vkq; }
      }
    }
  }
  const order = a.map((_, i) => i).sort((i, j) => a[j][j] - a[i][i]);
  return { values: order.map((i) => a[i][i]), vectors: order.map((i) => v.map((row) => row[i])) };
}

export interface PcaResult {
  /** Coordinates of every sample on PC1 and PC2. */
  scores: [number, number][];
  /** Share of variance explained by PC1 and PC2, between 0 and 1. */
  explained: [number, number];
  /** Number of proteins measured in every sample, which is all that PCA uses. */
  proteinsUsed: number;
}

/**
 * PCA of samples on the proteins that were measured in all of them. Imputed values are not used,
 * so the picture is not shaped by the imputation settings. Returns null when too few proteins qualify.
 */
export function samplePca(proteins: Protein[]): PcaResult | null {
  const n = proteins[0]?.raw.length ?? 0;
  const complete = proteins.filter((q) => q.raw.every((v) => !Number.isNaN(v)));
  if (n < 3 || complete.length < 10) return null;
  const X = complete.map((q) => { const m = mean(q.raw); return q.raw.map((v) => v - m); });
  const G = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => { let s = 0; for (const row of X) s += row[i] * row[j]; return s; }));
  const { values, vectors } = symmetricEigen(G);
  const total = values.reduce((s, x) => s + Math.max(0, x), 0) || 1;
  const comp = (m: number) => vectors[m].map((c) => c * Math.sqrt(Math.max(0, values[m])));
  const c1 = comp(0), c2 = vectors.length > 1 ? comp(1) : new Array<number>(n).fill(0);
  return { scores: c1.map((x, i) => [x, c2[i]] as [number, number]), explained: [Math.max(0, values[0]) / total, Math.max(0, values[1] ?? 0) / total], proteinsUsed: complete.length };
}

export type CheckStatus = "ok" | "warn" | "info";
export interface Check { id: string; status: CheckStatus; title: string; detail: string }

/** Median log2 intensity of each sample over its measured proteins. */
export function sampleMedians(proteins: Protein[]): number[] {
  const n = proteins[0]?.raw.length ?? 0;
  return Array.from({ length: n }, (_, j) => median(observed(proteins.map((q) => q.raw[j]))));
}

/** Warn when bait and control samples were loaded or measured at clearly different overall levels. */
export function loadingCheck(proteins: Protein[], nBait: number): Check {
  const med = sampleMedians(proteins);
  const g1 = mean(med.slice(0, nBait)), g2 = mean(med.slice(nBait));
  const d = g1 - g2;
  if (!Number.isFinite(d)) return { id: "loading", status: "info", title: "Loading balance", detail: "Not enough values to compare the overall signal of bait and control." };
  return Math.abs(d) > 1
    ? { id: "loading", status: "warn", title: "Loading balance", detail: `The typical protein is ${Math.abs(d).toFixed(1)} log2 units ${d > 0 ? "higher" : "lower"} in bait than in control (${(2 ** Math.abs(d)).toFixed(1)} fold). Unequal loading or acquisition shifts every point on the volcano.` }
    : { id: "loading", status: "ok", title: "Loading balance", detail: `Bait and control samples have similar overall signal (difference ${d.toFixed(2)} log2).` };
}

/**
 * Warn when one sample has far fewer quantified proteins than the others of its group.
 * Control samples are expected to have fewer values overall (the tested proteins were picked by their bait
 * signal), so a low share alone is not a problem. A sample that lags its own group is.
 */
export function coverageCheck(coverage: { column: string; observed: number; total: number; group?: "bait" | "control" }[]): Check {
  const title = "Quantified proteins per sample";
  if (!coverage.length) return { id: "coverage", status: "info", title, detail: "" };
  const frac = coverage.map((c) => (c.total ? c.observed / c.total : 0));
  const groupMedian = (g: string | undefined) => median(frac.filter((_, k) => coverage[k].group === g));
  let worstGap = 0, wi = -1;
  frac.forEach((f, k) => { const gap = groupMedian(coverage[k].group) - f; if (gap > worstGap) { worstGap = gap; wi = k; } });
  const lowest = Math.min(...frac);
  if (wi >= 0 && worstGap > 0.2 && frac[wi] < 0.7) {
    return { id: "coverage", status: "warn", title, detail: `${coverage[wi].column} has a value for ${(frac[wi] * 100).toFixed(0)}% of the tested proteins, well below the other ${coverage[wi].group ?? ""} samples (typically ${(groupMedian(coverage[wi].group) * 100).toFixed(0)}%). Heavy missingness makes filled-in values dominate that sample.` };
  }
  return { id: "coverage", status: "ok", title, detail: `No sample lags its group. The lowest share of proteins with a value is ${(lowest * 100).toFixed(0)}%.` };
}

export interface BaitCheck extends Check { found: boolean }

const norm = (s: string) => s.trim().toLowerCase();

/** Is the bait the most enriched protein? A bait that is not enriched usually means a failed pull down or the wrong group. */
export function baitRecovery(proteins: Protein[], presence: PresenceOnly[], key: string): BaitCheck {
  const title = "Bait recovery";
  const ranked = [...proteins].sort((a, b) => b.x - a.x);
  const top = ranked[0];
  const k = norm(key);
  if (!k) {
    return { id: "bait", status: "info", found: false, title, detail: top ? `Enter the bait's gene name in step 4 to check that it is enriched. The most enriched protein is ${top.gene}.` : "Enter the bait's gene name in step 4 to check that it is enriched." };
  }
  const match = (g: string, id: string) => norm(g) === k || norm(id) === k || norm(id).split("-")[0] === k;
  const idx = ranked.findIndex((q) => match(q.gene, q.id));
  if (idx >= 0) {
    const q = ranked[idx];
    // A bait is often not the single top hit (its partners can be enriched more), so only a clearly low rank is flagged.
    const limit = Math.max(5, Math.ceil(ranked.length * 0.01));
    return idx < limit
      ? { id: "bait", status: "ok", found: true, title, detail: `${q.gene} is among the most enriched proteins (number ${idx + 1} of ${ranked.length}, difference ${q.x.toFixed(1)} log2).` }
      : { id: "bait", status: "warn", found: true, title, detail: `${q.gene} ranks number ${idx + 1} of ${ranked.length} by enrichment (difference ${q.x.toFixed(1)} log2). The top protein is ${top.gene}. A bait far from the top can mean a failed pull down or mislabelled groups.` };
  }
  const p = presence.find((q) => match(q.gene, q.id));
  if (p) return { id: "bait", status: "ok", found: true, title, detail: `${p.gene} was seen in ${p.nObs} of ${p.nBait} bait replicates and in no control, as expected for the bait.` };
  return { id: "bait", status: "warn", found: false, title, detail: `${key} was not found among the quantified proteins. Check the spelling, or whether it was removed as a contaminant or as absent from bait.` };
}
