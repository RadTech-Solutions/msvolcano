// Experimental features. Not part of the 2016 publication, not peer reviewed or benchmarked.
import { mean, rng } from "./stats";
import type { Protein } from "./volcano";

/** Benjamini and Hochberg (1995) adjusted p values (q values), same order as the input. */
export function benjaminiHochberg(p: number[]): number[] {
  const n = p.length;
  const order = p.map((v, i) => [v, i] as const).sort((a, b) => b[0] - a[0]);
  const q = new Array<number>(n);
  let prev = 1;
  order.forEach(([v, i], k) => {
    prev = Math.min(prev, (v * n) / (n - k));
    q[i] = prev;
  });
  return q;
}

/** Student t statistic style score d = difference / (se + s0), after Tusher et al. 2001 and Tyanova et al. 2016. */
function scores(values: number[], nBait: number, s0: number): number {
  const n1 = nBait, n2 = values.length - nBait;
  let m1 = 0, m2 = 0;
  for (let i = 0; i < n1; i++) m1 += values[i];
  for (let i = n1; i < values.length; i++) m2 += values[i];
  m1 /= n1; m2 /= n2;
  let ss = 0;
  for (let i = 0; i < n1; i++) ss += (values[i] - m1) ** 2;
  for (let i = n1; i < values.length; i++) ss += (values[i] - m2) ** 2;
  const se = Math.sqrt((ss / (n1 + n2 - 2)) * (1 / n1 + 1 / n2));
  return (m1 - m2) / (se + s0);
}

export interface PermutationResult {
  /** Indices (into the input array) of proteins called enriched in bait at the requested FDR. */
  significant: Set<number>;
  threshold: number;
  fdr: number;
  permutations: number;
}

/**
 * Permutation FDR for enrichment in bait. Sample labels are shuffled, the score d is recomputed for every
 * protein, and the threshold is the smallest d where (median false calls) / (observed calls) <= fdr.
 */
export function permutationFdr(proteins: Protein[], nBait: number, s0: number, fdr: number, permutations: number, seed: number): PermutationResult | null {
  if (!proteins.length) return null;
  const nCols = proteins[0].values.length;
  const d = proteins.map((q) => scores(q.values, nBait, s0));
  const rand = rng(seed);
  const perm: number[][] = [];
  const idx = Array.from({ length: nCols }, (_, i) => i);
  const seen = new Set<string>();
  for (let t = 0, guard = 0; t < permutations && guard < permutations * 20; guard++) {
    for (let i = nCols - 1; i > 0; i--) { const j = Math.floor(rand() * (i + 1)); [idx[i], idx[j]] = [idx[j], idx[i]]; }
    const key = [...idx.slice(0, nBait)].sort((a, b) => a - b).join(",");
    if (key === Array.from({ length: nBait }, (_, i) => i).join(",") || seen.has(key)) continue; // skip the true labelling and repeats
    seen.add(key);
    const dp = proteins.map((q) => scores(idx.map((k) => q.values[k]), nBait, s0)).sort((a, b) => a - b);
    perm.push(dp);
    t++;
  }
  if (!perm.length) return null;
  const countGe = (sorted: number[], thr: number) => {
    let lo = 0, hi = sorted.length;
    while (lo < hi) { const mid = (lo + hi) >> 1; if (sorted[mid] >= thr) hi = mid; else lo = mid + 1; }
    return sorted.length - lo;
  };
  const candidates = d.filter((v) => v > 0).sort((a, b) => a - b);
  let threshold = Infinity;
  for (const thr of candidates) {
    const observed = d.filter((v) => v >= thr).length;
    const falseCalls = perm.map((sorted) => countGe(sorted, thr)).sort((a, b) => a - b);
    const medianFalse = falseCalls[Math.floor(falseCalls.length / 2)];
    if (medianFalse / observed <= fdr) { threshold = thr; break; }
  }
  const significant = new Set<number>();
  d.forEach((v, i) => { if (v >= threshold) significant.add(i); });
  return { significant, threshold, fdr, permutations: perm.length };
}

/** Pearson correlation matrix between replicate columns (imputed log2 values). */
export function replicateCorrelation(proteins: Protein[]): number[][] {
  const nCols = proteins[0]?.values.length ?? 0;
  const cols = Array.from({ length: nCols }, (_, j) => proteins.map((q) => q.values[j]));
  const ms = cols.map(mean);
  const cor = (a: number, b: number) => {
    let sab = 0, saa = 0, sbb = 0;
    for (let i = 0; i < cols[a].length; i++) {
      const x = cols[a][i] - ms[a], y = cols[b][i] - ms[b];
      sab += x * y; saa += x * x; sbb += y * y;
    }
    return sab / Math.sqrt(saa * sbb);
  };
  return cols.map((_, a) => cols.map((__, b) => (a === b ? 1 : cor(a, b))));
}
