import { describe, expect, it } from "vitest";
import { benjaminiHochberg, permutationFdr, replicateCorrelation } from "../src/engine/experimental";
import { quantile, randNormal, rng } from "../src/engine/stats";
import type { Protein } from "../src/engine/volcano";

const prot = (values: number[], row = 0): Protein => ({ row, gene: "g" + row, id: "i" + row, x: 0, y: 0, baitMean: 0, se: 0, values, raw: values, imputed: values.map(() => false), significant: false });

describe("benjaminiHochberg", () => {
  it("matches a hand computed example", () => {
    // Hand computed: sorted p 0.005, 0.01, 0.03, 0.04 give 0.02, 0.02, 0.04, 0.04
    const q = benjaminiHochberg([0.01, 0.04, 0.03, 0.005]);
    [0.02, 0.04, 0.04, 0.02].forEach((v, i) => expect(q[i]).toBeCloseTo(v, 12));
  });
  it("never exceeds 1 and keeps order", () => {
    const q = benjaminiHochberg([0.9, 0.5, 0.99]);
    expect(Math.max(...q)).toBeLessThanOrEqual(1);
  });
});

describe("quantile", () => {
  it("matches R type 7", () => {
    expect(quantile([1, 2, 3, 4], 0.5)).toBe(2.5);
    expect(quantile([1, 2, 3, 4], 0.01)).toBeCloseTo(1.03, 12);
    expect(quantile([], 0.5)).toBeNaN();
  });
});

describe("permutationFdr", () => {
  // 3 bait vs 3 control; ten strongly enriched proteins and 90 flat ones with noise
  let seed = 3;
  const noise = () => { seed = (seed * 16807) % 2147483647; return (seed / 2147483647 - 0.5) * 0.6; };
  const proteins: Protein[] = [];
  for (let i = 0; i < 100; i++) {
    const base = 20 + noise();
    const shift = i < 10 ? 6 : 0;
    proteins.push(prot([base + shift + noise(), base + shift + noise(), base + shift + noise(), base + noise(), base + noise(), base + noise()], i));
  }
  it("calls the enriched proteins and few others", () => {
    const r = permutationFdr(proteins, 3, 0.1, 0.05, 19, 1)!;
    expect(r.permutations).toBeGreaterThan(5);
    for (let i = 0; i < 10; i++) expect(r.significant.has(i)).toBe(true);
    expect(r.significant.size).toBeLessThanOrEqual(14);
  });
  it("is reproducible", () => {
    const a = permutationFdr(proteins, 3, 0.1, 0.05, 19, 5)!, b = permutationFdr(proteins, 3, 0.1, 0.05, 19, 5)!;
    expect([...a.significant]).toEqual([...b.significant]);
  });
});

describe("permutationFdr calibration", () => {
  it("makes almost no calls on pure noise", () => {
    const r = rng(11);
    const noise: Protein[] = Array.from({ length: 500 }, (_, i) => prot(Array.from({ length: 8 }, () => 20 + randNormal(r)), i));
    const res = permutationFdr(noise, 4, 0.1, 0.05, 250, 1)!;
    expect(res.significant.size).toBeLessThanOrEqual(5);
  });
});

describe("replicateCorrelation", () => {
  it("is 1 on the diagonal and 1 for identical columns", () => {
    const ps = [prot([1, 1, 5]), prot([2, 2, 4]), prot([3, 3, 9])];
    const c = replicateCorrelation(ps);
    expect(c[0][0]).toBe(1);
    expect(c[0][1]).toBeCloseTo(1, 12);
  });
});
