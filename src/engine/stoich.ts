import type { Protein } from "./volcano";

export interface PeptideInfo {
  gene: string;
  tryptic: number;
}

/** Fallback count for proteins missing from the table (as in the 2016 tool), avoids division by zero. */
export const MISSING_PEPTIDE_COUNT = 0.01;

const MIN_LEN = 7;
const MAX_LEN = 30;

/**
 * Number of theoretical tryptic peptides of 7 to 30 residues:
 * cleave after K or R unless followed by P, no missed cleavages.
 * Keep identical to the copy in scripts/build-peptides.mjs.
 */
export function countTrypticPeptides(seq: string): number {
  let count = 0;
  let start = 0;
  for (let i = 0; i < seq.length; i++) {
    const c = seq[i];
    const end = i === seq.length - 1;
    if (end || ((c === "K" || c === "R") && seq[i + 1] !== "P")) {
      const len = i + 1 - start;
      if (len >= MIN_LEN && len <= MAX_LEN) count++;
      start = i + 1;
    }
  }
  return count;
}

/** Parse a public/peptides/<taxid>.tsv table (header: accession, gene, tryptic_peptides). */
export function parsePeptideTable(text: string): Map<string, PeptideInfo> {
  const map = new Map<string, PeptideInfo>();
  const lines = text.split(/\r?\n/);
  for (let i = 0; i < lines.length; i++) {
    const line = lines[i];
    if (!line) continue;
    const f = line.split("\t");
    if (i === 0 && f[0] === "accession") continue;
    const tryptic = Number(f[2]);
    if (!f[0] || !Number.isFinite(tryptic)) continue;
    map.set(f[0], { gene: f[1] ?? "", tryptic });
  }
  return map;
}

function lookup(peptides: Map<string, PeptideInfo>, id: string): number {
  const direct = peptides.get(id);
  if (direct) return direct.tryptic;
  const stripped = id.replace(/-\d+$/, "");
  const iso = peptides.get(stripped);
  return iso ? iso.tryptic : MISSING_PEPTIDE_COUNT;
}

function findBait(proteins: Protein[], baitKey: string): Protein | undefined {
  const key = baitKey.trim();
  const lower = key.toLowerCase();
  return (
    proteins.find((p) => p.id === key) ??
    proteins.find((p) => p.id.replace(/-\d+$/, "") === key) ??
    proteins.find((p) => p.gene.toLowerCase() === lower)
  );
}

/** Intensity above control in linear space: 2^baitMean - 2^(baitMean - diff). */
function ip(p: Protein): number {
  return 2 ** p.baitMean - 2 ** (p.baitMean - p.x);
}

/**
 * Stoichiometry relative to the bait (iBAQ like, Singh et al. 2016):
 * sIp = Ip / tryptic peptides, st = sIp / sIp(bait).
 * Proteins absent from the table use a count of 0.01.
 * Proteins with non positive Ip (not enriched over control) get NaN.
 * Throws if the bait is not found or has non positive Ip.
 * The map is keyed by Protein.row.
 */
export function computeStoichiometry(
  proteins: Protein[],
  peptides: Map<string, PeptideInfo>,
  baitKey: string,
): Map<number, number> {
  const bait = findBait(proteins, baitKey);
  if (!bait) throw new Error(`Bait "${baitKey}" not found among the quantified proteins`);
  const baitIp = ip(bait);
  if (!(baitIp > 0)) throw new Error(`Bait "${baitKey}" is not enriched over control (non positive Ip)`);
  const baitS = baitIp / lookup(peptides, bait.id);
  const out = new Map<number, number>();
  for (const p of proteins) {
    const v = ip(p);
    out.set(p.row, v > 0 ? v / lookup(peptides, p.id) / baitS : NaN);
  }
  return out;
}
