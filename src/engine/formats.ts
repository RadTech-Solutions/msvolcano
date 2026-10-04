import { isMissingToken, parseNumber, usesDecimalComma } from "./numbers";
import type { Table } from "./parse";

export type FormatId = "maxquant" | "fragpipe" | "diann" | "spectronaut" | "proteome-discoverer" | "perseus" | "generic";

export interface FormatInfo {
  id: FormatId;
  /** Human readable name shown in the interface. */
  label: string;
  /** One or two words for tight spaces. */
  short: string;
  /** Columns that hold one quantity per sample. */
  sampleColumns: string[];
  geneCol: number;
  idCol: number;
  /** True for reverse (decoy) or contaminant rows. */
  isFlagged: (row: string[]) => boolean;
  /** Sample name without software specific decoration, used for grouping replicates. */
  clean: (column: string) => string;
  /** Shown to the user after detection. */
  note?: string;
}

const find = (t: Table, re: RegExp) => t.columns.findIndex((c) => re.test(c));
const CONTAMINANT_ID = /^(contam_|cont_|con__|rev_|rev__)/i;
const stripPath = (c: string) => c.split(/[\\/]/).pop() ?? c;
const stripExt = (c: string) => c.replace(/\.(raw|d|mzml|mzxml|wiff2?|dia|tdf|htrms)$/i, "");

function numericColumn(t: Table, i: number): boolean {
  const sample = t.rows.slice(0, 200);
  const cells = sample.map((r) => r[i]);
  const decimalComma = usesDecimalComma(cells);
  const vals = cells.filter((v) => !isMissingToken(v));
  return vals.length >= Math.max(1, sample.length * 0.2) && vals.every((v) => Number.isFinite(parseNumber(v, decimalComma)));
}

function maxquant(t: Table): FormatInfo | null {
  let sampleColumns = t.columns.filter((c) => /^LFQ intensity /i.test(c));
  if (!sampleColumns.length && find(t, /^majority.protein.ids?$/i) >= 0) sampleColumns = t.columns.filter((c) => /lfq/i.test(c) && !/maxlfq/i.test(c));
  if (!sampleColumns.length) return null;
  const flagCols = [find(t, /contaminant/i), find(t, /^reverse/i)].filter((i) => i >= 0);
  const majority = find(t, /^majority.protein.ids?$/i);
  return {
    id: "maxquant",
    label: "MaxQuant proteinGroups (LFQ)", short: "MaxQuant",
    sampleColumns,
    geneCol: find(t, /^gene.?names?$/i),
    idCol: majority >= 0 ? majority : find(t, /^protein.ids?$/i),
    isFlagged: (r) => flagCols.some((i) => r[i] === "+"),
    clean: (c) => c.replace(/^LFQ intensity\s*/i, ""),
  };
}

function fragpipe(t: Table): FormatInfo | null {
  if (find(t, /^Protein ID$/) < 0 || find(t, /^Gene$/) < 0) return null;
  let sampleColumns = t.columns.filter((c) => / MaxLFQ Intensity$/.test(c));
  let suffix = " MaxLFQ Intensity";
  let note = "Using the MaxLFQ intensity columns.";
  if (!sampleColumns.length) {
    sampleColumns = t.columns.filter((c) => / Intensity$/.test(c) && !/ (MaxLFQ |MaxLFQ Unique |MaxLFQ Total |Unique |Total |Razor )Intensity$/.test(c));
    suffix = " Intensity";
    note = "No MaxLFQ columns found, so the plain Intensity columns are used. Enable MaxLFQ in FragPipe for better quantification.";
  }
  if (!sampleColumns.length) return null;
  const proteinCol = find(t, /^Protein$/);
  const idCol = find(t, /^Protein ID$/);
  return {
    id: "fragpipe", label: "FragPipe combined_protein.tsv", short: "FragPipe", sampleColumns,
    geneCol: find(t, /^Gene$/), idCol,
    isFlagged: (r) => CONTAMINANT_ID.test(r[proteinCol] ?? "") || CONTAMINANT_ID.test(r[idCol] ?? ""),
    clean: (c) => (c.endsWith(suffix) ? c.slice(0, -suffix.length) : c),
    note,
  };
}

// Exact names of the descriptive columns, plus the count columns newer versions add. Anchored so that a run
// whose file path happens to contain "Sequences" or start with "N." is never mistaken for metadata.
const DIANN_META = /^(Protein\.Group|Protein\.Ids|Protein\.Names|Genes|First\.Protein\.Description|N\.(All|Proteotypic)\.Sequences?|Number of .*|N\.Sequences?)$/i;

function diann(t: Table): FormatInfo | null {
  if (find(t, /^Protein\.Group$/) < 0 || (find(t, /^Genes$/) < 0 && find(t, /^Protein\.Names$/) < 0)) return null;
  const sampleColumns = t.columns.filter((c, i) => !DIANN_META.test(c) && numericColumn(t, i));
  if (!sampleColumns.length) return null;
  const idCol = find(t, /^Protein\.Group$/);
  return {
    id: "diann", label: "DIA-NN report.pg_matrix.tsv (beta)", short: "DIA-NN (beta)", sampleColumns,
    geneCol: find(t, /^Genes$/), idCol,
    isFlagged: (r) => CONTAMINANT_ID.test(r[idCol] ?? ""),
    clean: (c) => stripExt(stripPath(c)),
    note: "Written from the DIA-NN documentation and not yet tested on many real files. Check the groups in step 2.",
  };
}

function spectronaut(t: Table): FormatInfo | null {
  const idCol = find(t, /^PG\.ProteinGroups$/);
  if (idCol < 0) return null;
  const sampleColumns = t.columns.filter((c) => /PG\.Quantity$/.test(c));
  if (!sampleColumns.length) return null;
  return {
    id: "spectronaut", label: "Spectronaut protein group report (beta)", short: "Spectronaut (beta)", sampleColumns,
    geneCol: find(t, /^PG\.Genes$/), idCol,
    isFlagged: (r) => CONTAMINANT_ID.test(r[idCol] ?? ""),
    clean: (c) => stripExt(c.replace(/^\[\d+\]\s*/, "").replace(/\.?PG\.Quantity$/, "")),
    note: "Written from the Spectronaut documentation and not yet tested on many real files. Check the groups in step 2.",
  };
}

function proteomeDiscoverer(t: Table): FormatInfo | null {
  const idCol = find(t, /^Accession$/);
  if (idCol < 0) return null;
  let sampleColumns = t.columns.filter((c) => /^Abundances \(Normalized\):/i.test(c));
  let which = "normalized abundances";
  if (!sampleColumns.length) { sampleColumns = t.columns.filter((c) => /^Abundance:/i.test(c)); which = "raw abundances"; }
  if (!sampleColumns.length) return null;
  const contCol = find(t, /^Contaminant$/i);
  return {
    id: "proteome-discoverer", label: "Proteome Discoverer protein export (beta)", short: "Proteome Discoverer (beta)", sampleColumns,
    geneCol: find(t, /^Gene( Symbol)?$/i), idCol,
    isFlagged: (r) => contCol >= 0 && /^(true|\+|1|yes)$/i.test((r[contCol] ?? "").trim()),
    clean: (c) => c.replace(/^Abundances?\s*(\([^)]*\))?:\s*/i, "").replace(/^F\d+:\s*/i, ""),
    note: `Using the ${which}. Written from documentation and not yet tested on many real files. Check the groups in step 2.`,
  };
}

function perseus(t: Table): FormatInfo | null {
  if (!t.types) return null;
  const sampleColumns = t.columns.filter((_, i) => t.types![i] === "E");
  if (sampleColumns.length < 2) return null;
  return {
    id: "perseus", label: "Perseus matrix", short: "Perseus", sampleColumns,
    geneCol: find(t, /^gene/i),
    idCol: find(t, /^(majority protein ids?|protein ids?|accession|uniprot)/i),
    isFlagged: () => false,
    clean: (c) => c.replace(/^LFQ intensity\s*/i, ""),
    note: "The columns marked as Main in Perseus are used as samples. Check that only intensity columns are marked Main, and that missing values were not already filled in inside Perseus.",
  };
}

/** Last resort: treat mostly numeric columns as samples. The user checks the groups in step 2. */
function generic(t: Table): FormatInfo | null {
  if (!t.rows.length) return null;
  const sampleColumns = t.columns.filter((_, i) => numericColumn(t, i));
  if (sampleColumns.length < 4) return null;
  const firstText = t.columns.findIndex((_, i) => !numericColumn(t, i));
  const id = find(t, /^(protein|accession|uniprot|id)/i);
  return {
    id: "generic",
    label: "Generic table (numeric columns as samples)", short: "Plain table",
    sampleColumns,
    geneCol: find(t, /^gene/i),
    idCol: id >= 0 ? id : firstText,
    isFlagged: () => false,
    clean: (c) => stripExt(stripPath(c)),
    note: "Format not recognised. Numeric columns are treated as sample intensities, so check the groups in step 2.",
  };
}

const DETECTORS: ((t: Table) => FormatInfo | null)[] = [perseus, fragpipe, diann, spectronaut, proteomeDiscoverer, maxquant, generic];

/** Detect the software that produced a protein table. Returns null when no sample columns can be found. */
export function detectFormat(t: Table): FormatInfo | null {
  if (t.forceFormat === "generic") return generic(t);
  for (const d of DETECTORS) {
    const f = d(t);
    if (f) return f;
  }
  return null;
}

/** A short reason when a table looks like a known export that this tool cannot read. */
export function explainUnsupported(t: Table): string {
  const cols = t.columns.join("\n");
  if (/^R\.FileName$/m.test(cols) && /^PG\./m.test(cols)) return "This looks like a Spectronaut long format report. Export the protein group pivot table (one column per run) instead.";
  if (/^Intensity\b/m.test(cols) && !/lfq/i.test(cols)) return "This looks like a MaxQuant file without LFQ columns. Run MaxQuant again with LFQ switched on.";
  return "";
}
