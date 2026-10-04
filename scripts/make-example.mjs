// README: every value produced by this script is SIMULATED. Nothing here comes from a real
// experiment. Accessions carry a SYN prefix and gene names are fictional, so the example
// cannot be mistaken for real data. Output is deterministic (seeded RNG).
//
// Usage: node scripts/make-example.mjs
// Writes public/example/proteinGroups_example.txt, a MaxQuant style proteinGroups.txt for an
// imaginary pull-down of a bait called BAITX (4 bait replicates vs 4 control replicates).
import { mkdirSync, writeFileSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

let state = 20240607;
function rand() {
  // mulberry32
  state = (state + 0x6d2b79f5) | 0;
  let t = Math.imul(state ^ (state >>> 15), 1 | state);
  t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
  return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
}
function normal() {
  const u = 1 - rand(), v = rand();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}
const uniform = (a, b) => a + (b - a) * rand();
const pick = (arr) => arr[Math.floor(rand() * arr.length)];

const N = 4;
const rows = [];
let acc = 1000;
const nextId = () => `SYN${String(acc++).padStart(5, "0")}`;

/** log2 intensities for bait and control, then raw values with zeros for missing. */
function addProtein({ gene, name, base, diff, sdRep = 0.35, ctrlMissing = 0, baitMissing = 0, reverse = false, contaminant = false, idPrefix = "" }) {
  const id = nextId();
  const bait = [], ctrl = [];
  const baitLevel = base + diff / 2, ctrlLevel = base - diff / 2;
  for (let i = 0; i < N; i++) {
    bait.push(baitLevel + sdRep * normal());
    ctrl.push(ctrlLevel + sdRep * normal());
  }
  const raw = (v) => Math.round(2 ** v);
  const baitRaw = bait.map(raw), ctrlRaw = ctrl.map(raw);
  if (ctrlMissing >= N) ctrlRaw.fill(0);
  else for (let k = 0; k < ctrlMissing; k++) ctrlRaw[(k + Math.floor(rand() * N)) % N] = 0;
  for (let k = 0; k < baitMissing; k++) baitRaw[(k + Math.floor(rand() * N)) % N] = 0;
  const peptides = Math.max(2, Math.round(2 + (base - 18) * 1.2 + 3 * rand()));
  const unique = Math.max(1, peptides - Math.floor(rand() * 3));
  const prefix = reverse ? "REV__" : contaminant ? "CON__" : "";
  rows.push({
    ids: `${prefix}${idPrefix}${id}`,
    gene: reverse ? "" : gene,
    name: reverse ? "" : name,
    peptides, unique,
    coverage: (Math.min(80, 2 + peptides * 2.2 + 8 * rand())).toFixed(1),
    reverse: reverse ? "+" : "",
    contaminant: contaminant ? "+" : "",
    site: rand() < 0.02 ? "+" : "",
    bait: baitRaw, ctrl: ctrlRaw,
  });
}

// Fictional bait and true interactors: strong enrichment, low variance.
addProtein({ gene: "BAITX", name: "Bait protein X (simulated)", base: 31, diff: 9, sdRep: 0.25 });
const trueNames = ["RLN1A", "RLN1B", "TPC2", "TPC4", "ASM3", "ASM7", "KDR5", "CHP9", "NUF2L", "SRB4", "PMC1", "HELX3"];
trueNames.forEach((g, i) => {
  addProtein({
    gene: g, name: `Complex subunit ${g} (simulated)`,
    base: uniform(24, 30), diff: 4.5 + (i / (trueNames.length - 1)) * 4, sdRep: 0.3,
    ctrlMissing: i % 4 === 3 ? 2 : 0,
  });
});

// Weak or borderline binders.
const weak = [
  ["WKB1", 3.4, 0.7], ["WKB2", 3.2, 0.9], ["WKB3", 2.6, 0.5], ["WKB4", 3.8, 1.1], ["WKB5", 2.2, 0.4],
  ["WKB6", 3.1, 0.6], ["WKB7", 1.9, 0.5], ["WKB8", 4.0, 1.3], ["WKB9", 2.9, 0.8], ["WKB10", 3.5, 1.0],
];
for (const [g, diff, sdRep] of weak) {
  addProtein({ gene: g, name: `Weak associated protein ${g} (simulated)`, base: uniform(22, 28), diff, sdRep });
}

// Proteins seen in most bait replicates and in no control at all (missing by design, not by chance).
// The app lists these separately instead of testing them on filled-in control values.
for (let i = 1; i <= 6; i++) {
  addProtein({
    gene: `BND${i}`, name: `Bait only protein BND${i} (simulated)`, base: uniform(24, 29), diff: 7, sdRep: 0.3,
    ctrlMissing: 4, baitMissing: i % 3 === 0 ? 1 : 0,
  });
}

// Background: no enrichment, noise around zero, intensities 1e6 to 1e10.
const stems = ["ZNF", "MRP", "RPL", "RPS", "HSP", "DNJ", "TUB", "ACT", "PRX", "GST", "EIF", "SEC", "VPS", "NDU", "COX", "ATP", "PSM", "UBE", "CDK", "MAP"];
const used = new Set(["BAITX", ...trueNames, ...weak.map((w) => w[0]), ...Array.from({ length: 6 }, (_, i) => `BND${i + 1}`)]);
let bg = 0;
while (bg < 554) {
  const gene = `${pick(stems)}${Math.floor(1 + rand() * 40)}${pick(["", "", "A", "B", "L"])}`;
  if (used.has(gene)) continue;
  used.add(gene);
  const base = uniform(20, 33);
  const lowAbundance = base < 23;
  addProtein({
    gene, name: `Simulated protein ${gene}`, base, diff: 0.5 * normal(), sdRep: 0.45,
    // A handful of low intensity proteins lack control values (left censored), some lack bait values.
    ctrlMissing: lowAbundance && rand() < 0.12 ? 1 + Math.floor(rand() * 2) : 0,
    baitMissing: lowAbundance && rand() < 0.1 ? 1 : 0,
  });
  bg++;
}

// Reverse hits and contaminants (the app removes both).
for (let i = 0; i < 8; i++) addProtein({ gene: "", name: "", base: uniform(20, 26), diff: 0.4 * normal(), reverse: true });
const con = ["KRT1", "KRT2", "KRT9", "KRT10", "TRYP", "ALBU", "KRT14", "CASB"];
con.forEach((g, i) => addProtein({
  gene: g, name: `Contaminant ${g} (simulated)`, base: uniform(24, 31), diff: i < 3 ? 4.5 : 0.3 * normal(),
  sdRep: 0.3, contaminant: true,
}));

// Shuffle so the special rows are not grouped at the top.
for (let i = rows.length - 1; i > 0; i--) {
  const j = Math.floor(rand() * (i + 1));
  [rows[i], rows[j]] = [rows[j], rows[i]];
}

const header = [
  "Protein IDs", "Majority protein IDs", "Gene names", "Protein names", "Peptides", "Unique peptides",
  "Sequence coverage [%]", "Only identified by site", "Reverse", "Potential contaminant",
  ...Array.from({ length: N }, (_, i) => `LFQ intensity BAIT_${i + 1}`),
  ...Array.from({ length: N }, (_, i) => `LFQ intensity CTRL_${i + 1}`),
];
const lines = [header.join("\t")];
for (const r of rows) {
  lines.push([
    r.ids, r.ids, r.gene, r.name, r.peptides, r.unique, r.coverage, r.site, r.reverse, r.contaminant,
    ...r.bait, ...r.ctrl,
  ].join("\t"));
}
const out = resolve(dirname(fileURLToPath(import.meta.url)), "../public/example/proteinGroups_example.txt");
mkdirSync(dirname(out), { recursive: true });
writeFileSync(out, lines.join("\n") + "\n");
console.log(`wrote ${rows.length} rows to ${out}`);
