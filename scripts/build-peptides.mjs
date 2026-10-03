// Builds public/peptides/<taxid>.tsv and index.json from UniProt reviewed FASTA.
// Usage: node scripts/build-peptides.mjs [taxid ...]
import { gunzipSync } from "node:zlib";
import { mkdirSync, writeFileSync, readFileSync, existsSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";

const ORGANISMS = [
  [9606, "Human"], [10090, "Mouse"], [10116, "Rat"], [559292, "Yeast (S. cerevisiae S288C)"],
  [7227, "Drosophila melanogaster"], [3702, "Arabidopsis thaliana"], [83333, "E. coli K12"],
  [7955, "Zebrafish"], [6239, "C. elegans"],
];
const out = join(dirname(fileURLToPath(import.meta.url)), "..", "public", "peptides");

// Keep identical to countTrypticPeptides in src/engine/stoich.ts.
function countTrypticPeptides(seq) {
  let count = 0, start = 0;
  for (let i = 0; i < seq.length; i++) {
    const c = seq[i];
    if (i === seq.length - 1 || ((c === "K" || c === "R") && seq[i + 1] !== "P")) {
      const len = i + 1 - start;
      if (len >= 7 && len <= 30) count++;
      start = i + 1;
    }
  }
  return count;
}

async function download(taxid) {
  const url = `https://rest.uniprot.org/uniprotkb/stream?query=(organism_id:${taxid})+AND+(reviewed:true)&format=fasta&compressed=true`;
  let err;
  for (let attempt = 1; attempt <= 4; attempt++) {
    try {
      const res = await fetch(url);
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      let buf = Buffer.from(await res.arrayBuffer());
      if (buf[0] === 0x1f && buf[1] === 0x8b) buf = gunzipSync(buf);
      return { text: buf.toString("utf8"), release: res.headers.get("x-uniprot-release") };
    } catch (e) {
      err = e;
      await new Promise((r) => setTimeout(r, 2000 * attempt));
    }
  }
  throw err;
}

function parseFasta(fasta) {
  const rows = [];
  for (const rec of fasta.split(/^>/m)) {
    if (!rec.trim()) continue;
    const nl = rec.indexOf("\n");
    const header = rec.slice(0, nl);
    const seq = rec.slice(nl + 1).replace(/\s+/g, "");
    const acc = header.split("|")[1];
    const gene = /\bGN=(\S+)/.exec(header)?.[1] ?? "";
    if (acc) rows.push(`${acc}\t${gene}\t${countTrypticPeptides(seq)}`);
  }
  return rows;
}

mkdirSync(out, { recursive: true });
const indexPath = join(out, "index.json");
const prev = existsSync(indexPath) ? JSON.parse(readFileSync(indexPath, "utf8")) : [];
const wanted = process.argv.slice(2).map(Number);
const todo = ORGANISMS.filter(([t]) => !wanted.length || wanted.includes(t));
const index = new Map(prev.map((e) => [e.taxid, e]));
const today = new Date().toISOString().slice(0, 10);
for (const [taxid, name] of todo) {
  try {
    const { text, release } = await download(taxid);
    const rows = parseFasta(text);
    writeFileSync(join(out, `${taxid}.tsv`), "accession\tgene\ttryptic_peptides\n" + rows.join("\n") + "\n");
    index.set(taxid, { taxid, name, proteins: rows.length, source: `UniProt reviewed, ${release ?? today}` });
    console.log(`ok   ${taxid} ${name}: ${rows.length}`);
  } catch (e) {
    console.error(`FAIL ${taxid} ${name}: ${e.message}`);
  }
}
const order = ORGANISMS.map(([t]) => t);
writeFileSync(indexPath, JSON.stringify([...index.values()].sort((a, b) => order.indexOf(a.taxid) - order.indexOf(b.taxid)), null, 1) + "\n");
