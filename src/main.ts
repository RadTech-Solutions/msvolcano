import Plotly from "plotly.js-cartesian-dist-min";
import "./style.css";
import { aboutHtml } from "./about";
import { parseDelimited, sniffSeparator, type Table } from "./engine/parse";
import { computeStoichiometry, parsePeptideTable } from "./engine/stoich";
import { benjaminiHochberg, permutationFdr, replicateCorrelation, type PermutationResult } from "./engine/experimental";
import { computeVolcano, cutoffCurve, defaultParams, lfqColumns, toCsv, type Params, type VolcanoResult } from "./engine/volcano";

const MAX_LABELS = 150;

document.querySelector<HTMLDivElement>("#app")!.innerHTML = `
<header>
  <h1>msVolcano</h1>
  <p>Volcano plots for label free interactomics data from MaxQuant. Everything runs in your browser; your file is never uploaded.</p>
  <nav><a href="#tool" id="nav-tool">Tool</a><a href="#about" id="nav-about">About, citation and license</a></nav>
</header>
<section id="about-view" hidden>${aboutHtml}</section>
<div class="layout" id="tool-view">
  <aside>
    <h2>Input</h2>
    <input id="file" type="file" accept=".txt,.tsv,.csv" />
    <div class="actions"><button id="example" type="button">Load simulated example</button></div>
    <label>Statistical test</label>
    <select id="test"><option value="student">Student t-test</option><option value="welch">Welch t-test</option></select>
    <h2>LFQ columns</h2>
    <label>Bait (at least 2)</label><div id="bait" class="cols"></div>
    <label>Control (at least 2)</label><div id="control" class="cols"></div>
    <h2>Cutoff and imputation</h2>
    <div id="sliders"></div>
    <h2>Stoichiometry (optional)</h2>
    <label><input id="stoich" type="checkbox" style="width:auto"> Estimate stoichiometry relative to the bait</label>
    <div id="stoichOpts" hidden>
      <label>Organism</label><select id="organism"></select>
      <label>Bait protein (gene name or UniProt accession)</label><input id="stoichBait" type="text" placeholder="BAITX" />
      <p class="note">iBAQ style, as in the 2016 paper. Tryptic peptide counts come from UniProt Swiss-Prot and follow the standard rule of no cleavage before proline, so values can differ slightly from the original tool.</p>
    </div>
    <details id="exp" class="exp">
      <summary>Experimental <span class="badge">Experimental</span></summary>
      <p class="note">New in version 2, not part of the 2016 publication, not peer reviewed or benchmarked. Use with care and compare with an established tool.</p>
      <label>Cutoff method</label>
      <select id="mode">
        <option value="hyperbola">Hyperbolic curve (as published)</option>
        <option value="perm">s0 score with permutation FDR (Perseus style)</option>
        <option value="bh">Benjamini-Hochberg q value</option>
      </select>
      <div id="expOpts"></div>
      <p class="note">Permutation FDR assumes samples are interchangeable. With few replicates there are only a handful of label shuffles (4 vs 4 gives 69), and it can call more proteins than the published curve or Benjamini-Hochberg when many background proteins differ slightly between groups.</p>
      <label>Missing value imputation</label>
      <select id="imputation">
        <option value="normal">Shifted normal (as published)</option>
        <option value="mindet">MinDet: low quantile of each column</option>
        <option value="minprob">MinProb: random draw around low quantile</option>
      </select>
    </details>
    <h2>Plot</h2>
    <div id="plotOpts"></div>
    <label>Label proteins (gene names, separated by ;)</label>
    <input id="manual" type="text" placeholder="Wdr5;Mll2" />
    <label>Plot title</label><input id="title" type="text" />
    <label>Bait name</label><input id="bait-name" type="text" />
    <div class="actions">
      <button id="png" disabled>PNG</button>
      <button id="svg" disabled>SVG</button>
      <button id="csv" disabled>Hits CSV</button>
    </div>
  </aside>
  <main>
    <div id="msg" class="msg">Choose a MaxQuant proteinGroups.txt (or a Perseus export with LFQ columns).</div>
    <div id="plot"></div>
    <details id="qc" class="exp" hidden><summary>Replicate correlation <span class="badge">Experimental</span></summary><div id="qcTable"></div></details>
    <div id="hits"></div>
  </main>
</div>
<footer>
  If you use msVolcano, please cite Singh et al., <a href="https://doi.org/10.1002/pmic.201600167">Proteomics 2016, 16(18):2491</a>.
  Free for noncommercial use under the PolyForm Noncommercial License 1.0.0.
</footer>`;

function route() {
  const about = location.hash === "#about";
  document.getElementById("about-view")!.hidden = !about;
  document.getElementById("tool-view")!.hidden = about;
  document.getElementById("nav-about")!.classList.toggle("on", about);
  document.getElementById("nav-tool")!.classList.toggle("on", !about);
  if (!about) window.dispatchEvent(new Event("resize"));
}
window.addEventListener("hashchange", route);
route();

const $ = <T extends HTMLElement>(id: string) => document.getElementById(id) as T;

interface Slider { key: keyof typeof state.num; label: string; min: number; max: number; step: number }
const state = {
  table: null as Table | null,
  fileName: "",
  fitAxes: false,
  perm: null as PermutationResult | null,
  stoich: new Map<number, number>(),
  stoichError: "",
  result: null as VolcanoResult | null,
  params: null as Params | null,
  num: { minFoldChange: defaultParams.minFoldChange, curvature: defaultParams.curvature, shift: defaultParams.shift, shrink: defaultParams.shrink, xMax: 15, yMax: 5, labelSize: 11, s0: 0.1, fdr: 0.05, perms: 250, qcut: 0.05 },
};

const sliders: Slider[] = [
  { key: "minFoldChange", label: "minFoldChange", min: 0, max: 10, step: 0.1 },
  { key: "curvature", label: "curvature", min: 0, max: 10, step: 0.1 },
  { key: "shift", label: "shift", min: 0, max: 5, step: 0.1 },
  { key: "shrink", label: "shrink", min: 0, max: 5, step: 0.1 },
];
const expSliders: Slider[] = [
  { key: "s0", label: "s0 (permutation mode)", min: 0, max: 2, step: 0.05 },
  { key: "fdr", label: "FDR (permutation mode)", min: 0.01, max: 0.25, step: 0.01 },
  { key: "perms", label: "Permutations", min: 20, max: 1000, step: 10 },
  { key: "qcut", label: "q value cutoff (BH mode)", min: 0.001, max: 0.25, step: 0.001 },
];
const plotSliders: Slider[] = [
  { key: "xMax", label: "x limit", min: 1, max: 25, step: 0.5 },
  { key: "yMax", label: "y limit", min: 1, max: 25, step: 0.5 },
  { key: "labelSize", label: "Label size", min: 6, max: 24, step: 1 },
];

const setters: Partial<Record<keyof typeof state.num, (v: number) => void>> = {};

function mountSliders(host: HTMLElement, list: Slider[]) {
  for (const s of list) {
    const wrap = document.createElement("div");
    wrap.innerHTML = `<label>${s.label}</label><div class="row"><input type="range" min="${s.min}" max="${s.max}" step="${s.step}" value="${state.num[s.key]}"><input type="number" min="${s.min}" max="${s.max}" step="${s.step}" value="${state.num[s.key]}"></div>`;
    const [range, num] = wrap.querySelectorAll("input");
    const set = (v: string) => { const n = Number(v); if (Number.isFinite(n)) { state.num[s.key] = n; range.value = num.value = String(n); schedule(["shift", "shrink", "s0", "fdr", "perms"].includes(s.key) ? "compute" : "draw"); } };
    setters[s.key] = (v) => { state.num[s.key] = v; range.value = num.value = String(v); };
    range.addEventListener("input", () => set(range.value));
    num.addEventListener("change", () => set(num.value));
    host.append(wrap);
  }
}
mountSliders($("sliders"), sliders);
mountSliders($("plotOpts"), plotSliders);
mountSliders($("expOpts"), expSliders);

function checkboxes(host: HTMLElement, names: string[], preset: string[]) {
  host.innerHTML = "";
  for (const n of names) {
    const l = document.createElement("label");
    const c = document.createElement("input");
    c.type = "checkbox"; c.value = n; c.checked = preset.includes(n);
    c.addEventListener("change", () => schedule("compute"));
    l.append(c, document.createTextNode(n.replace(/^LFQ intensity /i, "")));
    host.append(l);
  }
}
const checked = (host: HTMLElement) => Array.from(host.querySelectorAll<HTMLInputElement>("input:checked")).map((c) => c.value);

const peptideCache = new Map<number, Promise<ReturnType<typeof parsePeptideTable>>>();
const loadPeptides = (taxid: number) => {
  if (!peptideCache.has(taxid)) {
    peptideCache.set(taxid, fetch(`peptides/${taxid}.tsv`).then((r) => { if (!r.ok) throw new Error(`peptide table ${taxid}: HTTP ${r.status}`); return r.text(); }).then(parsePeptideTable));
  }
  return peptideCache.get(taxid)!;
};
fetch("peptides/index.json").then((r) => r.json()).then((list: { taxid: number; name: string; proteins: number }[]) => {
  $("organism").innerHTML = list.map((o) => `<option value="${o.taxid}">${o.name} (${o.proteins} proteins)</option>`).join("");
}).catch(() => { /* stoichiometry stays unavailable without the tables */ });
$("stoich").addEventListener("change", () => { $("stoichOpts").hidden = !($("stoich") as HTMLInputElement).checked; schedule("compute"); });
["organism", "stoichBait"].forEach((id) => $(id).addEventListener("change", () => schedule("compute")));

function setMsg(text: string, err = false) { const m = $("msg"); m.textContent = text; m.className = err ? "msg err" : "msg"; }

let timer = 0;
let needsCompute = false;
// A pending compute also redraws, so a later draw request must not cancel it.
function schedule(kind: "compute" | "draw") {
  if (kind === "compute") needsCompute = true;
  clearTimeout(timer);
  timer = window.setTimeout(() => {
    const full = needsCompute;
    needsCompute = false;
    if (full) compute(); else draw();
  }, 120);
}
["test", "imputation", "mode", "manual", "title", "bait-name"].forEach((id) => $(id).addEventListener("input", () => schedule(["test", "imputation", "mode"].includes(id) ? "compute" : "draw")));

function loadText(text: string, name: string, preset?: { bait: RegExp; control: RegExp; baitName: string }) {
  const table = parseDelimited(text, sniffSeparator(text));
  const lfq = lfqColumns(table);
  if (lfq.length < 4) { setMsg(`Found ${lfq.length} columns with "LFQ" in the name. At least 4 are needed (2 bait, 2 control).`, true); return; }
  state.table = table; state.fileName = name; state.fitAxes = true; state.result = null;
  checkboxes($("bait"), lfq, preset ? lfq.filter((c) => preset.bait.test(c)) : []);
  checkboxes($("control"), lfq, preset ? lfq.filter((c) => preset.control.test(c)) : []);
  ($("bait-name") as HTMLInputElement).value = preset?.baitName ?? "";
  setMsg(`${table.rows.length} rows, ${lfq.length} LFQ columns. Select bait and control columns.`);
  if (preset) compute();
}

$("file").addEventListener("change", async (e) => {
  const f = (e.target as HTMLInputElement).files?.[0];
  if (!f) return;
  setMsg(`Reading ${f.name}...`);
  loadText(await f.text(), f.name);
});

$("example").addEventListener("click", async () => {
  try {
    const res = await fetch("example/proteinGroups_example.txt");
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    loadText(await res.text(), "proteinGroups_example.txt", { bait: /BAIT_/, control: /CTRL_/, baitName: "BAITX" });
  } catch (err) { setMsg(`Could not load the example: ${(err as Error).message}`, true); }
});

function currentParams(): Params {
  return {
    bait: checked($("bait")), control: checked($("control")), welch: ($("test") as HTMLSelectElement).value === "welch",
    shift: state.num.shift, shrink: state.num.shrink, minFoldChange: state.num.minFoldChange, curvature: state.num.curvature, seed: defaultParams.seed,
    imputation: ($("imputation") as HTMLSelectElement).value as Params["imputation"],
  };
}

function compute() {
  if (!state.table) return;
  const p = currentParams();
  if (p.bait.length < 2 || p.control.length < 2) { setMsg("Select at least two bait and two control columns."); return; }
  const overlap = p.bait.filter((b) => p.control.includes(b));
  if (overlap.length) { setMsg(`A column is selected as both bait and control: ${overlap[0]}`, true); return; }
  try {
    state.result = computeVolcano(state.table, p);
    state.params = p;
    state.perm = null;
    const mode = ($("mode") as HTMLSelectElement).value;
    if (mode === "perm") {
      state.perm = permutationFdr(state.result.proteins, p.bait.length, state.num.s0, state.num.fdr, state.num.perms, p.seed);
    }
    renderQc(state.result);
    void updateStoich();
    if (state.fitAxes) {
      // Fit the axes to the data once per file; the user can still adjust the limits afterwards.
      const xs = state.result.proteins.map((q) => Math.abs(q.x)), ys = state.result.proteins.map((q) => q.y);
      setters.xMax?.(Math.min(25, Math.max(5, Math.ceil(Math.max(...xs)) + 1)));
      setters.yMax?.(Math.min(25, Math.max(3, Math.ceil(Math.max(...ys)) + 1)));
      state.fitAxes = false;
    }
    draw();
  } catch (err) { setMsg((err as Error).message, true); }
}

async function updateStoich() {
  state.stoich = new Map(); state.stoichError = "";
  if (!($("stoich") as HTMLInputElement).checked || !state.result) return;
  const key = ($("stoichBait") as HTMLInputElement).value.trim();
  if (!key) { state.stoichError = "Enter the bait protein to estimate stoichiometry."; draw(); return; }
  try {
    const peptides = await loadPeptides(Number(($("organism") as HTMLSelectElement).value));
    state.stoich = computeStoichiometry(state.result.proteins, peptides, key);
  } catch (err) { state.stoichError = (err as Error).message; }
  draw();
}

function draw() {
  const r = state.result, table = state.table;
  if (!r || !table) return;
  const p = { ...state.params!, minFoldChange: state.num.minFoldChange, curvature: state.num.curvature };
  // Significance depends on the cutoff only, so recompute the flag without redoing the statistics.
  const mode = ($("mode") as HTMLSelectElement).value;
  const qv = mode === "bh" ? benjaminiHochberg(r.proteins.map((q) => 10 ** -q.y)) : [];
  const prot = r.proteins.map((q, i) => ({
    ...q,
    q: qv[i],
    significant:
      mode === "perm" ? !!state.perm?.significant.has(i)
      : mode === "bh" ? q.x > 0 && qv[i] < state.num.qcut
      : q.x > p.minFoldChange && p.curvature / (q.x - p.minFoldChange) < q.y,
  }));
  const hits = prot.filter((q) => q.significant).sort((a, b) => b.x - a.x);
  const { xMax, yMax, labelSize } = state.num;
  const rest = prot.filter((q) => !q.significant);
  const curve = cutoffCurve(p.minFoldChange, p.curvature, xMax, yMax);
  const manual = ($("manual") as HTMLInputElement).value.split(";").map((s) => s.trim().toLowerCase()).filter(Boolean);
  const manualHits = manual.length ? prot.filter((q) => manual.some((m) => q.gene.toLowerCase().startsWith(m))) : [];
  const labelled = hits.length <= MAX_LABELS ? hits.filter((q) => !manualHits.includes(q)) : [];

  const hitSize = (q: { row: number }) => (state.stoich.size ? Math.min(30, 5 + 3 * Math.log2((state.stoich.get(q.row) ?? 0) + 1.5)) : 6);
  const trace = (name: string, pts: typeof prot, color: string, size: number | number[], text?: boolean) => ({
    type: "scattergl", mode: text ? "markers+text" : "markers", name,
    x: pts.map((q) => q.x), y: pts.map((q) => q.y),
    text: pts.map((q) => q.gene), textposition: "top right", textfont: { size: labelSize, color },
    customdata: pts.map((q) => q.id), hovertemplate: "%{text}<br>%{customdata}<br>diff %{x:.2f}<br>-log10 p %{y:.2f}<extra></extra>",
    marker: { color, size, opacity: text ? 0.9 : 0.55 },
    ...(text ? {} : { textfont: undefined, mode: "markers" }),
  });
  const line = (pts: { x: number; y: number }[]) => ({ type: "scatter", mode: "lines", x: pts.map((q) => q.x), y: pts.map((q) => q.y), line: { color: "#888", dash: "dash", width: 1 }, hoverinfo: "skip", showlegend: false });

  const title = ($("title") as HTMLInputElement).value;
  const baitName = ($("bait-name") as HTMLInputElement).value;
  const css = getComputedStyle(document.documentElement);
  const fg = css.getPropertyValue("--fg").trim(), grid = css.getPropertyValue("--line").trim();
  const layout = {
    font: { color: fg },
    title: { text: title + (baitName ? `<br><sup>Bait: ${baitName}</sup>` : ""), x: 0.5 },
    xaxis: { title: "t-test difference (log2 LFQ)", range: [-xMax, xMax], zeroline: false, gridcolor: grid, linecolor: grid },
    yaxis: { title: "-log10 p value", range: [0, yMax], zeroline: false, gridcolor: grid, linecolor: grid },
    showlegend: false, hovermode: "closest", margin: { t: 60, r: 20, b: 55, l: 60 },
    paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
  };
  Plotly.react($("plot"), [
    trace("other", rest, "#1f5fbf", 5),
    trace("hits", hits, "#d62728", hits.map(hitSize), false),
    { ...trace("labels", labelled, "#d62728", 6, true), marker: { color: "#d62728", size: labelled.map(hitSize) } },
    { ...trace("selected", manualHits, "#2ca02c", 9, true), marker: { color: "#2ca02c", size: 9, line: { color: "#fff", width: 1 } } },
    ...(mode === "hyperbola" ? [line(curve.right), line(curve.left)] : []),
  ], layout, { responsive: true, displaylogo: false });

  const stNote = state.stoichError ? ` Stoichiometry: ${state.stoichError}` : "";
  const note = hits.length > MAX_LABELS ? ` Too many hits to label (max ${MAX_LABELS}); raise minFoldChange or curvature.` : "";
  const d = r.dropped;
  const modeNote = mode === "perm" ? (state.perm ? ` Permutation FDR ${state.perm.fdr} (${state.perm.permutations} permutations, d threshold ${state.perm.threshold.toFixed(2)}), experimental.` : " Not enough replicates for permutation FDR.")
    : mode === "bh" ? ` Benjamini-Hochberg q < ${state.num.qcut}, experimental.` : "";
  setMsg(`${hits.length} significant of ${prot.length} proteins. Removed ${d.contaminantOrReverse} contaminant/reverse, ${d.absentInBait} absent in bait.${r.logTransformed ? " Values were log2 transformed." : ""}${modeNote}${note}${stNote}`);
  $("hits").innerHTML = hits.length
    ? `<table><thead><tr><th>Gene</th><th>Protein ID</th><th>Difference</th><th>-log10 p</th>${mode === "bh" ? "<th>q value</th>" : ""}${state.stoich.size ? "<th>Stoichiometry</th>" : ""}</tr></thead><tbody>${hits.slice(0, 200).map((q) => `<tr><td>${esc(q.gene)}</td><td>${esc(q.id)}</td><td>${q.x.toFixed(2)}</td><td>${q.y.toFixed(2)}</td>${mode === "bh" ? `<td>${q.q.toExponential(1)}</td>` : ""}${state.stoich.size ? `<td>${fmtSt(state.stoich.get(q.row))}</td>` : ""}</tr>`).join("")}</tbody></table>`
    : "";
  for (const id of ["png", "svg", "csv"]) ($(id) as HTMLButtonElement).disabled = false;
  (window as unknown as { __hits: typeof hits }).__hits = hits;
}

function renderQc(r: VolcanoResult) {
  const p = state.params!;
  const names = [...p.bait, ...p.control].map((n) => n.replace(/^LFQ intensity /i, "")).map((n) => (n.length > 18 ? n.slice(0, 8) + "..." + n.slice(-7) : n));
  const c = replicateCorrelation(r.proteins);
  const shade = (v: number) => `background:rgba(31,95,191,${Math.max(0, (v - 0.5) * 2).toFixed(2)})`;
  $("qcTable").innerHTML = `<table><thead><tr><th></th>${names.map((n) => `<th>${esc(n)}</th>`).join("")}</tr></thead><tbody>${c.map((row, i) => `<tr><th>${esc(names[i])}</th>${row.map((v) => `<td style="${shade(v)}">${v.toFixed(2)}</td>`).join("")}</tr>`).join("")}</tbody></table><p class="note">Pearson correlation of imputed log2 LFQ values. Replicates of one condition should correlate more strongly with each other than with the other group.</p>`;
  $("qc").hidden = false;
}

const fmtSt = (v: number | undefined) => (v === undefined || !Number.isFinite(v) ? "" : v.toPrecision(3));

const esc = (s: string) => s.replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]!);
const stem = () => state.fileName.replace(/\.[^.]+$/, "") || "msvolcano";

$("png").addEventListener("click", () => Plotly.downloadImage($("plot"), { format: "png", filename: stem(), width: 1400, height: 1000 }));
$("svg").addEventListener("click", () => Plotly.downloadImage($("plot"), { format: "svg", filename: stem(), width: 1400, height: 1000 }));
$("csv").addEventListener("click", () => {
  if (!state.table || !state.params) return;
  const hits = (window as unknown as { __hits: VolcanoResult["proteins"] }).__hits;
  const url = URL.createObjectURL(new Blob([toCsv(state.table, hits, state.params, state.stoich)], { type: "text/csv" }));
  const a = Object.assign(document.createElement("a"), { href: url, download: `${stem()}_hits.csv` });
  a.click(); URL.revokeObjectURL(url);
});
