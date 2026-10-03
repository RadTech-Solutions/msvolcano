import "@fontsource-variable/geist";
import "@fontsource-variable/geist-mono";
import Plotly from "plotly.js-cartesian-dist-min";
import "./style.css";
import { aboutHtml } from "./about";
import { benjaminiHochberg, permutationFdr, replicateCorrelation, type PermutationResult } from "./engine/experimental";
import { parseDelimited, sniffSeparator, type Table } from "./engine/parse";
import { computeStoichiometry, parsePeptideTable } from "./engine/stoich";
import { computeVolcano, cutoffCurve, defaultParams, lfqColumns, toCsv, type Params, type Protein, type VolcanoResult } from "./engine/volcano";
import { groupColumns, suggestRole, type Group, type Role } from "./groups";

const MAX_LABELS = 150;

const svg = (body: string, extra = "") => `<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" ${extra}>${body}</svg>`;
const ICON = {
  sun: svg('<circle cx="12" cy="12" r="4"/><path d="M12 2v2M12 20v2M4.9 4.9l1.4 1.4M17.7 17.7l1.4 1.4M2 12h2M20 12h2M4.9 19.1l1.4-1.4M17.7 6.3l1.4-1.4"/>'),
  moon: svg('<path d="M21 12.8A9 9 0 1 1 11.2 3a7 7 0 0 0 9.8 9.8z"/>'),
  auto: svg('<rect x="3" y="4" width="18" height="12" rx="2"/><path d="M8 20h8M12 16v4"/>'),
  github: svg('<path d="M9 19c-4.3 1.4-4.3-2.5-6-3m12 5v-3.5c0-1 .1-1.4-.5-2 2.8-.3 5.5-1.4 5.5-6a4.6 4.6 0 0 0-1.3-3.2 4.2 4.2 0 0 0-.1-3.2s-1.1-.3-3.5 1.3a12.3 12.3 0 0 0-6.2 0C6.5 2.8 5.4 3.1 5.4 3.1a4.2 4.2 0 0 0-.1 3.2A4.6 4.6 0 0 0 4 9.5c0 4.6 2.7 5.7 5.5 6-.6.6-.6 1.2-.5 2V21"/>'),
  upload: svg('<path d="M12 16V4M7 9l5-5 5 5M5 20h14"/>'),
  download: svg('<path d="M12 4v12M7 11l5 5 5-5M5 20h14"/>'),
  copy: svg('<rect x="9" y="9" width="11" height="11" rx="2"/><path d="M5 15V6a2 2 0 0 1 2-2h9"/>'),
  chev: svg('<path d="m9 6 6 6-6 6"/>', 'class="chev"'),
  lock: svg('<rect x="4" y="10" width="16" height="10" rx="2"/><path d="M8 10V7a4 4 0 0 1 8 0v3"/>'),
  alert: svg('<path d="M12 9v4M12 17h.01"/><path d="M10.3 3.9 1.8 18a2 2 0 0 0 1.7 3h17a2 2 0 0 0 1.7-3L13.7 3.9a2 2 0 0 0-3.4 0z"/>'),
  info: svg('<circle cx="12" cy="12" r="10"/><path d="M12 16v-4M12 8h.01"/>'),
  flask: svg('<path d="M9 3h6M10 3v6L4.5 19a2 2 0 0 0 1.8 3h11.4a2 2 0 0 0 1.8-3L14 9V3"/>'),
};

const logo = `<svg viewBox="0 0 32 32" aria-hidden="true"><rect width="32" height="32" rx="8" fill="var(--color-accent)"/><path d="M5 25 13 12h6l8 13z" fill="var(--color-on-accent)"/><circle cx="16" cy="7" r="2" fill="var(--plot-hit)"/><circle cx="21.5" cy="9.5" r="1.5" fill="var(--plot-hit)"/><circle cx="10.5" cy="9" r="1.3" fill="var(--plot-hit)"/></svg>`;

function heroArt(): string {
  // Decorative volcano plot drawn from a fixed pseudo random sequence.
  let s = 7;
  const r = () => { s = (s * 16807) % 2147483647; return s / 2147483647; };
  const dots: string[] = [];
  for (let i = 0; i < 150; i++) {
    const x = 200 + (r() + r() + r() - 1.5) * 70;
    const y = 232 - Math.pow(r(), 2.2) * 60 - Math.abs(x - 200) * 0.1;
    dots.push(`<circle cx="${x.toFixed(1)}" cy="${y.toFixed(1)}" r="2.6" fill="currentColor" opacity=".55"/>`);
  }
  const hits = [[286, 60], [320, 96], [262, 120], [350, 70], [300, 150], [334, 140]].map(([x, y]) => `<circle cx="${x}" cy="${y}" r="4.2" fill="var(--plot-hit)"/>`).join("");
  return `<svg class="hero-art" viewBox="0 0 400 260" aria-hidden="true"><path d="M20 244H388M200 244V12" stroke="currentColor" stroke-width="1" opacity=".5"/>
    <path d="M232 14C236 70 250 118 292 150S372 190 388 196M168 14C164 70 150 118 108 150S28 190 12 196" stroke="currentColor" stroke-width="1.4" stroke-dasharray="5 5" fill="none" opacity=".8" transform="translate(0 0)"/>${dots.join("")}${hits}</svg>`;
}

/* ---------- Skeleton ---------- */

interface Slider { key: keyof typeof state.num; label: string; min: number; max: number; step: number; help?: string }
const sliderGroups: Record<string, Slider[]> = {
  cutoff: [
    { key: "minFoldChange", label: "Minimum enrichment (log2)", min: 0, max: 10, step: 0.1, help: "Proteins must be enriched in bait by at least this much. 3 means 8 times more." },
    { key: "curvature", label: "Curvature", min: 0, max: 10, step: 0.1, help: "How tightly the curve hugs the axes. Raise it if you get too many weak hits." },
  ],
  impute: [
    { key: "shift", label: "Fill-in shift", min: 0, max: 5, step: 0.1, help: "Missing values are replaced by small random numbers this many standard deviations below your typical signal." },
    { key: "shrink", label: "Fill-in spread", min: 0, max: 5, step: 0.1, help: "How varied those random numbers are. The defaults suit most data." },
  ],
  plot: [
    { key: "xMax", label: "x limit", min: 1, max: 25, step: 0.5 },
    { key: "yMax", label: "y limit", min: 1, max: 25, step: 0.5 },
    { key: "labelSize", label: "Label size", min: 6, max: 24, step: 1 },
    { key: "labelCount", label: "Labels shown", min: 0, max: 150, step: 1, help: "Strongest interactors are labelled first. The table lists all of them." },
  ],
  exp: [
    { key: "s0", label: "s0 (fudge factor)", min: 0, max: 2, step: 0.05, help: "Stabilises the score for proteins with very small variance." },
    { key: "fdr", label: "FDR", min: 0.01, max: 0.25, step: 0.01, help: "Share of false calls you accept." },
    { key: "perms", label: "Permutations", min: 20, max: 1000, step: 10 },
    { key: "qcut", label: "q value cutoff", min: 0.001, max: 0.25, step: 0.001, help: "Used by the Benjamini-Hochberg method." },
  ],
};

const step = (id: string, n: string, title: string, hint: string, body: string, opts: { open?: boolean; exp?: boolean; optional?: boolean } = {}) => `
<details class="step${opts.exp ? " exp" : ""}" id="${id}"${opts.open ? " open" : ""}>
  <summary><span class="badge-n">${n}</span><span class="step-title">${title}</span>${opts.exp ? '<span class="pill exp">Experimental</span>' : opts.optional ? '<span class="pill">Optional</span>' : ""}<span class="hint" id="${id}-hint">${hint}</span>${ICON.chev}</summary>
  <div class="body">${body}</div>
</details>`;

document.querySelector<HTMLDivElement>("#app")!.innerHTML = `
<header class="topbar">
  <h1 class="sr">msVolcano, volcano plots for interactomics</h1>
  <a class="brand" href="#tool">${logo}<span>msVolcano</span><small>v2</small></a>
  <nav class="tabs" aria-label="Main"><a href="#tool" id="nav-tool">Analyze</a><a href="#about" id="nav-about">About and citation</a></nav>
  <span class="spacer"></span>
  <div class="seg" role="group" aria-label="Color theme" id="themeSeg">
    <button type="button" data-theme-set="light" aria-pressed="false" title="Light" aria-label="Light">${ICON.sun}<span class="lbl">Light</span></button>
    <button type="button" data-theme-set="system" aria-pressed="true" title="Match system" aria-label="Match system">${ICON.auto}<span class="lbl">Auto</span></button>
    <button type="button" data-theme-set="dark" aria-pressed="false" title="Dark" aria-label="Dark">${ICON.moon}<span class="lbl">Dark</span></button>
  </div>
  <a class="icon-link" href="https://github.com/uksurd88/msvolcano" aria-label="Source code on GitHub" title="Source code on GitHub">${ICON.github}</a>
</header>

<div class="workspace" id="tool-view">
  <aside class="sidebar" aria-label="Analysis settings">
    ${step("s1", "1", "Add your data", "MaxQuant file", `
      <label class="drop" id="drop" for="file">${ICON.upload}<strong>Drop proteinGroups.txt here</strong><span>or click to browse (.txt, .tsv, .csv)</span></label>
      <input id="file" class="sr" type="file" accept=".txt,.tsv,.csv" />
      <div class="filechip" id="filechip" hidden></div>
      <div class="row-actions"><button class="btn small" id="example" type="button">${ICON.flask}Try an example dataset (simulated)</button><button class="btn small" id="replace" type="button" hidden>Replace file</button></div>
      <details style="font-size:.82rem;color:var(--color-ink-2)"><summary style="cursor:pointer;font-weight:600">What file do I need?</summary><p style="margin-top:6px">The <code>proteinGroups.txt</code> from your MaxQuant output folder. Switch on <b>LFQ</b> in MaxQuant before running it, so the file has columns named "LFQ intensity ...". Other formats are not supported yet.</p></details>
      <p class="help" style="font-size:.78rem;color:var(--color-ink-3);display:flex;gap:6px;align-items:center"><span style="width:14px;height:14px;display:inline-flex">${ICON.lock}</span>Stays on this computer. Works offline once loaded.</p>`, { open: true })}
    ${step("s2", "2", "Choose bait and control", "Pick two groups", `
      <p class="help" style="font-size:.8rem;color:var(--color-ink-3)">Replicates are grouped by name. Mark one group as bait and one as control. Untick a replicate to leave it out. Several bait groups are pooled into one.</p>
      <div class="groups" id="groups"></div>`)}
    ${step("s3", "3", "Set the cutoff", "Hyperbolic curve", `
      <div class="field"><label for="test">Statistical test</label><select id="test"><option value="student">Student t-test (equal variance)</option><option value="welch">Welch t-test (unequal variance)</option></select></div>
      <div id="sl-cutoff"></div>
      <details style="margin-top:4px"><summary style="cursor:pointer;font-weight:600;font-size:.85rem">Advanced: missing values</summary><div id="sl-impute" style="display:grid;gap:var(--space-3);margin-top:var(--space-3)"></div></details>`)}
    ${step("s4", "4", "Plot and labels", "Axes, names, titles", `
      <div id="sl-plot"></div>
      <div class="field"><label for="manual">Highlight proteins</label><input id="manual" type="text" placeholder="Wdr5;Mll2" autocomplete="off" /><span class="help">Gene names separated by semicolons. Matches names that start with your text.</span></div>
      <div class="field"><label for="title">Plot title</label><input id="title" type="text" autocomplete="off" /></div>
      <div class="field"><label for="bait-name">Bait name</label><input id="bait-name" type="text" autocomplete="off" /></div>`)}
    ${step("s5", "+", "Stoichiometry", "", `
      <label class="check"><input id="stoich" type="checkbox" /> <span>Estimate abundance relative to the bait</span></label>
      <div id="stoichOpts" hidden style="display:grid;gap:12px">
        <div class="field"><label for="organism">Organism</label><select id="organism"></select></div>
        <div class="field"><label for="stoichBait">Bait protein</label><input id="stoichBait" type="text" placeholder="Gene name or UniProt accession" autocomplete="off" /></div>
        <p class="help" style="font-size:.78rem;color:var(--color-ink-3)">iBAQ style, as in the 2016 paper. Peptide counts come from UniProt Swiss-Prot and follow the standard rule of no cleavage before proline, so values can differ from the original tool.</p>
      </div>`, { optional: true })}
    ${step("s6", "+", "Experimental methods", "", `
      <p class="help" style="font-size:.78rem;color:var(--color-ink-3)">New in version 2, not part of the 2016 publication and not peer reviewed. Compare with an established tool before relying on them.</p>
      <div class="field"><label for="mode">Cutoff method</label><select id="mode"><option value="hyperbola">Hyperbolic curve (as published)</option><option value="perm">Permutation FDR (Perseus style)</option><option value="bh">Benjamini-Hochberg q value</option></select></div>
      <div id="sl-exp"></div>
      <p class="help" style="font-size:.78rem;color:var(--color-ink-3)">Permutation FDR assumes samples are interchangeable. With few replicates there are only a handful of label shuffles (4 vs 4 gives 69), and it can call more proteins than the published curve when many background proteins differ slightly between groups.</p>
      <div class="field"><label for="imputation">Missing value imputation</label><select id="imputation"><option value="normal">Shifted normal (as published)</option><option value="mindet">MinDet: low quantile of each column</option><option value="minprob">MinProb: random draw around low quantile</option></select></div>`, { exp: true })}
    <div class="sidebar-foot"><button class="link-btn" id="reset" type="button">Reset cutoff and plot settings</button><button class="link-btn" id="share" type="button">Copy settings link</button></div>
  </aside>

  <main class="main">
    <section class="hero" id="hero">
      <div style="display:grid;gap:var(--space-4)">
        <h2>See what binds your bait</h2>
        <p class="lede">Turn a MaxQuant label free interactomics run into a volcano plot and a ranked list of interactors. Choose bait and control, tune the cutoff, export the figure.</p>
        <div class="cta">
          <button class="btn primary" id="hero-browse" type="button">${ICON.upload}Choose a file</button>
          <button class="btn" id="hero-example" type="button">${ICON.flask}Try an example dataset</button>
        </div>
        <p class="privacy">${ICON.lock}<span>Everything runs in your browser. Your data is never uploaded. You can also drop a file anywhere on this page.</span></p>
      </div>
      ${heroArt()}
      <div class="how">
        <div><b>Add proteinGroups.txt</b><span>Reverse hits and contaminants are removed for you.</span></div>
        <div><b>Pick bait and control</b><span>Replicates are grouped by name. One click per group.</span></div>
        <div><b>Tune and export</b><span>Move the curve until the real interactors separate. Save PNG, SVG or CSV.</span></div>
      </div>
    </section>
    <div class="callout err" id="errbox" role="alert" hidden></div>
    <div class="ghost" id="prompt" hidden><b id="prompt-title">Almost there</b><span id="prompt-text"></span><div class="checklist" id="checklist"></div></div>
    <section id="results" hidden>
      <div class="stats" id="stats"></div>
      <div class="card" style="margin-top:var(--space-4)">
        <header>
          <h2 id="plot-title">Volcano plot</h2><span class="spacer"></span>
          <details class="menu" id="exportMenu">
            <summary class="btn small primary">${ICON.download}Export</summary>
            <div class="pop">
              <button id="png" type="button">${ICON.download}Figure as PNG</button>
              <button id="svg" type="button">${ICON.download}Figure as SVG</button>
              <button id="csv" type="button">${ICON.download}Interactors as CSV</button>
              <button id="copy" type="button">${ICON.copy}Copy gene list</button>
            </div>
          </details>
        </header>
        <div id="plot" role="group" aria-label="Interactive volcano plot. Hover points to see gene names. The interactors are listed in the table below."></div>
        <div class="legend-key" id="legend"><span><i style="background:var(--plot-hit)"></i>Significant</span><span><i style="background:var(--plot-point)"></i>Other proteins</span><span><i class="d" style="background:var(--plot-pick)"></i>Highlighted</span><span id="legend-cutoff"><i class="l"></i>Cutoff</span></div>
        <div class="notes" id="notes"></div>
      </div>
      <div class="card" style="margin-top:var(--space-4)">
        <div class="subtabs" role="tablist" aria-label="Results views">
          <button role="tab" id="tab-hits" aria-selected="true" aria-controls="panel-hits" tabindex="0" type="button">Interactors <span id="hit-count" class="num"></span></button>
          <button role="tab" id="tab-qc" aria-selected="false" aria-controls="panel-qc" tabindex="-1" type="button">Replicates <span class="pill exp" style="margin-left:4px">Experimental</span><span id="qc-dot"></span></button>
        </div>
        <div id="panel-hits" role="tabpanel" aria-labelledby="tab-hits">
          <div class="toolbar">
            <input id="hit-filter" type="search" placeholder="Filter by gene or ID" aria-label="Filter interactors" autocomplete="off" />
            <span class="spacer" style="flex:1"></span>
            <button class="btn small" id="copy2" type="button">${ICON.copy}Copy genes</button>
          </div>
          <div class="tablewrap" id="hits"></div>
        </div>
        <div id="panel-qc" role="tabpanel" aria-labelledby="tab-qc" hidden>
          <div class="tablewrap" style="padding:10px 14px" id="qcTable"></div>
        </div>
      </div>
    </section>
  </main>
</div>

<main id="about-view" hidden>${aboutHtml}</main>
<footer class="foot">If you use msVolcano, please cite <a href="https://doi.org/10.1002/pmic.201600167">Singh, Hein and Stewart, Proteomics 2016</a>. Free for noncommercial use under the PolyForm Noncommercial License 1.0.0.<small>Version ${__APP_VERSION__}</small></footer>
<div class="mbar" id="mbar" hidden><b id="mbar-count"></b><span style="flex:1"></span><button class="btn small" id="mbar-adjust" type="button">Adjust cutoff</button><button class="btn small" id="mbar-plot" type="button">Plot</button></div>
<div class="toast" id="toast" role="status" aria-live="polite" aria-atomic="true"></div>`;

const $ = <T extends HTMLElement>(id: string) => document.getElementById(id) as T;

/* ---------- State ---------- */

const DEFAULT_NUM = { minFoldChange: defaultParams.minFoldChange, curvature: defaultParams.curvature, shift: defaultParams.shift, shrink: defaultParams.shrink, xMax: 15, yMax: 5, labelSize: 11, labelCount: 12, s0: 0.1, fdr: 0.05, perms: 250, qcut: 0.05 };
// Exports always use this palette so a figure saved in dark mode is still print ready.
const LIGHT_PLOT = { ink: "#1f2430", grid: "#e9ecf1", axis: "#5b6577", bg: "#ffffff", point: "#7b879c", hit: "#d94a2b", pick: "#0e8a86", curve: "#5b6577" };

interface Saved { n?: Record<string, number>; r?: Record<string, Role>; x?: string[]; m?: string; t?: string; i?: string; man?: string; ti?: string; b?: string; so?: boolean; sg?: string; sb?: string }

const state = {
  table: null as Table | null,
  fileName: "",
  groups: [] as Group[],
  roles: new Map<string, Role>(),
  excluded: new Set<string>(),
  suggested: new Set<string>(),
  groupFilter: "",
  fitAxes: false,
  result: null as VolcanoResult | null,
  params: null as Params | null,
  perm: null as PermutationResult | null,
  stoich: new Map<number, number>(),
  stoichError: "",
  hits: [] as (Protein & { q?: number })[],
  sort: { key: "x", dir: -1 as 1 | -1 },
  filter: "",
  pending: null as Saved | null,
  rowsShown: 25,
  num: { ...DEFAULT_NUM },
};

/* ---------- Theme ---------- */

type ThemeChoice = "light" | "dark" | "system";
function applyTheme(choice: ThemeChoice, persist: boolean) {
  const root = document.documentElement;
  if (choice === "system") delete root.dataset.theme; else root.dataset.theme = choice;
  document.querySelectorAll<HTMLButtonElement>("[data-theme-set]").forEach((b) => b.setAttribute("aria-pressed", String(b.dataset.themeSet === choice)));
  if (persist) { try { if (choice === "system") localStorage.removeItem("msv-theme"); else localStorage.setItem("msv-theme", choice); } catch { /* storage may be blocked */ } }
  if (state.result) void draw();
}
document.querySelectorAll<HTMLButtonElement>("[data-theme-set]").forEach((b) => b.addEventListener("click", () => applyTheme(b.dataset.themeSet as ThemeChoice, true)));
{
  let saved: ThemeChoice = "system";
  try { const t = localStorage.getItem("msv-theme"); if (t === "light" || t === "dark") saved = t; } catch { /* ignore */ }
  applyTheme(saved, false);
}
window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", () => { if (!document.documentElement.dataset.theme && state.result) void draw(); });

/* ---------- Routing ---------- */

function route() {
  const about = location.hash === "#about";
  $("about-view").hidden = !about;
  $("tool-view").hidden = about;
  const [tool, ab] = [$("nav-tool"), $("nav-about")];
  if (about) { ab.setAttribute("aria-current", "page"); tool.removeAttribute("aria-current"); } else { tool.setAttribute("aria-current", "page"); ab.removeAttribute("aria-current"); }
  if (!about) window.dispatchEvent(new Event("resize"));
  window.scrollTo({ top: 0 });
}
window.addEventListener("hashchange", route);
route();

/* ---------- Helpers ---------- */

const esc = (s: string) => s.replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]!);
const cssVar = (name: string) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const stem = () => (state.fileName.replace(/\.[^.]+$/, "") || "msvolcano") + (($("bait-name") as HTMLInputElement).value.trim() ? `_${($("bait-name") as HTMLInputElement).value.trim().replace(/[^\w.-]+/g, "_")}` : "");
let toastTimer = 0;
function toast(msg: string) {
  const t = $("toast"); t.textContent = msg; t.classList.add("on");
  clearTimeout(toastTimer); toastTimer = window.setTimeout(() => t.classList.remove("on"), 4000);
}
function setError(msg: string, offerExample = false) {
  const e = $("errbox");
  e.hidden = !msg;
  e.innerHTML = msg ? `${ICON.alert}<span>${esc(msg)}${offerExample ? ' <button class="link-btn" id="err-example" type="button">Try the example instead</button>' : ""}</span><button class="link-btn dismiss" id="err-dismiss" type="button" aria-label="Dismiss message">Dismiss</button>` : "";
  if (msg) {
    $("err-dismiss").addEventListener("click", () => setError(""));
    if (offerExample) $("err-example").addEventListener("click", () => { setError(""); void loadExample(); });
  }
}
function markStep(id: string, done: boolean, hint?: string) {
  $(id).classList.toggle("done", done);
  const b = $(id).querySelector<HTMLElement>(".badge-n")!;
  if (!b.dataset.n) b.dataset.n = b.textContent ?? "";
  b.innerHTML = done ? svg('<path d="m5 12 5 5 9-10"/>', 'width="14" height="14"') : b.dataset.n;
  if (hint !== undefined) $(`${id}-hint`).textContent = hint;
}
const fmtNum = (v: number) => (Math.abs(v) >= 100 ? v.toFixed(0) : String(Math.round(v * 10) / 10));
const parseNum = (v: string) => Number(v.trim().replace(",", "."));

/* ---------- Sliders ---------- */

const setters: Partial<Record<keyof typeof state.num, (v: number) => void>> = {};
const RECOMPUTE = ["shift", "shrink", "s0", "fdr", "perms"];
function mountSliders(hostId: string, list: Slider[]) {
  const host = $(hostId);
  for (const s of list) {
    const wrap = document.createElement("div");
    wrap.className = "slider";
    const id = `sl-${s.key}`;
    wrap.innerHTML = `<div class="top"><label for="${id}">${s.label}</label><input type="text" inputmode="decimal" aria-label="${s.label} value" value="${state.num[s.key]}"></div><input id="${id}" type="range" min="${s.min}" max="${s.max}" step="${s.step}" value="${state.num[s.key]}">${s.help ? `<div class="help">${s.help}</div>` : ""}`;
    const num = wrap.querySelector<HTMLInputElement>('input[type="text"]')!;
    const range = wrap.querySelector<HTMLInputElement>('input[type="range"]')!;
    const help = wrap.querySelector<HTMLElement>(".help");
    const paint = () => range.style.setProperty("--fill", `${((Number(range.value) - s.min) / (s.max - s.min)) * 100}%`);
    const apply = (v: number) => {
      state.num[s.key] = v; range.value = num.value = String(v); paint();
      if (s.key === "minFoldChange" && help) help.textContent = `Proteins must be enriched in bait by at least this much. ${fmtNum(v)} means ${fmtNum(2 ** v)} times more.`;
    };
    setters[s.key] = apply;
    apply(state.num[s.key]);
    range.addEventListener("input", () => { apply(Number(range.value)); schedule(RECOMPUTE.includes(s.key) ? "compute" : "draw"); });
    num.addEventListener("change", () => {
      const n = parseNum(num.value);
      if (Number.isFinite(n)) { apply(Math.min(s.max, Math.max(s.min, n))); schedule(RECOMPUTE.includes(s.key) ? "compute" : "draw"); } else num.value = String(state.num[s.key]);
    });
    host.append(wrap);
  }
}
for (const [k, list] of Object.entries(sliderGroups)) mountSliders(`sl-${k}`, list);

/* ---------- Scheduling ---------- */

let timer = 0;
let needsCompute = false;
// A pending compute also redraws, so a later draw request must not cancel it.
function schedule(kind: "compute" | "draw") {
  if (kind === "compute") needsCompute = true;
  clearTimeout(timer);
  timer = window.setTimeout(() => {
    const full = needsCompute;
    needsCompute = false;
    if (full) compute(); else void draw();
  }, 120);
}
for (const id of ["test", "imputation", "mode"]) $(id).addEventListener("input", () => schedule("compute"));
for (const id of ["manual", "title", "bait-name"]) $(id).addEventListener("input", () => schedule("draw"));

/* ---------- Data loading ---------- */

const paint = () => new Promise<void>((r) => requestAnimationFrame(() => setTimeout(r, 0)));

function diagnose(table: Table): string {
  const cols = table.columns.join("\n");
  if (/Protein\.Group/i.test(cols) && !/lfq/i.test(cols)) return "This looks like a DIA-NN report.";
  if (/PG\.ProteinGroups/i.test(cols)) return "This looks like a Spectronaut report.";
  if (/MaxLFQ Intensity/i.test(cols) && /Combined|Protein Probability/i.test(cols)) return "This looks like a FragPipe report.";
  if (/^Intensity\b/im.test(cols) && !/lfq/i.test(cols)) return "This looks like a MaxQuant file without LFQ columns. Run MaxQuant again with LFQ switched on.";
  return "";
}

let prevChip = { html: "", hidden: true };
// A failed load must not leave the user guessing which file the plot belongs to.
function fail(msg: string) {
  const chip = $("filechip");
  chip.innerHTML = prevChip.html; chip.hidden = prevChip.hidden;
  setError(`${state.table ? `Still showing ${state.fileName}. ` : ""}${msg}`, true);
}

async function loadFile(f: File) {
  location.hash === "#about" && (location.hash = "#tool");
  setError("");
  if (f.size > 150e6) toast("This is a large file. Reading it can take a while.");
  prevChip = { html: $("filechip").innerHTML, hidden: $("filechip").hidden };
  $("filechip").hidden = false;
  $("filechip").innerHTML = `<b>${esc(f.name)}</b><span>Reading...</span>`;
  await paint();
  try { loadText(await f.text(), f.name); } catch (err) { fail(`Could not read the file: ${(err as Error).message}`); }
}

function loadText(text: string, name: string, preset?: { baitName: string }) {
  setError("");
  const table = parseDelimited(text, sniffSeparator(text));
  const lfq = lfqColumns(table);
  if (table.rows.length === 0 || table.columns.length < 2) {
    fail("That file has no table in it. msVolcano needs a tab separated MaxQuant proteinGroups.txt.");
    return;
  }
  if (lfq.length < 4) {
    const why = diagnose(table);
    fail(`${why ? why + " " : ""}Found ${lfq.length} columns with "LFQ" in the name. msVolcano needs at least 4 (two bait and two control replicates) from a MaxQuant proteinGroups.txt.`);
    return;
  }
  const groups = groupColumns(lfq);
  if (groups.length < 2) {
    fail(`Only one group of LFQ columns was found (${groups[0].name}). You need at least two: one bait and one control.`);
    return;
  }
  state.table = table; state.fileName = name; state.fitAxes = true; state.result = null;
  state.groups = groups;
  // Only the first obvious bait and the first obvious control are pre-selected, so unrelated groups are never pooled silently.
  const seen = { bait: false, control: false };
  state.roles = new Map(groups.map((g) => { let r = suggestRole(g.name); if (r !== "off") { if (seen[r]) r = "off"; else seen[r] = true; } return [g.name, r]; }));
  state.excluded = new Set();
  state.suggested = new Set(groups.filter((g) => state.roles.get(g.name) !== "off").map((g) => g.name));
  state.groupFilter = "";
  const chip = $("filechip");
  chip.hidden = false;
  chip.innerHTML = `<b>${esc(name)}</b><span class="num">${table.rows.length.toLocaleString()} proteins</span>`;
  ($("bait-name") as HTMLInputElement).value = preset?.baitName ?? ($("bait-name") as HTMLInputElement).value;
  $("replace").hidden = false;
  $("tool-view").classList.add("has-data");
  markStep("s1", true, `${table.rows.length.toLocaleString()} proteins`);
  ($("s1") as HTMLDetailsElement).open = false;
  ($("s2") as HTMLDetailsElement).open = true;
  if (state.pending) applyPendingRoles();
  renderGroups();
  compute();
}

$("file").addEventListener("change", (e) => { const input = e.target as HTMLInputElement; const f = input.files?.[0]; if (f) void loadFile(f).finally(() => { input.value = ""; }); });
const drop = $("drop");
["dragenter", "dragover"].forEach((ev) => drop.addEventListener(ev, (e) => { e.preventDefault(); drop.classList.add("over"); }));
["dragleave", "drop"].forEach((ev) => drop.addEventListener(ev, () => drop.classList.remove("over")));
// Dropping anywhere on the page works.
window.addEventListener("dragover", (e) => { e.preventDefault(); document.body.classList.add("dragging"); });
window.addEventListener("dragleave", (e) => { if (!e.relatedTarget) document.body.classList.remove("dragging"); });
window.addEventListener("drop", (e) => { e.preventDefault(); document.body.classList.remove("dragging"); const f = e.dataTransfer?.files?.[0]; if (f) void loadFile(f); });
$("hero-browse").addEventListener("click", () => $("file").click());
$("replace").addEventListener("click", () => $("file").click());

async function loadExample() {
  try {
    const res = await fetch("example/proteinGroups_example.txt");
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const text = await res.text();
    loadText(text, "proteinGroups_example.txt", { baitName: "BAITX" });
    state.roles.set("BAIT", "bait"); state.roles.set("CTRL", "control"); state.suggested.clear();
    renderGroups(); compute();
  } catch (err) { setError(`Could not load the example: ${(err as Error).message}`); }
}
$("example").addEventListener("click", () => void loadExample());
$("hero-example").addEventListener("click", () => void loadExample());

/* ---------- Groups ---------- */

function selection() {
  const pick = (role: Role) => state.groups.filter((g) => state.roles.get(g.name) === role).flatMap((g) => g.columns.filter((c) => !state.excluded.has(c)));
  return { bait: pick("bait"), control: pick("control") };
}

function renderGroups() {
  const host = $("groups");
  const f = state.groupFilter.trim().toLowerCase();
  const filterBox = state.groups.length > 6 ? `<input type="search" class="group-filter" id="group-filter" placeholder="Filter ${state.groups.length} groups" aria-label="Filter groups" value="${esc(state.groupFilter)}" autocomplete="off" />` : "";
  host.innerHTML = filterBox + state.groups.map((g, i) => {
    const role = state.roles.get(g.name) ?? "off";
    if (f && !g.name.toLowerCase().includes(f) && role === "off") return "";
    const reps = g.columns.length > 1 && role !== "off"
      ? `<div class="reps-list">${g.columns.map((c, j) => `<label class="rep" title="${esc(c)}"><input type="checkbox" data-col="${esc(c)}" ${state.excluded.has(c) ? "" : "checked"}>${j + 1}</label>`).join("")}</div>` : "";
    const seg = (r: Role, label: string) => `<button type="button" data-g="${i}" data-role="${r}" aria-pressed="${role === r}">${label}</button>`;
    const tag = state.suggested.has(g.name) ? '<span class="tag-suggest" title="Guessed from the name. Please check it.">suggested</span>' : "";
    return `<div class="group ${role}"><div class="head"><span class="gname" title="${esc(g.name)}">${esc(g.name)}</span>${tag}<span class="reps">${g.columns.length} replicate${g.columns.length > 1 ? "s" : ""}</span></div><div class="seg" role="group" aria-label="Role of ${esc(g.name)}">${seg("bait", "Bait")}${seg("control", "Control")}${seg("off", "Off")}</div>${reps}</div>`;
  }).join("");
  host.querySelectorAll<HTMLButtonElement>("button[data-role]").forEach((b) => b.addEventListener("click", () => {
    const g = state.groups[Number(b.dataset.g)];
    state.roles.set(g.name, b.dataset.role as Role);
    state.suggested.delete(g.name);
    renderGroups(); schedule("compute");
  }));
  const gf = host.querySelector<HTMLInputElement>("#group-filter");
  gf?.addEventListener("input", () => { state.groupFilter = gf.value; const pos = gf.selectionStart ?? gf.value.length; renderGroups(); const again = $("group-filter") as HTMLInputElement; again.focus(); again.setSelectionRange(pos, pos); });
  host.querySelectorAll<HTMLInputElement>("input[data-col]").forEach((c) => c.addEventListener("change", () => {
    if (c.checked) state.excluded.delete(c.dataset.col!); else state.excluded.add(c.dataset.col!);
    schedule("compute");
  }));
}

/* ---------- Stoichiometry ---------- */

const peptideCache = new Map<number, Promise<ReturnType<typeof parsePeptideTable>>>();
const loadPeptides = (taxid: number) => {
  if (!peptideCache.has(taxid)) {
    peptideCache.set(taxid, fetch(`peptides/${taxid}.tsv`).then((r) => { if (!r.ok) throw new Error(`peptide table ${taxid}: HTTP ${r.status}`); return r.text(); }).then(parsePeptideTable));
  }
  return peptideCache.get(taxid)!;
};
let pendingOrg: string | undefined;
fetch("peptides/index.json").then((r) => r.json()).then((list: { taxid: number; name: string; proteins: number }[]) => {
  $("organism").innerHTML = list.map((o) => `<option value="${o.taxid}">${o.name} (${o.proteins.toLocaleString()} proteins)</option>`).join("");
  if (pendingOrg) ($("organism") as HTMLSelectElement).value = pendingOrg;
}).catch(() => { /* stoichiometry stays unavailable without the tables */ });
$("stoich").addEventListener("change", () => { $("stoichOpts").hidden = !($("stoich") as HTMLInputElement).checked; schedule("compute"); });
for (const id of ["organism", "stoichBait"]) $(id).addEventListener("change", () => schedule("compute"));

async function updateStoich() {
  state.stoich = new Map(); state.stoichError = "";
  if (!($("stoich") as HTMLInputElement).checked || !state.result) return;
  const key = ($("stoichBait") as HTMLInputElement).value.trim();
  if (!key) { state.stoichError = "Enter the bait protein to estimate stoichiometry."; void draw(); return; }
  try {
    const peptides = await loadPeptides(Number(($("organism") as HTMLSelectElement).value));
    state.stoich = computeStoichiometry(state.result.proteins, peptides, key);
  } catch (err) { state.stoichError = (err as Error).message; }
  void draw();
}

/* ---------- Compute ---------- */

function currentParams(): Params {
  const { bait, control } = selection();
  return {
    bait, control, welch: ($("test") as HTMLSelectElement).value === "welch",
    shift: state.num.shift, shrink: state.num.shrink, minFoldChange: state.num.minFoldChange, curvature: state.num.curvature, seed: defaultParams.seed,
    imputation: ($("imputation") as HTMLSelectElement).value as Params["imputation"],
  };
}

function showState(kind: "hero" | "prompt" | "results") {
  const wasResults = !$("results").hidden;
  $("tool-view").classList.toggle("has-results", kind === "results");
  // On a phone the settings sit below the plot, so fold the group list away when the first result appears.
  if (kind === "results" && !wasResults && window.matchMedia("(max-width: 960px)").matches) ($("s2") as HTMLDetailsElement).open = false;
  $("hero").hidden = kind !== "hero";
  $("prompt").hidden = kind !== "prompt";
  $("results").hidden = kind !== "results";
  if (kind !== "results") $("mbar").hidden = true;
}

function showPrompt(title: string, text: string, items: [string, boolean][]) {
  $("prompt-title").textContent = title;
  $("prompt-text").textContent = text;
  $("checklist").innerHTML = items.map(([t, ok]) => `<span class="${ok ? "ok" : ""}">${esc(t)}</span>`).join("");
  showState("prompt");
}

function compute() {
  if (!state.table) { showState("hero"); return; }
  const p = currentParams();
  setError("");
  const baitGroups = state.groups.filter((g) => state.roles.get(g.name) === "bait");
  const ctrlGroups = state.groups.filter((g) => state.roles.get(g.name) === "control");
  const nb = baitGroups.length, nc = ctrlGroups.length;
  const names = (gs: Group[]) => gs.map((g) => g.name).join(", ");
  const items: [string, boolean][] = [
    [nb ? `Bait: ${names(baitGroups)}${nb > 1 ? " (pooled)" : ""}` : "Bait: not chosen yet", nb > 0],
    [nc ? `Control: ${names(ctrlGroups)}${nc > 1 ? " (pooled)" : ""}` : "Control: not chosen yet", nc > 0],
  ];
  markStep("s2", false, nb && nc ? `${p.bait.length} vs ${p.control.length} replicates` : "Pick two groups");
  if (!nb || !nc) { state.result = null; showPrompt("Choose bait and control", "In step 2, mark one group as bait and one as control. The plot appears as soon as both are set.", items); return; }
  if (p.bait.length < 2 || p.control.length < 2) { state.result = null; showPrompt("Each side needs two replicates", "Tick more replicates, or pick a group that has at least two.", items); return; }
  try {
    state.result = computeVolcano(state.table, p);
    state.params = p;
    state.perm = null;
    if (($("mode") as HTMLSelectElement).value === "perm") {
      state.perm = permutationFdr(state.result.proteins, p.bait.length, state.num.s0, state.num.fdr, state.num.perms, p.seed);
    }
    markStep("s2", true, `${p.bait.length} vs ${p.control.length} replicates${nb > 1 || nc > 1 ? ", pooled" : ""}`);
    renderQc(state.result);
    void updateStoich();
    if (state.fitAxes) {
      // Fit the axes to the data once per file; the user can still adjust the limits afterwards.
      const xs = state.result.proteins.map((q) => Math.abs(q.x)), ys = state.result.proteins.map((q) => q.y);
      setters.xMax?.(Math.min(25, Math.max(5, Math.ceil(Math.max(...xs)) + 1)));
      setters.yMax?.(Math.min(25, Math.max(3, Math.ceil(Math.max(...ys)) + 1)));
      state.fitAxes = false;
    }
    void draw();
  } catch (err) { setError((err as Error).message); }
}

/* ---------- Draw ---------- */

type Box = [number, number, number, number];
const overlaps = (a: Box, b: Box) => a[0] < b[2] && a[2] > b[0] && a[1] < b[3] && a[3] > b[1];
const LABEL_POS = ["top right", "bottom right", "top left", "bottom left", "middle right", "middle left", "top center", "bottom center"];

/** Greedy label placement in pixel space. Returns the Plotly text position per point, or "" when no free spot exists. */
function makePlacer(W: number, H: number, xMax: number, yMax: number, fs: number, markers: { x: number; y: number }[]) {
  const toPx = (x: number, y: number): [number, number] => [((x + xMax) / (2 * xMax)) * W, (1 - y / yMax) * H];
  const dots: Box[] = markers.map((m) => { const [px, py] = toPx(m.x, m.y); return [px - 5, py - 5, px + 5, py + 5]; });
  const placed: Box[] = [];
  return (pts: { x: number; y: number; gene: string }[]): string[] => pts.map((p) => {
    const [px, py] = toPx(p.x, p.y);
    const w = p.gene.length * fs * 0.6 + 4, h = fs + 4, g = 8;
    for (const pos of LABEL_POS) {
      const [v, hz] = pos.split(" ");
      const x0 = hz === "right" ? px + g : hz === "left" ? px - g - w : px - w / 2;
      const y0 = v === "top" ? py - g - h : v === "bottom" ? py + g : py - h / 2;
      const box: Box = [x0, y0, x0 + w, y0 + h];
      if (box[2] > W + 36 || box[0] < -36 || box[1] < -10 || box[3] > H + 10) continue;
      if (placed.some((b) => overlaps(b, box)) || dots.some((d) => overlaps(d, box))) continue;
      placed.push(box);
      return pos;
    }
    return "";
  });
}

const MODE_LABEL: Record<string, string> = { hyperbola: "Hyperbolic curve", perm: "Permutation FDR", bh: "Benjamini-Hochberg" };

function modeSummary(mode: string): string {
  const imp = ($("imputation") as HTMLSelectElement).value;
  const base = mode === "perm" ? `s0 ${state.num.s0}, FDR ${state.num.fdr} (experimental)` : mode === "bh" ? `BH q < ${state.num.qcut} (experimental)` : `min enrichment ${state.num.minFoldChange}, curvature ${state.num.curvature}`;
  const impText = imp === "normal" ? `fill-in shift ${state.num.shift}, spread ${state.num.shrink}` : `${imp} imputation (experimental)`;
  return `msVolcano ${__APP_VERSION__} | ${($("test") as HTMLSelectElement).value === "welch" ? "Welch" : "Student"} t-test | ${base} | ${impText}`;
}

async function draw(forceLight = false) {
  const r = state.result, table = state.table;
  if (!r || !table || !state.params) return;
  showState("results");
  const p = { ...state.params, minFoldChange: state.num.minFoldChange, curvature: state.num.curvature };
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
  state.hits = hits;
  const rest = prot.filter((q) => !q.significant);
  const { xMax, yMax, labelSize } = state.num;
  const narrow = window.innerWidth < 600;
  const curve = cutoffCurve(p.minFoldChange, p.curvature, xMax, yMax);
  const manual = ($("manual") as HTMLInputElement).value.split(";").map((s) => s.trim().toLowerCase()).filter(Boolean);
  const manualHits = manual.length ? prot.filter((q) => manual.some((m) => q.gene.toLowerCase().startsWith(m))) : [];
  const labelLimit = Math.min(MAX_LABELS, narrow ? Math.min(5, state.num.labelCount) : state.num.labelCount);
  const labelled = hits.filter((q) => !manualHits.includes(q)).slice(0, labelLimit);

  const c = forceLight
    ? LIGHT_PLOT
    : { ink: cssVar("--plot-ink"), grid: cssVar("--plot-grid"), axis: cssVar("--plot-axis"), bg: cssVar("--plot-bg"), point: cssVar("--plot-point"), hit: cssVar("--plot-hit"), pick: cssVar("--plot-pick"), curve: cssVar("--plot-curve") };
  const marginB = narrow ? 120 : 96, marginT = ($("title") as HTMLInputElement).value ? 56 : 20;
  const plotH = Math.max(200, $("plot").clientHeight - marginT - marginB);
  const hitSize = (q: { row: number }) => (state.stoich.size ? Math.min(30, 5 + 3 * Math.log2((state.stoich.get(q.row) ?? 0) + 1.5)) : 8);
  const placer = makePlacer($("plot").clientWidth - 60 - (narrow ? 40 : 28), plotH, xMax, yMax, labelSize, [...hits, ...manualHits]);
  // Highlighted proteins claim their label spot first, then the strongest interactors.
  const manualPos = placer(manualHits);
  const labelPos = placer(labelled);
  const hiddenLabels = labelPos.filter((x) => !x).length;
  const trace = (pts: typeof prot, color: string, size: number | number[], text: boolean, opacity: number, symbol = "circle", pos: string[] = []) => ({
    type: "scattergl", mode: text ? "markers+text" : "markers",
    x: pts.map((q) => q.x), y: pts.map((q) => q.y),
    text: pts.map((q, i) => (pos[i] ? q.gene : "")), textposition: pos.map((x) => x || "top right"), textfont: { size: labelSize, color: c.ink }, cliponaxis: false,
    customdata: pts.map((q) => [q.id, q.gene]), hovertemplate: "<b>%{customdata[1]}</b><br>%{customdata[0]}<br>difference %{x:.2f}<br>-log10 p %{y:.2f}<extra></extra>",
    marker: { color, size, opacity, symbol, line: { color: c.bg, width: symbol === "circle" && color === c.point ? 0 : 1.2 } },
  });
  const line = (pts: { x: number; y: number }[]) => ({ type: "scatter", mode: "lines", x: pts.map((q) => q.x), y: pts.map((q) => q.y), line: { color: c.curve, dash: "dash", width: 1.3 }, hoverinfo: "skip", showlegend: false });

  const title = ($("title") as HTMLInputElement).value;
  const baitName = ($("bait-name") as HTMLInputElement).value;
  const axisFont = { color: c.ink, size: 13 };
  const footGap = narrow ? 56 : 62; // distance below the x axis, past the tick labels and axis title
  const layout = {
    font: { color: c.ink, family: "Geist Variable, system-ui, sans-serif", size: 13 },
    title: title ? { text: esc(title), x: 0.5, font: { size: 16 } } : undefined,
    xaxis: { title: { text: narrow ? "Δ mean log2 LFQ" : "Difference in mean log2 LFQ (bait minus control)", standoff: 10, font: axisFont }, range: [-xMax, xMax], zeroline: false, gridcolor: c.grid, showline: true, linecolor: c.axis, ticks: "outside", tickcolor: c.axis, tickfont: { color: c.ink } },
    yaxis: { title: { text: "-log10 p value", standoff: 8, font: axisFont }, range: [0, yMax], zeroline: false, gridcolor: c.grid, showline: true, linecolor: c.axis, ticks: "outside", tickcolor: c.axis, tickfont: { color: c.ink } },
    showlegend: false, hovermode: "closest", margin: { t: marginT, r: narrow ? 40 : 28, b: marginB, l: 60 },
    paper_bgcolor: c.bg, plot_bgcolor: c.bg,
    hoverlabel: { bgcolor: c.bg, bordercolor: c.grid, font: { color: c.ink } },
    annotations: [{ xref: "paper", yref: "paper", x: 0, y: -(footGap / plotH), xanchor: "left", yanchor: "top", showarrow: false, text: esc(`${baitName ? `Bait: ${baitName}. ` : ""}${modeSummary(mode)}`).replace(/ \| /g, narrow ? "<br>" : " | "), font: { size: 10, color: c.axis } }],
  };
  const plotEl = $("plot");
  plotEl.setAttribute("aria-label", `Volcano plot of ${prot.length} proteins with ${hits.length} significant interactors${hits.length ? `, led by ${hits.slice(0, 3).map((h) => h.gene).join(", ")}` : ""}. The full list is in the table below.`);
  await Plotly.react(plotEl, [
    trace(rest, c.point, 6, false, 0.6),
    trace(hits, c.hit, hits.map(hitSize), false, 0.95),
    trace(labelled, c.hit, labelled.map(hitSize), true, 0.95, "circle", labelPos),
    trace(manualHits, c.pick, 12, true, 1, "diamond", manualPos),
    ...(mode === "hyperbola" ? [line(curve.right), line(curve.left)] : []),
  ], layout, { responsive: true, displaylogo: false, modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d", "toImage"] });
  if (forceLight) return;

  // Summary
  const d = r.dropped;
  const removed = d.contaminantOrReverse + d.absentInBait + d.untestable;
  $("stats").innerHTML = `
    <div class="stat hit"><span>Significant interactors</span><b>${hits.length.toLocaleString()}</b><small>${prot.length ? ((hits.length / prot.length) * 100).toFixed(1) : "0"}% of ${prot.length.toLocaleString()} proteins tested</small></div>
    <div class="stat"><span>Proteins tested</span><b>${prot.length.toLocaleString()}</b><small>from ${table.rows.length.toLocaleString()} rows</small></div>
    <div class="stat"><span>Removed</span><b>${removed.toLocaleString()}</b><small>${d.contaminantOrReverse} contaminant or reverse, ${d.absentInBait} absent in bait${d.untestable ? `, ${d.untestable} untestable` : ""}</small></div>`;
  $("plot-title").textContent = baitName ? `Volcano plot: ${baitName}` : "Volcano plot";
  markStep("s3", false, `${MODE_LABEL[mode]} · ${hits.length} hits`);
  $("legend-cutoff").hidden = mode !== "hyperbola";
  $("mbar-count").textContent = `${hits.length} interactors`;
  $("mbar").hidden = false;
  const notes: string[] = [];
  if (r.logTransformed) notes.push("LFQ values looked unlogged and were log2 transformed.");
  if (mode === "perm") notes.push(state.perm ? `Experimental: permutation FDR ${state.perm.fdr}, ${state.perm.permutations} permutations, score threshold ${state.perm.threshold.toFixed(2)}. Proteins can be called inside the funnel because each one is scored on its own variability, not by the fixed curve.` : "Experimental: not enough replicates for permutation FDR.");
  if (mode === "bh") notes.push(`Experimental: Benjamini-Hochberg q below ${state.num.qcut}, enriched side only.`);
  if (($("imputation") as HTMLSelectElement).value !== "normal") notes.push("Experimental imputation method in use.");
  if (hits.length > labelLimit) notes.push(`Showing labels for the top ${labelLimit} of ${hits.length} interactors.${narrow ? " On phones labels are limited to 5." : ' Change "Labels shown" in step 4, or use the table.'}`);
  if (hiddenLabels) notes.push(`${hiddenLabels} label${hiddenLabels > 1 ? "s were" : " was"} left off to avoid overlaps. Hover a point to see its name, or use the table.`);
  const pooledB = state.groups.filter((g) => state.roles.get(g.name) === "bait").length, pooledC = state.groups.filter((g) => state.roles.get(g.name) === "control").length;
  if (pooledB > 1 || pooledC > 1) notes.push(`Pooled groups: ${pooledB > 1 ? `${pooledB} bait groups` : ""}${pooledB > 1 && pooledC > 1 ? " and " : ""}${pooledC > 1 ? `${pooledC} control groups` : ""} are treated as one. Check that this is what you want in step 2.`);
  if (state.stoichError) notes.push(`Stoichiometry: ${state.stoichError}`);
  $("notes").innerHTML = notes.map((n) => `<div>${esc(n)}</div>`).join("");
  renderHits();
}

/* ---------- Tables ---------- */

const fmtSt = (v: number | undefined) => (v === undefined || !Number.isFinite(v) ? "" : v.toPrecision(3));

function renderHits() {
  const hits = state.hits;
  $("hit-count").textContent = `(${hits.length})`;
  const mode = ($("mode") as HTMLSelectElement).value;
  const showQ = mode === "bh";
  const showSt = state.stoich.size > 0;
  const f = state.filter.trim().toLowerCase();
  const rows = hits.filter((q) => !f || q.gene.toLowerCase().includes(f) || q.id.toLowerCase().includes(f));
  const key = state.sort.key;
  const val = (q: (typeof hits)[number]): number | string => key === "gene" ? q.gene.toLowerCase() : key === "id" ? q.id : key === "y" ? q.y : key === "q" ? (q.q ?? 0) : key === "st" ? (state.stoich.get(q.row) ?? -Infinity) : q.x;
  rows.sort((a, b) => { const A = val(a), B = val(b); return (A < B ? -1 : A > B ? 1 : 0) * state.sort.dir; });
  const th = (k: string, label: string, tip: string, num = false) => `<th class="${num ? "num" : ""}" title="${esc(tip)}" aria-sort="${key === k ? (state.sort.dir === 1 ? "ascending" : "descending") : "none"}"><button type="button" data-sort="${k}">${label}${key === k ? `<span aria-hidden="true">${state.sort.dir === 1 ? " ↑" : " ↓"}</span>` : ""}</button></th>`;
  if (!hits.length) {
    const advice = mode === "hyperbola" ? "Lower the minimum enrichment or curvature in step 3." : mode === "perm" ? "Raise the FDR or lower s0 in Experimental methods." : "Raise the q value cutoff in Experimental methods.";
    $("hits").innerHTML = `<div class="callout" style="margin:14px;border:0">${ICON.info}<span>No proteins pass the cutoff. ${advice}</span></div>`; return;
  }
  $("hits").innerHTML = `<table><thead><tr>${th("gene", "Gene", "Gene name")}${th("id", "Protein ID", "First majority protein ID")}${th("x", "log2 difference", "Mean log2 LFQ in bait minus control", true)}${th("y", "-log10 p", "Higher means stronger evidence", true)}${showQ ? th("q", "q value", "Benjamini-Hochberg adjusted p value", true) : ""}${showSt ? th("st", "Stoichiometry", "Relative to the bait, which is 1", true) : ""}</tr></thead><tbody>${rows.slice(0, state.rowsShown).map((q) => `<tr><td class="gene">${esc(q.gene)}</td><td>${esc(q.id)}</td><td class="num">${q.x.toFixed(2)}</td><td class="num">${q.y.toFixed(2)}</td>${showQ ? `<td class="num">${(q.q ?? 0).toExponential(1)}</td>` : ""}${showSt ? `<td class="num">${fmtSt(state.stoich.get(q.row))}</td>` : ""}</tr>`).join("")}</tbody></table>${rows.length > state.rowsShown ? `<div class="notes" style="display:flex;gap:12px;align-items:center">Showing ${state.rowsShown} of ${rows.length}. <button class="btn small" id="more" type="button">Show 50 more</button></div>` : ""}`;
  $("more")?.addEventListener("click", () => { state.rowsShown += 50; renderHits(); });
  $("hits").querySelectorAll<HTMLButtonElement>("button[data-sort]").forEach((b) => b.addEventListener("click", () => {
    const k = b.dataset.sort!;
    state.sort = { key: k, dir: state.sort.key === k ? (state.sort.dir === 1 ? -1 : 1) : k === "gene" || k === "id" ? 1 : -1 };
    renderHits();
  }));
}
$("hit-filter").addEventListener("input", (e) => { state.filter = (e.target as HTMLInputElement).value; state.rowsShown = 25; renderHits(); });

function renderQc(r: VolcanoResult) {
  const p = state.params!;
  const shorten = (n: string) => { const s = n.replace(/^LFQ intensity /i, ""); return s.length > 18 ? `${s.slice(0, 8)}...${s.slice(-7)}` : s; };
  const names = [...p.bait, ...p.control].map(shorten);
  const c = replicateCorrelation(r.proteins);
  const nb = p.bait.length;
  // Lowest correlation between two replicates of the same group.
  let worst = 1, worstPair = "";
  const scan = (lo: number, hi: number) => { for (let i = lo; i < hi; i++) for (let j = i + 1; j < hi; j++) if (c[i][j] < worst) { worst = c[i][j]; worstPair = `${names[i]} and ${names[j]}`; } };
  scan(0, nb); scan(nb, names.length);
  const bad = worst < 0.8;
  const verdict = bad
    ? `Check your replicates: ${worstPair} agree only moderately (r = ${worst.toFixed(2)}). One may be an outlier. You can untick it in step 2.`
    : `Replicates agree well. The lowest correlation within a group is ${worst.toFixed(2)}.`;
  // Scale the tint from the lowest value in the matrix to 1 so small differences between replicates stay visible.
  const minOff = Math.min(...c.flatMap((row, i) => row.filter((_, j) => j !== i)));
  const cell = (v: number) => `style="background:color-mix(in oklch, var(--color-accent) ${Math.round(Math.max(0, Math.min(1, (v - minOff) / Math.max(1e-6, 1 - minOff))) * 40)}%, transparent)"`;
  $("qcTable").innerHTML = `<div class="qc-verdict ${bad ? "warn" : ""}" style="margin:-10px -14px 10px">${bad ? ICON.alert.replace("<svg", '<svg width="18" height="18"') : ""}<span>${esc(verdict)}</span></div><table class="qc"><thead><tr><th></th>${names.map((n) => `<th scope="col">${esc(n)}</th>`).join("")}</tr></thead><tbody>${c.map((row, i) => `<tr><th scope="row">${esc(names[i])}</th>${row.map((v) => `<td ${cell(v)}>${v.toFixed(2)}</td>`).join("")}</tr>`).join("")}</tbody></table><p style="font-size:.82rem;color:var(--color-ink-2);margin-top:10px">Pearson correlation of imputed log2 LFQ values. Replicates of one condition should agree with each other more than with the other group. A deeper tint means a higher correlation, scaled from the lowest value shown to 1.</p>`;
  $("qc-dot").innerHTML = bad ? '<span class="dot-warn" role="img" aria-label="Replicate warning"></span>' : "";
}

/* ---------- Tabs ---------- */

function selectTab(which: "hits" | "qc", focus = false) {
  const hitsTab = which === "hits";
  $("tab-hits").setAttribute("aria-selected", String(hitsTab)); $("tab-hits").tabIndex = hitsTab ? 0 : -1;
  $("tab-qc").setAttribute("aria-selected", String(!hitsTab)); $("tab-qc").tabIndex = hitsTab ? -1 : 0;
  $("panel-hits").hidden = !hitsTab; $("panel-qc").hidden = hitsTab;
  if (focus) $(hitsTab ? "tab-hits" : "tab-qc").focus();
}
$("tab-hits").addEventListener("click", () => selectTab("hits"));
$("tab-qc").addEventListener("click", () => selectTab("qc"));
document.querySelector(".subtabs")!.addEventListener("keydown", (e) => {
  const k = (e as KeyboardEvent).key;
  if (k === "ArrowRight" || k === "End") selectTab("qc", true);
  else if (k === "ArrowLeft" || k === "Home") selectTab("hits", true);
});

/* ---------- Export ---------- */

const closeMenu = () => (($("exportMenu") as HTMLDetailsElement).open = false);
document.addEventListener("click", (e) => { if (!$("exportMenu").contains(e.target as Node)) closeMenu(); });
document.addEventListener("keydown", (e) => { if (e.key === "Escape") closeMenu(); });

async function exportImage(format: "png" | "svg") {
  closeMenu();
  await draw(true);
  try { await Plotly.downloadImage($("plot"), { format, filename: stem(), width: 1400, height: 1000 }); } finally { await draw(); }
}
$("png").addEventListener("click", () => void exportImage("png"));
$("svg").addEventListener("click", () => void exportImage("svg"));
$("csv").addEventListener("click", () => {
  closeMenu();
  if (!state.table || !state.params) return;
  const url = URL.createObjectURL(new Blob([toCsv(state.table, state.hits, state.params, state.stoich)], { type: "text/csv" }));
  const a = Object.assign(document.createElement("a"), { href: url, download: `${stem()}_hits.csv` });
  a.click(); URL.revokeObjectURL(url);
});
async function copyGenes() {
  closeMenu();
  const text = state.hits.map((q) => q.gene).join("\n");
  try { await navigator.clipboard.writeText(text); toast(`Copied ${state.hits.length} gene names`); } catch { toast("Copy is blocked by the browser"); }
}
$("copy").addEventListener("click", () => void copyGenes());
$("copy2").addEventListener("click", () => void copyGenes());
$("mbar-adjust").addEventListener("click", () => { const st = $("s3") as HTMLDetailsElement; st.open = true; st.scrollIntoView({ behavior: "smooth", block: "start" }); });
$("mbar-plot").addEventListener("click", () => $("results").scrollIntoView({ behavior: "smooth", block: "start" }));

/* ---------- Reset and share ---------- */

const selVal = (id: string) => ($(id) as HTMLSelectElement | HTMLInputElement).value;

function encodeSettings(): string {
  const saved: Saved = {
    n: { ...state.num }, r: Object.fromEntries([...state.roles].filter(([, v]) => v !== "off")), x: [...state.excluded],
    m: selVal("mode"), t: selVal("test"), i: selVal("imputation"), man: selVal("manual"), ti: selVal("title"), b: selVal("bait-name"),
    so: ($("stoich") as HTMLInputElement).checked, sg: selVal("organism"), sb: selVal("stoichBait"),
  };
  return btoa(unescape(encodeURIComponent(JSON.stringify(saved)))).replace(/\+/g, "-").replace(/\//g, "_").replace(/=+$/, "");
}
function decodeSettings(raw: string): Saved | null {
  try { return JSON.parse(decodeURIComponent(escape(atob(raw.replace(/-/g, "+").replace(/_/g, "/"))))) as Saved; } catch { return null; }
}
function applySavedControls(s: Saved) {
  for (const [k, v] of Object.entries(s.n ?? {})) if (k in state.num && Number.isFinite(v)) setters[k as keyof typeof state.num]?.(v);
  const set = (id: string, v: string | undefined) => { if (v !== undefined) ($(id) as HTMLSelectElement | HTMLInputElement).value = v; };
  set("mode", s.m); set("test", s.t); set("imputation", s.i); set("manual", s.man); set("title", s.ti); set("bait-name", s.b); set("stoichBait", s.sb);
  if (s.so) { ($("stoich") as HTMLInputElement).checked = true; $("stoichOpts").hidden = false; }
  if (s.sg) { pendingOrg = s.sg; if ($("organism").children.length) ($("organism") as HTMLSelectElement).value = s.sg; }
}
function applyPendingRoles() {
  const s = state.pending;
  if (!s) return;
  let matched = 0;
  for (const g of state.groups) if (s.r && s.r[g.name]) { state.roles.set(g.name, s.r[g.name]); state.suggested.delete(g.name); matched++; }
  for (const c of s.x ?? []) state.excluded.add(c);
  if (matched) toast("Your saved settings were applied.");
  state.pending = null;
}
$("share").addEventListener("click", async () => {
  const url = `${location.origin}${location.pathname}#s=${encodeSettings()}`;
  try { await navigator.clipboard.writeText(url); toast("Settings link copied. Open it and add the same file to restore this view."); } catch { toast("Copy is blocked by the browser"); }
});
$("reset").addEventListener("click", () => {
  for (const [k, v] of Object.entries(DEFAULT_NUM)) setters[k as keyof typeof state.num]?.(v);
  for (const id of ["manual", "title", "stoichBait"]) ($(id) as HTMLInputElement).value = "";
  ($("mode") as HTMLSelectElement).value = "hyperbola"; ($("test") as HTMLSelectElement).value = "student"; ($("imputation") as HTMLSelectElement).value = "normal";
  ($("stoich") as HTMLInputElement).checked = false; $("stoichOpts").hidden = true;
  state.fitAxes = true;
  compute();
  toast("Cutoff, plot and method settings reset. Your bait and control choices were kept.");
});
{
  const m = /^#s=(.+)$/.exec(location.hash);
  const saved = m ? decodeSettings(m[1]) : null;
  if (saved) {
    state.pending = saved;
    applySavedControls(saved);
    toast("Settings restored. Add your file to continue.");
  }
}
