import "@fontsource-variable/geist";
import "@fontsource-variable/geist-mono";
import Plotly from "plotly.js-cartesian-dist-min";
import "./style.css";
import { aboutHtml } from "./about";
import { benjaminiHochberg, permutationFdr, type PermutationResult } from "./engine/experimental";
import { baitRecovery, coverageCheck, loadingCheck, observedCorrelation, samplePca, type Check, type PcaResult } from "./engine/qc";
import { parseDelimited, sniffSeparator, type Table } from "./engine/parse";
import { computeStoichiometry, parsePeptideTable } from "./engine/stoich";
import { detectFormat, explainUnsupported } from "./engine/formats";
import { computeVolcano, cutoffCurve, defaultParams, lfqColumns, presenceCsv, shiftSensitivity, toCsv, type Params, type Protein, type VolcanoResult } from "./engine/volcano";
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
      <label class="drop" id="drop" for="file">${ICON.upload}<strong>Drop your protein table here</strong><span>MaxQuant, FragPipe, DIA-NN and more. Click to browse.</span></label>
      <input id="file" class="sr" type="file" accept=".txt,.tsv,.csv" />
      <div class="filechip" id="filechip" hidden></div>
      <p class="help" id="fmt-note" style="font-size:.78rem;color:var(--color-ink-3)" hidden></p>
      <div class="field" id="fmt-force-wrap" hidden><label for="fmt-force">Not read correctly?</label><select id="fmt-force"><option value="auto">Detect the format automatically</option><option value="generic">Treat as a plain table (numeric columns are samples)</option></select></div>
      <div class="row-actions"><button class="btn small" id="example" type="button">${ICON.flask}Try an example dataset (simulated)</button><button class="btn small" id="replace" type="button" hidden>Replace file</button></div>
      <details style="font-size:.82rem;color:var(--color-ink-2)"><summary style="cursor:pointer;font-weight:600">What file do I need?</summary><p style="margin-top:6px">A protein table with one intensity column per sample. Supported: MaxQuant <code>proteinGroups.txt</code> (with <b>LFQ</b> switched on), FragPipe <code>combined_protein.tsv</code>, DIA-NN <code>report.pg_matrix.tsv</code>, Spectronaut protein group pivot, Proteome Discoverer protein export, Perseus matrices, or any table of numbers. DIA-NN, Spectronaut and Proteome Discoverer readers are beta.</p></details>
      <p class="help" style="font-size:.78rem;color:var(--color-ink-3);display:flex;gap:6px;align-items:center"><span style="width:14px;height:14px;display:inline-flex">${ICON.lock}</span>Stays on this computer. Works offline once loaded.</p>`, { open: true })}
    ${step("s2", "2", "Choose bait and control", "Pick two groups", `
      <p class="help" style="font-size:.8rem;color:var(--color-ink-3)">Replicates are grouped by name. Mark one group as bait and one as control. Untick a replicate to leave it out. Several bait groups are pooled into one.</p>
      <div class="groups" id="groups"></div>`)}
    ${step("s3", "3", "Set the cutoff", "Hyperbolic curve", `
      <div class="field"><label for="test">Statistical test</label><select id="test"><option value="student">Student t-test (equal variance)</option><option value="welch">Welch t-test (unequal variance)</option></select></div>
      <div id="sl-cutoff"></div>
      <details style="margin-top:4px"><summary style="cursor:pointer;font-weight:600;font-size:.85rem">Missing values and bait-only proteins <span class="pill new">New in v2</span></summary>
        <div style="display:grid;gap:var(--space-3);margin-top:var(--space-3)">
          <div class="field"><label for="presence">Seen in bait, never in a control</label><select id="presence"><option value="separate">List separately (recommended)</option><option value="include">Fill in and test (as published)</option></select><span class="help">A protein found in at least 2 bait replicates and not detected in any control cannot be ranked fairly by a t-test, which would use filled-in control values. It gets its own list with replicate counts. Proteins measured once in total are left out. Pick "Fill in and test" to match the 2016 tool.</span></div>
          <div class="field"><label for="imputation">Filling in other missing values</label><select id="imputation"><option value="normal">Shifted normal (as published)</option><option value="mindet">MinDet: low quantile of each column (experimental)</option><option value="minprob">MinProb: random draw around low quantile (experimental)</option><option value="none">Do not fill in: test the measured values only</option></select><span class="help">Points that include a filled-in value are hollow. "Do not fill in" needs two measured values per group, so proteins with fewer drop out of the volcano.</span></div>
          <div class="field"><label for="logmode">Intensity scale</label><select id="logmode"><option value="auto">Detect automatically</option><option value="raw">Raw intensities (not logged)</option><option value="log2">Already log2</option></select><span class="help">msVolcano works in log2. It reads the note under the plot to tell you what it decided. Set this if the guess is wrong, for example for log10 data or spectral counts.</span></div>
          <div id="sl-impute" style="display:grid;gap:var(--space-3)"></div>
        </div></details>`)}
    ${step("s4", "4", "Plot and labels", "Axes, names, titles", `
      <div id="sl-plot"></div>
      <div class="field"><label for="manual">Highlight proteins</label><input id="manual" type="text" placeholder="Wdr5;Mll2" autocomplete="off" /><span class="help">Gene names separated by semicolons. Matches names that start with your text.</span></div>
      <div class="field"><label for="title">Plot title</label><input id="title" type="text" autocomplete="off" /></div>
      <div class="field"><label for="bait-name">Bait name (gene)</label><input id="bait-name" type="text" autocomplete="off" /><span class="help">Shown in the title. Also used to check that the bait is enriched.</span></div>`)}
    ${step("s5", "+", "Stoichiometry", "", `
      <label class="check"><input id="stoich" type="checkbox" /> <span>Estimate abundance relative to the bait</span></label>
      <div id="stoichOpts" hidden style="display:grid;gap:12px">
        <div class="field"><label for="organism">Organism</label><select id="organism"></select></div>
        <div class="field"><label for="stoichBait">Bait protein</label><input id="stoichBait" type="text" placeholder="Defaults to the bait name above" autocomplete="off" /></div>
        <p class="help" style="font-size:.78rem;color:var(--color-ink-3)">iBAQ style, as in the 2016 paper. Peptide counts come from UniProt Swiss-Prot and follow the standard rule of no cleavage before proline, so values can differ from the original tool.</p>
      </div>`, { optional: true })}
    ${step("s6", "+", "Experimental methods", "", `
      <p class="help" style="font-size:.78rem;color:var(--color-ink-3)">New in version 2, not part of the 2016 publication and not peer reviewed. Compare with an established tool before relying on them.</p>
      <div class="field"><label for="mode">Cutoff method</label><select id="mode"><option value="hyperbola">Hyperbolic curve (as published)</option><option value="perm">Permutation FDR (Perseus style)</option><option value="bh">Benjamini-Hochberg q value</option></select></div>
      <div id="sl-exp"></div>
      <p class="help" style="font-size:.78rem;color:var(--color-ink-3)">Permutation FDR assumes samples are interchangeable. With few replicates there are only a handful of label shuffles (4 vs 4 gives 69), and it can call more proteins than the published curve when many background proteins differ slightly between groups.</p>`, { exp: true })}
    <div class="sidebar-foot"><button class="link-btn" id="reset" type="button">Reset cutoff and plot settings</button><button class="link-btn" id="share" type="button">Copy settings link</button></div>
  </aside>

  <main class="main">
    <section class="hero" id="hero">
      <div style="display:grid;gap:var(--space-4)">
        <h2>See what binds your bait</h2>
        <p class="lede">Turn a label free interactomics run (MaxQuant, FragPipe, DIA-NN and more) into a volcano plot and a ranked list of interactors. Choose bait and control, tune the cutoff, check quality, export the figure.</p>
        <div class="cta">
          <button class="btn primary" id="hero-browse" type="button">${ICON.upload}Choose a file</button>
          <button class="btn" id="hero-example" type="button">${ICON.flask}Try an example dataset</button>
        </div>
        <p class="privacy">${ICON.lock}<span>Everything runs in your browser. Your data is never uploaded. You can also drop a file anywhere on this page.</span></p>
      </div>
      ${heroArt()}
      <div class="how">
        <div><b>Add your protein table</b><span>Reverse hits and contaminants are removed for you.</span></div>
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
              <button id="csvp" type="button">${ICON.download}Present-only list as CSV</button>
              <button id="copy" type="button">${ICON.copy}Copy interactor genes</button>
              <button id="copyall" type="button">${ICON.copy}Copy interactors and present-only genes</button>
            </div>
          </details>
        </header>
        <div id="plot" role="group" aria-label="Interactive volcano plot. Hover points to see gene names. The interactors are listed in the table below."></div>
        <div class="legend-key" id="legend"><span><i style="background:var(--plot-hit)"></i>Significant</span><span><i style="background:var(--plot-point)"></i>Other proteins</span><span><i class="h"></i><i class="h g"></i>Hollow: includes a filled-in value</span><button class="link-btn" id="legend-presence" type="button" hidden></button><span><i class="d" style="background:var(--plot-pick)"></i>Highlighted</span><span id="legend-cutoff"><i class="l"></i>Cutoff</span></div>
        <div class="notes" id="notes"></div>
      </div>
      <div class="card" style="margin-top:var(--space-4)">
        <div class="subtabs" role="tablist" aria-label="Results views">
          <button role="tab" id="tab-hits" aria-selected="true" aria-controls="panel-hits" tabindex="0" type="button">Interactors <span id="hit-count" class="num"></span></button>
          <button role="tab" id="tab-presence" aria-selected="false" aria-controls="panel-presence" tabindex="-1" type="button">Present only in bait <span id="presence-count" class="num"></span></button>
          <button role="tab" id="tab-qc" aria-selected="false" aria-controls="panel-qc" tabindex="-1" type="button">Quality checks <span id="qc-dot"></span></button>
        </div>
        <div id="panel-hits" role="tabpanel" aria-labelledby="tab-hits">
          <div class="toolbar">
            <input id="hit-filter" type="search" placeholder="Filter by gene or ID" aria-label="Filter interactors" autocomplete="off" />
            <span class="spacer" style="flex:1"></span>
            <button class="btn small" id="copy2" type="button">${ICON.copy}Copy genes</button>
          </div>
          <div class="tablewrap" id="hits" tabindex="0" role="region" aria-label="Interactors table, scrollable"></div>
        </div>
        <div id="panel-presence" role="tabpanel" aria-labelledby="tab-presence" hidden>
          <div class="callout info" style="margin:14px;border-radius:var(--radius-m)">${ICON.info}<span><b>Seen in bait, never in a control.</b> Quantified in at least 2 bait replicates and not detected in any control replicate, so no p value is computed. A protein can also be missing from controls by chance. Ranked by replicate count, then signal. Often strong candidates: check the counts and known contaminants. These are not part of the interactors list.</span></div>
          <div class="toolbar"><span class="spacer" style="flex:1"></span><button class="btn small" id="copy3" type="button">${ICON.copy}Copy genes</button><button class="btn small" id="csvp2" type="button">${ICON.download}CSV</button></div>
          <div class="tablewrap" id="presenceTable" tabindex="0" role="region" aria-label="Present only in bait table, scrollable"></div>
        </div>
        <div id="panel-qc" role="tabpanel" aria-labelledby="tab-qc" hidden>
          <div id="qcBody"></div>
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

interface Saved { pr?: string; lm?: string; n?: Record<string, number>; r?: Record<string, Role>; x?: string[]; m?: string; t?: string; i?: string; man?: string; ti?: string; b?: string; so?: boolean; sg?: string; sb?: string }

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
  fmt: null as { label: string; note?: string } | null,
  clean: ((c: string) => c) as (c: string) => string,
  rowsShown: 25,
  pca: null as PcaResult | null,
  baseChecks: [] as Check[],
  pcaNames: [] as string[],
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
for (const id of ["test", "imputation", "mode", "presence", "logmode"]) $(id).addEventListener("input", () => schedule("compute"));
for (const id of ["manual", "title", "bait-name"]) $(id).addEventListener("input", () => schedule("draw"));

/* ---------- Data loading ---------- */

const paint = () => new Promise<void>((r) => requestAnimationFrame(() => setTimeout(r, 0)));

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
    const why = explainUnsupported(table);
    fail(`${why ? why + " " : ""}Found ${lfq.length} sample columns. msVolcano needs at least 4 (two bait and two control replicates). It reads MaxQuant, FragPipe, DIA-NN, Spectronaut and Proteome Discoverer tables, Perseus matrices, and plain tables of numbers.`);
    return;
  }
  const detected = detectFormat(table)!;
  const groups = groupColumns(lfq, detected.clean);
  if (groups.length < 2) {
    fail(`Only one group of LFQ columns was found (${groups[0].name}). You need at least two: one bait and one control.`);
    return;
  }
  state.table = table; state.fileName = name; state.clean = detected.clean; state.fmt = { label: detected.label, note: detected.note }; state.fitAxes = true; state.result = null;
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
  const fmt = detectFormat(table);
  const note = $("fmt-note");
  note.hidden = !fmt;
  note.textContent = fmt ? `Detected: ${fmt.label}.${fmt.note ? " " + fmt.note : ""}` : "";
  ($("bait-name") as HTMLInputElement).value = preset?.baitName ?? ($("bait-name") as HTMLInputElement).value;
  $("replace").hidden = false;
  $("tool-view").classList.add("has-data");
  markStep("s1", true, detected.short);
  $("fmt-force-wrap").hidden = false;
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
$("fmt-force").addEventListener("change", () => {
  if (!state.table) return;
  state.table.forceFormat = ($("fmt-force") as HTMLSelectElement).value === "generic" ? "generic" : undefined;
  const text = state.table; // keep the parsed table, only re-detect the columns
  const detected = detectFormat(text);
  if (!detected || detected.sampleColumns.length < 4) { toast("That reading did not find four sample columns."); state.table.forceFormat = undefined; ($("fmt-force") as HTMLSelectElement).value = "auto"; return; }
  const groups = groupColumns(detected.sampleColumns, detected.clean);
  state.groups = groups; state.clean = detected.clean; state.fmt = { label: detected.label, note: detected.note };
  const seen = { bait: false, control: false };
  state.roles = new Map(groups.map((g) => { let r = suggestRole(g.name); if (r !== "off") { if (seen[r]) r = "off"; else seen[r] = true; } return [g.name, r]; }));
  state.excluded = new Set(); state.suggested = new Set(groups.filter((g) => state.roles.get(g.name) !== "off").map((g) => g.name));
  $("fmt-note").textContent = `Detected: ${detected.label}.${detected.note ? " " + detected.note : ""}`; markStep("s1", true, detected.short);
  renderGroups(); compute();
});
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
    presence: ($("presence") as HTMLSelectElement).value as Params["presence"],
    presenceMin: 2,
    logMode: ($("logmode") as HTMLSelectElement).value as Params["logMode"],
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
    if (($("mode") as HTMLSelectElement).value === "perm" && p.imputation !== "none") {
      state.perm = permutationFdr(state.result.proteins, p.bait.length, state.num.s0, state.num.fdr, state.num.perms, p.seed);
    }
    markStep("s2", true, `${p.bait.length} vs ${p.control.length} replicates${nb > 1 || nc > 1 ? ", pooled" : ""}`);
    renderQc(state.result);
    renderPresence();
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
  const impText = imp === "normal" ? `fill-in shift ${state.num.shift}, spread ${state.num.shrink}` : imp === "none" ? "no fill-in" : `${imp} imputation (experimental)`;
  const pres = ($("presence") as HTMLSelectElement).value === "separate" ? "bait-only proteins listed separately" : "bait-only proteins filled in";
  return `msVolcano ${__APP_VERSION__} | ${($("test") as HTMLSelectElement).value === "welch" ? "Welch" : "Student"} t-test | ${base} | ${impText} | ${pres}`;
}

async function draw(forceLight = false) {
  const r = state.result, table = state.table;
  if (!r || !table || !state.params) return;
  showState("results");
  const p = { ...state.params, minFoldChange: state.num.minFoldChange, curvature: state.num.curvature };
  const mode = ($("mode") as HTMLSelectElement).value;
  void renderPca();
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
  const filledOf = (q: { imputed: boolean[] }) => q.imputed.filter(Boolean).length;
  // Points that rest on at least one filled-in value are drawn hollow so they are never mistaken for fully measured ones.
  const trace = (pts: typeof prot, color: string, size: number | number[], text: boolean, opacity: number, symbol = "circle", pos: string[] = [], hollowAware = true) => {
    const hollow = pts.map((q) => hollowAware && filledOf(q) > 0);
    return {
      type: "scattergl", mode: text ? "markers+text" : "markers",
      x: pts.map((q) => q.x), y: pts.map((q) => q.y),
      text: pts.map((q, i) => (pos[i] ? q.gene : "")), textposition: pos.map((x) => x || "top right"), textfont: { size: labelSize, color: c.ink }, cliponaxis: false,
      customdata: pts.map((q) => [q.id, q.gene, filledOf(q) ? `${filledOf(q)} of ${q.imputed.length} values filled in` : "all values measured"]),
      hovertemplate: "<b>%{customdata[1]}</b><br>%{customdata[0]}<br>difference %{x:.2f}<br>-log10 p %{y:.2f}<br>%{customdata[2]}<extra></extra>",
      marker: {
        color: pts.map((_, i) => (hollow[i] ? "rgba(0,0,0,0)" : color)), size, opacity, symbol,
        line: { color: pts.map((_, i) => (hollow[i] ? color : c.bg)), width: pts.map((_, i) => (hollow[i] ? 1.7 : symbol === "circle" && color === c.point ? 0 : 1.2)) },
      },
    };
  };
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
    trace(manualHits, c.pick, 12, true, 1, "diamond", manualPos, false),
    ...(mode === "hyperbola" ? [line(curve.right), line(curve.left)] : []),
  ], layout, { responsive: true, displaylogo: false, modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d", "toImage"] });
  if (forceLight) return;

  // Summary
  const d = r.dropped;
  const removed = d.contaminantOrReverse + d.absentInBait + d.untestable + d.tooFewValues;
  const nPresence = r.presenceOnly.length;
  const filledHits = hits.filter((q) => filledOf(q) > 0).length;
  $("stats").classList.toggle("four", nPresence > 0);
  const lp = $("legend-presence");
  lp.hidden = nPresence === 0;
  lp.textContent = `${nPresence} protein${nPresence === 1 ? "" : "s"} present only in bait, not plotted. Show list`;
  $("stats").innerHTML = `
    <div class="stat hit"><span>Significant interactors</span><b>${hits.length.toLocaleString()}</b><small>${prot.length ? ((hits.length / prot.length) * 100).toFixed(1) : "0"}% of ${prot.length.toLocaleString()} proteins tested${filledHits ? `. ${filledHits} rely on a filled-in value (hollow).` : ""}</small></div>
    ${nPresence ? `<button class="stat stat-btn" id="stat-presence" type="button"><span>Present only in bait</span><b>${nPresence.toLocaleString()}</b><small>not plotted, no p value. Open the list</small></button>` : ""}
    <div class="stat"><span>Proteins tested</span><b>${prot.length.toLocaleString()}</b><small>from ${table.rows.length.toLocaleString()} rows</small></div>
    <div class="stat"><span>Removed</span><b>${removed.toLocaleString()}</b><small>${d.contaminantOrReverse} contaminant or reverse, ${d.absentInBait} absent in bait${d.tooFewValues ? `, ${d.tooFewValues} with fewer than 2 measured values` : ""}${d.untestable ? `, ${d.untestable} untestable` : ""}</small></div>`;
  $("plot-title").textContent = baitName ? `Volcano plot: ${baitName}` : "Volcano plot";
  markStep("s3", false, `${MODE_LABEL[mode]} · ${hits.length} hits`);
  $("legend-cutoff").hidden = mode !== "hyperbola";
  const openPresence = () => { selectTab("presence"); $("tab-presence").scrollIntoView({ behavior: "smooth", block: "start" }); };
  $("legend-presence").onclick = openPresence;
  document.getElementById("stat-presence")?.addEventListener("click", openPresence);
  $("mbar-count").textContent = `${hits.length} interactors${nPresence ? ` + ${nPresence} present only` : ""}`;
  $("mbar").hidden = false;
  const notes: string[] = [];
  if (nPresence) notes.push(`${nPresence} protein${nPresence > 1 ? "s were" : " was"} seen in bait but in no control, so ${nPresence > 1 ? "they are" : "it is"} not on the volcano. See "Present only in bait". They are not counted in the ${hits.length} interactors.`);
  if (state.fmt) notes.push(`Input: ${state.fmt.label}.${state.fmt.note ? " " + state.fmt.note : ""}`);
  notes.push(r.scale.mode === "raw" ? `Intensities were read as raw values (median ${r.scale.median.toExponential(1)}) and log2 transformed${r.scale.decided === "auto" ? " (decided automatically)" : ""}.` : `Intensities were read as already log2 (median ${r.scale.median.toFixed(1)})${r.scale.decided === "auto" ? " (decided automatically)" : ""}.`);
  if (r.scale.warning) notes.push(`Check the scale: ${r.scale.warning}`);
  if (r.scale.decimalComma) notes.push("Numbers were read with a decimal comma.");
  if (state.params && (state.params.bait.length < 3 || state.params.control.length < 3)) notes.push("With fewer than 3 replicates in a group the p values are weak. Treat the hits as leads to confirm.");
  if (state.params && state.params.control.length < state.params.bait.length) notes.push("There are fewer control than bait replicates, so a protein can be missing from the controls by chance. Read the present-only list with that in mind.");
  if (mode === "perm") notes.push(state.perm ? `Experimental: permutation FDR ${state.perm.fdr}, ${state.perm.permutations} permutations, score threshold ${state.perm.threshold.toFixed(2)}. Proteins can be called inside the funnel because each one is scored on its own variability, not by the fixed curve.${state.perm.permutations < 10 ? ` Only ${state.perm.permutations} distinct label shuffles exist with this many replicates, so treat the result as exploratory.` : ""}` : "Experimental: not enough replicates for permutation FDR.");
  if (mode === "bh") notes.push(`Experimental: Benjamini-Hochberg q below ${state.num.qcut}, enriched side only.`);
  if (($("imputation") as HTMLSelectElement).value === "none") notes.push("Missing values were not filled in. Each test uses only the measured values, so proteins with fewer than two per group are not tested.");
  else if (($("imputation") as HTMLSelectElement).value !== "normal") notes.push("Experimental imputation method in use.");
  if (mode === "perm" && ($("imputation") as HTMLSelectElement).value === "none") notes.push("Permutation FDR needs filled-in values, so it did not run. Choose an imputation method in step 3.");
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
  const filled = (q: { imputed: boolean[] }) => q.imputed.filter(Boolean).length;
  const val = (q: (typeof hits)[number]): number | string => key === "gene" ? q.gene.toLowerCase() : key === "id" ? q.id : key === "y" ? q.y : key === "q" ? (q.q ?? 0) : key === "fi" ? filled(q) : key === "st" ? (state.stoich.get(q.row) ?? -Infinity) : q.x;
  rows.sort((a, b) => { const A = val(a), B = val(b); return (A < B ? -1 : A > B ? 1 : 0) * state.sort.dir; });
  const th = (k: string, label: string, tip: string, num = false) => `<th class="${num ? "num" : ""}" title="${esc(tip)}" aria-sort="${key === k ? (state.sort.dir === 1 ? "ascending" : "descending") : "none"}"><button type="button" data-sort="${k}">${label}${key === k ? `<span aria-hidden="true">${state.sort.dir === 1 ? " ↑" : " ↓"}</span>` : ""}</button></th>`;
  if (!hits.length) {
    const advice = mode === "hyperbola" ? "Lower the minimum enrichment or curvature in step 3." : mode === "perm" ? "Raise the FDR or lower s0 in Experimental methods." : "Raise the q value cutoff in Experimental methods.";
    $("hits").innerHTML = `<div class="callout" style="margin:14px;border:0">${ICON.info}<span>No proteins pass the cutoff. ${advice}</span></div>`; return;
  }
  $("hits").innerHTML = `<table><thead><tr>${th("gene", "Gene", "Gene name")}${th("id", "Protein ID", "First majority protein ID")}${th("x", "log2 difference", "Mean log2 LFQ in bait minus control", true)}${th("y", "-log10 p", "Higher means stronger evidence", true)}${th("fi", "Filled in", "Replicates where the value was missing and was filled in. 0 means every value was measured.", true)}${showQ ? th("q", "q value", "Benjamini-Hochberg adjusted p value", true) : ""}${showSt ? th("st", "Stoichiometry", "Relative to the bait, which is 1", true) : ""}</tr></thead><tbody>${rows.slice(0, state.rowsShown).map((q) => `<tr><td class="gene">${esc(q.gene)}</td><td>${esc(q.id)}</td><td class="num">${q.x.toFixed(2)}</td><td class="num">${q.y.toFixed(2)}</td><td class="num">${filled(q) ? `${filled(q)} of ${q.imputed.length}` : "0"}</td>${showQ ? `<td class="num">${(q.q ?? 0).toExponential(1)}</td>` : ""}${showSt ? `<td class="num">${fmtSt(state.stoich.get(q.row))}</td>` : ""}</tr>`).join("")}</tbody></table>${rows.length > state.rowsShown ? `<div class="notes" style="display:flex;gap:12px;align-items:center">Showing ${state.rowsShown} of ${rows.length}. <button class="btn small" id="more" type="button">Show 50 more</button></div>` : ""}`;
  $("more")?.addEventListener("click", () => { state.rowsShown += 50; renderHits(); });
  $("hits").querySelectorAll<HTMLButtonElement>("button[data-sort]").forEach((b) => b.addEventListener("click", () => {
    const k = b.dataset.sort!;
    state.sort = { key: k, dir: state.sort.key === k ? (state.sort.dir === 1 ? -1 : 1) : k === "gene" || k === "id" ? 1 : -1 };
    renderHits();
  }));
}
$("hit-filter").addEventListener("input", (e) => { state.filter = (e.target as HTMLInputElement).value; state.rowsShown = 25; renderHits(); });

function renderPresence() {
  const r = state.result;
  if (!r) return;
  const rows = r.presenceOnly;
  $("presence-count").textContent = `(${rows.length})`;
  if (!rows.length) {
    $("presenceTable").innerHTML = `<div class="callout" style="margin:14px;border:0">${ICON.info}<span>${($("presence") as HTMLSelectElement).value === "include" ? "Bait-only proteins are being filled in and tested, as in the published method. Switch step 3, Advanced, to list them separately." : "No protein was seen in two or more bait replicates and in no control."}</span></div>`;
    return;
  }
  const shown = rows.slice(0, 200);
  $("presenceTable").innerHTML = `<table><thead><tr><th>Gene</th><th>Protein ID</th><th class="num" title="Bait replicates with a value">Seen in</th><th class="num" title="Mean log2 intensity over the bait replicates that have a value">Mean log2 in bait</th></tr></thead><tbody>${shown.map((q) => `<tr><td class="gene">${esc(q.gene)}</td><td>${esc(q.id)}</td><td class="num">${q.nObs} of ${q.nBait}</td><td class="num">${q.meanLog2.toFixed(2)}</td></tr>`).join("")}</tbody></table>${rows.length > shown.length ? `<div class="notes">Showing the first ${shown.length} of ${rows.length}. Export the CSV for the full list.</div>` : ""}`;
}

const checkIcon = (st: Check["status"]) => (st === "ok" ? svg('<path d="m5 12 5 5 9-10"/>', 'width="16" height="16"') : st === "warn" ? svg('<path d="M12 9v4M12 17h.01"/><circle cx="12" cy="12" r="10"/>', 'width="16" height="16"') : svg('<circle cx="12" cy="12" r="10"/><path d="M12 16v-4M12 8h.01"/>', 'width="16" height="16"'));

function renderQc(r: VolcanoResult) {
  const p = state.params!;
  const shorten = (n: string) => { const s = state.clean(n); return s.length > 18 ? `${s.slice(0, 8)}...${s.slice(-7)}` : s; };
  const names = [...p.bait, ...p.control].map(shorten);
  state.pcaNames = names;
  const nb = p.bait.length;
  const c = observedCorrelation(r.proteins);

  // Replicate agreement: the lowest correlation between two replicates of the same group, from measured values only.
  let worst = 1, worstPair = "";
  const scan = (lo: number, hi: number) => { for (let i = lo; i < hi; i++) for (let j = i + 1; j < hi; j++) if (Number.isFinite(c[i][j]) && c[i][j] < worst) { worst = c[i][j]; worstPair = `${names[i]} and ${names[j]}`; } };
  scan(0, nb); scan(nb, names.length);
  const agree: Check = worst < 0.8
    ? { id: "agree", status: "warn", title: "Replicate agreement", detail: `${worstPair} agree only moderately (r = ${worst.toFixed(2)}). One may be an outlier. You can untick it in step 2.` }
    : { id: "agree", status: "ok", title: "Replicate agreement", detail: `Replicates of the same group agree well. The lowest correlation within a group is ${worst.toFixed(2)}.` };
  state.baseChecks = [agree, coverageCheck(r.coverage), loadingCheck(r.proteins, nb)];

  const minOff = Math.min(...c.flatMap((row, i) => row.filter((v, j) => j !== i && Number.isFinite(v))));
  const cell = (v: number) => `style="background:color-mix(in oklch, var(--color-accent) ${Number.isFinite(v) ? Math.round(Math.max(0, Math.min(1, (v - minOff) / Math.max(1e-6, 1 - minOff))) * 40) : 0}%, transparent)"`;
  const bars = r.coverage.map((k, i) => { const f = k.total ? k.observed / k.total : 0; return `<tr><th scope="row">${esc(names[i])}</th><td style="width:60%"><div class="bar"><i style="width:${(f * 100).toFixed(1)}%"></i></div></td><td class="num">${(f * 100).toFixed(0)}%</td></tr>`; }).join("");
  state.pca = samplePca(r.proteins);
  const shift = state.num.shift;
  $("qcBody").innerHTML = `
    <div class="qc-section">
      <div class="qc-head"><h3>Checks</h3><span class="pill new">New in v2</span></div>
      <ul class="checks" id="checks"></ul>
      <details class="qc-how"><summary>How these are judged</summary><ul><li>Replicate agreement: warns when two replicates of one group correlate below 0.8.</li><li>Quantified proteins: warns when one sample has a value for far fewer proteins than the other samples of its group (control samples normally have fewer).</li><li>Loading balance: warns when the typical protein differs by more than 1 log2 unit between bait and control.</li><li>Bait recovery: needs the bait name in step 4. Warns when the bait is outside the top 1% of proteins (at least the top 5) by enrichment.</li></ul><p>These are prompts to look closer, not statistical tests.</p></details>
    </div>
    <div class="qc-section">
      <div class="qc-head"><h3>Do the hits change if the filled-in values are lower or higher?</h3><span class="pill new">New in v2</span></div>
      <p class="qc-note">Re-runs the analysis with the fill-in shift ${shift} moved down and up by 0.6, and with two other random draws, and compares the interactor lists, using the hyperbolic curve.</p>
      <div class="row-actions" style="padding:0 0 4px"><button class="btn small" id="run-sens" type="button">Run the check</button></div>
      <div id="sens" style="padding:10px 0 0"></div>
    </div>
    <div class="qc-section">
      <div class="qc-head"><h3>Sample map (PCA)</h3><span class="pill new">New in v2</span></div>
      ${state.pca ? `<div id="pca"></div><p class="qc-note">Each dot is one sample, drawn to true scale. Bait and control should separate along PC1 (${(state.pca.explained[0] * 100).toFixed(0)}% of variance). PC2 explains only ${(state.pca.explained[1] * 100).toFixed(0)}%, so spread up and down matters much less. Uses only the ${state.pca.proteinsUsed.toLocaleString()} proteins measured in every sample, so filled-in values do not shape it. Proteins missing from the controls, which are often the specific interactors, are left out, so this map shows the background more than the interactors.</p>` : `<p class="qc-note">Not enough proteins were measured in every sample to draw a PCA (it needs at least 10).</p>`}
    </div>
    <div class="qc-section">
      <details class="qc-more"><summary>Replicate correlation and proteins with a value per sample</summary>
        <div class="qc-head" style="margin-top:12px"><h3>Replicate correlation</h3></div>
        <div class="tablewrap" tabindex="0" role="region" aria-label="Replicate correlation table, scrollable">
        <table class="qc"><thead><tr><th><span class="sr">Sample</span></th>${names.map((n) => `<th scope="col">${esc(n)}</th>`).join("")}</tr></thead><tbody>${c.map((row, i) => `<tr><th scope="row">${esc(names[i])}</th>${row.map((v) => `<td ${cell(v)}>${Number.isFinite(v) ? v.toFixed(2) : ""}</td>`).join("")}</tr>`).join("")}</tbody></table></div>
        <p class="qc-note">Pearson correlation of log2 intensities between samples, using proteins measured in both. A deeper tint means a higher correlation, scaled from the lowest value shown to 1.</p>
        <div class="qc-head" style="margin-top:14px"><h3>Proteins with a value, per sample</h3></div>
        <table class="cov"><tbody>${bars}</tbody></table>
      </details>
    </div>`;
  paintChecks();
  $("run-sens").addEventListener("click", runSensitivity);
}

/** The bait check depends on the bait name, so it is repainted whenever that field changes. */
function paintChecks() {
  const r = state.result;
  const list = document.getElementById("checks");
  if (!r || !list) return;
  const baitKey = (($("bait-name") as HTMLInputElement).value || ($("stoichBait") as HTMLInputElement).value).trim();
  const checks: Check[] = [...state.baseChecks, baitRecovery(r.proteins, r.presenceOnly, baitKey)];
  list.innerHTML = checks.map((k) => `<li class="check-row ${k.status}"><span class="ci" aria-hidden="true">${checkIcon(k.status)}</span><div><b>${esc(k.title)}</b><div>${esc(k.detail)}</div></div><span class="sr">${k.status === "ok" ? "Passed" : k.status === "warn" ? "Needs attention" : "Note"}</span></li>`).join("");
  const warn = checks.filter((k) => k.status === "warn").length;
  $("qc-dot").innerHTML = warn ? `<span class="dot-warn" role="img" aria-label="${warn} quality warning${warn > 1 ? "s" : ""}"></span>` : `<span class="all-pass" role="img" aria-label="All ${checks.length} checks passed">${svg('<path d="m5 12 5 5 9-10"/>', 'width="13" height="13"')}</span>`;
}
for (const id of ["bait-name", "stoichBait"]) $(id).addEventListener("input", () => paintChecks());

function runSensitivity() {
  const out = $("sens");
  if (!state.table || !state.params) return;
  const pr = { ...state.params, minFoldChange: state.num.minFoldChange, curvature: state.num.curvature };
  const base = pr.shift;
  const shifts = [Math.max(0, +(base - 0.6).toFixed(2)), base, +(base + 0.6).toFixed(2)];
  const rows = shiftSensitivity(state.table, pr, shifts, [pr.seed + 1, pr.seed + 2]);
  if (!rows) { out.innerHTML = `<div class="callout">${ICON.info}<span>Nothing is filled in with the current method, so there is nothing to test.</span></div>`; return; }
  const baseHits = rows.find((x) => x.base)!.hits;
  const worstChange = Math.max(...rows.map((x) => x.lost + x.gained));
  const stable = baseHits === 0 ? worstChange === 0 : worstChange / baseHits <= 0.1;
  out.innerHTML = `<div class="qc-verdict ${stable ? "" : "warn"}" style="margin:0 0 10px;border:1px solid var(--color-line);border-radius:var(--radius-s)"><span>${stable ? "Stable. Lower or higher filled-in values and other random draws change at most 10% of the interactor list (cutoff: enrichment " + pr.minFoldChange + ", curvature " + pr.curvature + ")." : `Sensitive. Lower or higher filled-in values or other random draws change up to ${baseHits ? Math.round((worstChange / baseHits) * 100) : 0}% of the interactor list. Treat borderline hits with caution and check them against the measured values.`}</span></div>
    <div class="tablewrap" tabindex="0" role="region" aria-label="Sensitivity results, scrollable"><table><thead><tr><th>Setting</th><th class="num">Interactors</th><th class="num">Shared with current</th><th class="num">Dropped</th><th class="num">Added</th></tr></thead><tbody>${rows.map((x) => `<tr><td>${esc(x.label)}</td><td class="num">${x.hits}</td><td class="num">${x.shared}</td><td class="num">${x.lost}</td><td class="num">${x.gained}</td></tr>`).join("")}</tbody></table></div>`;
}

async function renderPca() {
  const el = document.getElementById("pca");
  if (!el || !state.pca || !state.params) return;
  const css = (n: string) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
  const ink = css("--plot-ink"), grid = css("--plot-grid"), axis = css("--plot-axis"), bg = css("--plot-bg");
  const nb = state.params.bait.length;
  const mk = (idx: number[], name: string, color: string, symbol: string) => ({
    type: "scatter", mode: "markers+text", name, x: idx.map((i) => state.pca!.scores[i][0]), y: idx.map((i) => state.pca!.scores[i][1]),
    text: idx.map((i) => state.pcaNames[i]), textposition: idx.map((_, k) => (k % 2 ? "bottom center" : "top center")), textfont: { size: 11, color: ink },
    marker: { size: 11, color, symbol, line: { color: bg, width: 1.2 } },
    hovertemplate: "%{text}<extra>" + name + "</extra>",
  });
  const all = state.pcaNames.map((_, i) => i);
  const pct = (v: number) => `${(v * 100).toFixed(0)}%`;
  await Plotly.react(el, [mk(all.filter((i) => i < nb), "Bait", css("--plot-hit"), "circle"), mk(all.filter((i) => i >= nb), "Control", css("--plot-curve"), "square")], {
    font: { color: ink, family: "Geist Variable, system-ui, sans-serif", size: 12 },
    xaxis: { title: { text: `PC1 (${pct(state.pca.explained[0])} of variance)` }, zeroline: false, gridcolor: grid, showline: true, linecolor: axis, constrain: "domain" },
    yaxis: { title: { text: `PC2 (${pct(state.pca.explained[1])})` }, zeroline: false, gridcolor: grid, showline: true, linecolor: axis, scaleanchor: "x", scaleratio: 1, constrain: "domain" },
    showlegend: true, legend: { orientation: "h", y: 1.12 }, margin: { t: 30, r: 20, b: 50, l: 56 }, paper_bgcolor: bg, plot_bgcolor: bg, hovermode: "closest",
  }, { responsive: true, displaylogo: false, displayModeBar: false });
}

/* ---------- Tabs ---------- */

type TabId = "hits" | "presence" | "qc";
const TABS: TabId[] = ["hits", "presence", "qc"];
function selectTab(which: TabId, focus = false) {
  for (const t of TABS) {
    $(`tab-${t}`).setAttribute("aria-selected", String(t === which));
    $(`tab-${t}`).tabIndex = t === which ? 0 : -1;
    $(`panel-${t}`).hidden = t !== which;
  }
  if (focus) $(`tab-${which}`).focus();
  if (which === "qc") void renderPca();
}
for (const t of TABS) $(`tab-${t}`).addEventListener("click", () => selectTab(t));
document.querySelector(".subtabs")!.addEventListener("keydown", (e) => {
  const k = (e as KeyboardEvent).key;
  const cur = TABS.findIndex((t) => $(`tab-${t}`).getAttribute("aria-selected") === "true");
  if (k === "ArrowRight") selectTab(TABS[(cur + 1) % TABS.length], true);
  else if (k === "ArrowLeft") selectTab(TABS[(cur + TABS.length - 1) % TABS.length], true);
  else if (k === "Home") selectTab(TABS[0], true);
  else if (k === "End") selectTab(TABS[TABS.length - 1], true);
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
async function copyGenes(list: { gene: string }[] = state.hits) {
  closeMenu();
  const text = list.map((q) => q.gene).join("\n");
  try { await navigator.clipboard.writeText(text); toast(`Copied ${list.length} gene names`); } catch { toast("Copy is blocked by the browser"); }
}
function downloadPresence() {
  closeMenu();
  if (!state.table || !state.params || !state.result) return;
  const url = URL.createObjectURL(new Blob([presenceCsv(state.table, state.result.presenceOnly, state.params)], { type: "text/csv" }));
  const a = Object.assign(document.createElement("a"), { href: url, download: `${stem()}_present_only.csv` });
  a.click(); URL.revokeObjectURL(url);
}
$("csvp").addEventListener("click", downloadPresence);
$("csvp2").addEventListener("click", downloadPresence);
$("copy3").addEventListener("click", () => void copyGenes(state.result?.presenceOnly ?? []));
$("copyall").addEventListener("click", () => void copyGenes([...state.hits, ...(state.result?.presenceOnly ?? [])]));
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
    so: ($("stoich") as HTMLInputElement).checked, sg: selVal("organism"), sb: selVal("stoichBait"), pr: selVal("presence"), lm: selVal("logmode"),
  };
  return btoa(unescape(encodeURIComponent(JSON.stringify(saved)))).replace(/\+/g, "-").replace(/\//g, "_").replace(/=+$/, "");
}
function decodeSettings(raw: string): Saved | null {
  try { return JSON.parse(decodeURIComponent(escape(atob(raw.replace(/-/g, "+").replace(/_/g, "/"))))) as Saved; } catch { return null; }
}
function applySavedControls(s: Saved) {
  for (const [k, v] of Object.entries(s.n ?? {})) if (k in state.num && Number.isFinite(v)) setters[k as keyof typeof state.num]?.(v);
  const set = (id: string, v: string | undefined) => { if (v !== undefined) ($(id) as HTMLSelectElement | HTMLInputElement).value = v; };
  set("mode", s.m); set("test", s.t); set("imputation", s.i); set("presence", s.pr); set("logmode", s.lm); set("manual", s.man); set("title", s.ti); set("bait-name", s.b); set("stoichBait", s.sb);
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
  ($("mode") as HTMLSelectElement).value = "hyperbola"; ($("test") as HTMLSelectElement).value = "student"; ($("imputation") as HTMLSelectElement).value = "normal"; ($("presence") as HTMLSelectElement).value = "separate"; ($("logmode") as HTMLSelectElement).value = "auto";
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
