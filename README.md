# msVolcano

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23262108.svg)](https://doi.org/10.5281/zenodo.23262108)

Volcano plots for label free interactomics data from MaxQuant. Load a `proteinGroups.txt`, pick bait and control LFQ columns, tune the hyperbolic cutoff, and export the plot and the list of interactors. Everything runs in your browser, so your data stays on your machine.

Live: https://radtech-solutions.github.io/msvolcano/ (maintained by [RadTech](https://radtech.nl))

A browser rewrite of the Shiny app described in the original paper. The old server at TU Dresden BIOTEC is gone.

The site has two versions, chosen with the switch at the top of the page:

* **Version 1** (`?v=1`): the method published in Proteomics 2016. MaxQuant LFQ input, Student or Welch t-test, shifted normal fill-in, hyperbolic cutoff, stoichiometry. Settings are locked to the published behaviour.
* **Version 2** (`?v=2`): experimental, not validated. Everything below that is marked as new.

## What it does

* Reads MaxQuant (LFQ), FragPipe, DIA-NN, Spectronaut and Proteome Discoverer tables, Perseus matrices and plain numeric tables. The DIA-NN, Spectronaut and Proteome Discoverer readers are beta.
* Removes contaminants and reverse hits, log2 transforms unlogged intensities.
* Imputes missing values from a shifted normal distribution (`shift`, `shrink`), seeded so results repeat.
* Student or Welch t-test per protein, volcano plot with the hyperbolic cutoff `y = curvature / (x - minFoldChange)` from Keilhauer, Hein and Mann (2015).
* Labels significant proteins and any gene names you list. Exports PNG, SVG and a CSV of hits.

* Optional iBAQ style stoichiometry relative to the bait. Tryptic peptide tables for nine organisms are rebuilt from UniProt Swiss-Prot (`npm run` is not needed, the tables are in `public/peptides`; regenerate with `node scripts/build-peptides.mjs`). They follow the standard rule of no cleavage before proline, so values can differ from the 2016 tool.
* A simulated example dataset (`scripts/make-example.mjs`) to try the app. It contains no real data.

### New in version 2 (not part of the 2016 publication)

* Proteins seen in two or more bait replicates and in no control are listed separately with replicate counts instead of being tested on filled-in numbers. The published behaviour is one setting away.
* Points that rest on filled-in values are hollow, the table says how many values were filled in, there is a "do not fill in" mode, and a sensitivity check shows how much the hit list depends on the fill-in shift.
* Quality checks: replicate agreement from measured values, PCA of samples, per-sample coverage, loading balance, and bait recovery.

### Experimental features (not part of the 2016 publication)

Clearly marked in the app and on the About page: s0 score with permutation FDR (Perseus style), Benjamini-Hochberg q values, and MinDet and MinProb imputation. They are not peer reviewed or benchmarked. Tests cover the arithmetic and a pure noise calibration, not their scientific validity.

Not ported: loading files from an ftp URL. Ideas: limma style moderated t-test, SAINT or limma result import, CRAPome contaminant flags, CORUM complex enrichment.

## Run locally

```
npm install
npm run dev     # development server
npm test        # unit tests for the statistics and the cutoff
npm run build   # static site in dist/
```

The build is a static folder. Open it from any web server, or deploy it to GitHub Pages (the workflow in `.github/workflows/pages.yml` does this on every push to `main`).

## Cite

If you use msVolcano, cite:

Singh S, Hein MY, Stewart AF. msVolcano: A flexible web application for visualizing quantitative proteomics data. Proteomics 2016;16(18):2491. https://doi.org/10.1002/pmic.201600167

To cite the software, use the Zenodo DOI [10.5281/zenodo.23262108](https://doi.org/10.5281/zenodo.23262108), which always points to the latest release. GitHub also shows a "Cite this repository" button from `CITATION.cff`.

## License

Free for noncommercial use under the [PolyForm Noncommercial License 1.0.0](LICENSE). Research, teaching and nonprofit use are covered. If you copy or adapt the code, keep the `NOTICE` file with it. This is source available software, not OSI open source.

For commercial use, contact the author for a commercial license (see the [Google Scholar profile](https://scholar.google.com/citations?user=y62wfS8AAAAJ&hl=en)).
