# msVolcano

Volcano plots for label free interactomics data from MaxQuant. Load a `proteinGroups.txt`, pick bait and control LFQ columns, tune the hyperbolic cutoff, and export the plot and the list of interactors. Everything runs in your browser, so your data stays on your machine.

This is version 2, a rewrite of the Shiny app described in the original paper. The old server at TU Dresden BIOTEC is gone.

## What it does

* Removes contaminants and reverse hits, log2 transforms unlogged LFQ values.
* Imputes missing values from a shifted normal distribution (`shift`, `shrink`), seeded so results repeat.
* Student or Welch t-test per protein, volcano plot with the hyperbolic cutoff `y = curvature / (x - minFoldChange)` from Keilhauer, Hein and Mann (2015).
* Labels significant proteins and any gene names you list. Exports PNG, SVG and a CSV of hits.

Not yet ported: iBAQ style stoichiometry and loading files from an ftp URL.

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

GitHub shows a "Cite this repository" button from `CITATION.cff`.

## License

Free for noncommercial use under the [PolyForm Noncommercial License 1.0.0](LICENSE). Research, teaching and nonprofit use are covered. If you copy or adapt the code, keep the `NOTICE` file with it. This is source available software, not OSI open source.

For commercial use, contact the author for a commercial license.
