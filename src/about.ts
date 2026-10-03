// About page content. Plain HTML string, rendered next to the tool.
export const aboutHtml = `
<article class="about">
  <h2>About msVolcano</h2>

  <h3>What it does</h3>
  <p>msVolcano turns the output of a label free affinity purification or affinity enrichment mass spectrometry experiment (AP/MS, AE/MS) into a volcano plot and a list of candidate interactors. You give it a MaxQuant <code>proteinGroups.txt</code> that was processed with MaxLFQ, choose which LFQ columns are the bait replicates and which are the controls, and adjust a hyperbolic cutoff until the true interactors separate from the background binders.</p>

  <h3>Why it was built</h3>
  <p>In an interactomics experiment, replicates of the affinity enriched bait are compared with negative controls. Proteins that bind non specifically sit around zero on the plot, and enriched interactors move to the right. Choosing the threshold that separates the two is the critical step and usually needs some manual tuning. At the time (2016) this analysis meant specialist desktop software and some scripting, which was a hurdle for bench scientists and a burden for mass spectrometry core facilities. msVolcano put every downstream step behind one simple interface that needs no bioinformatics knowledge. It was developed in the Stewart lab at the Biotechnology Center (BIOTEC) of TU Dresden and published in <em>Proteomics</em>.</p>

  <h3>What happens to your data</h3>
  <ol>
    <li>Rows marked as reverse (decoy) hits or potential contaminants are removed.</li>
    <li>Proteins that are absent from all bait replicates are removed.</li>
    <li>LFQ values that look unlogged are transformed to log2.</li>
    <li>Missing values are imputed by drawing from a normal distribution shifted down from the group mean (<code>shift</code> times the standard deviation) and narrowed (<code>shrink</code> times the standard deviation), which mimics low abundance proteins that fell below detection. The draw is seeded, so the same input gives the same plot.</li>
    <li>A Student or Welch t-test is run for every protein, bait against control.</li>
    <li>The volcano plot shows the difference of the means (x) against the negative log10 p value (y). A protein is called significant when it lies to the right of the hyperbolic curve <code>y = curvature / (x - minFoldChange)</code>, which tightens the p value requirement for weakly enriched proteins and relaxes it for strongly enriched ones.</li>
  </ol>
  <p>Your file is read by your browser and processed on your computer. Nothing is uploaded to any server.</p>

  <h3>Stoichiometry</h3>
  <p>Optionally, msVolcano estimates how abundant each enriched protein is relative to the bait, as described in the 2016 paper. The intensity above the control is divided by the number of theoretical tryptic peptides (7 to 30 amino acids) of the protein, then divided by the same quantity for the bait. The peptide tables for nine organisms were rebuilt from UniProt Swiss-Prot for version 2. They count tryptic peptides with the standard rule that trypsin does not cleave before proline, whereas the original tool cleaved there too, so absolute values can differ from the original. Compare numbers only within one analysis.</p>

  <h3>Input</h3>
  <p>A tab separated MaxQuant <code>proteinGroups.txt</code>, or a comma separated file with the same column names. It needs at least two bait and two control columns whose names contain "LFQ". Gene names and majority protein IDs are used for labels when present. Use the example dataset (simulated, not real data) to try the tool.</p>

  <h3>Experimental features</h3>
  <p>Features marked <span class="badge">Experimental</span> are new in version 2. They are <strong>not part of the 2016 publication</strong>, have not been peer reviewed or benchmarked, and may change or disappear. They are offered as is, open for anyone to try, test and improve. Please check results from them against your own judgement and an established tool before drawing conclusions, and report problems on GitHub.</p>
  <ul>
    <li><strong>s0 score with permutation FDR (Perseus style).</strong> Scores each protein as the difference divided by its standard error plus a fudge factor s0, shuffles the sample labels to estimate how many false calls to expect, and picks the threshold that keeps the false discovery rate at your chosen level. Idea from Tusher et al. 2001, as used in Perseus (Tyanova et al. 2016). Only proteins enriched in the bait are called.</li>
    <li><strong>Benjamini-Hochberg q value.</strong> Multiple testing correction of the t-test p values (Benjamini and Hochberg 1995), with a q value cutoff you choose.</li>
    <li><strong>Alternative imputation.</strong> MinDet and MinProb replace missing values using a low quantile of each column, following the comparison by Lazar et al. 2016 for values missing because they fell below detection.</li>
    <li><strong>Replicate correlation.</strong> A quick quality check of how well replicates agree.</li>
  </ul>
  <p>Ideas not built yet: input from FragPipe, DIA-NN and Spectronaut, moderated (limma style) t-tests, Fisher exact test for proteins seen only in bait, complex enrichment with CORUM. Pull requests and issues are welcome.</p>

  <h3>How to cite</h3>
  <p>If you use msVolcano, please cite the paper:</p>
  <p class="cite">Singh S, Hein MY, Stewart AF. msVolcano: A flexible web application for visualizing quantitative proteomics data. <em>Proteomics</em> 2016;16(18):2491. <a href="https://doi.org/10.1002/pmic.201600167">doi:10.1002/pmic.201600167</a> (<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC5096246/">PMC5096246</a>, <a href="https://pubmed.ncbi.nlm.nih.gov/27440201/">PubMed</a>)</p>
  <p>The repository also has a <code>CITATION.cff</code>, so GitHub offers a "Cite this repository" button. A software DOI will be added when the first release is archived.</p>

  <h3>References</h3>
  <ol class="refs">
    <li>Keilhauer EC, Hein MY, Mann M. Accurate protein complex retrieval by affinity enrichment mass spectrometry (AE-MS) rather than affinity purification mass spectrometry (AP-MS). <em>Mol Cell Proteomics</em> 2015;14(1):120-135. (hyperbolic cutoff)</li>
    <li>Cox J, Hein MY, Luber CA, Paron I, Nagaraj N, Mann M. Accurate proteome-wide label-free quantification by delayed normalization and maximal peptide ratio extraction, termed MaxLFQ. <em>Mol Cell Proteomics</em> 2014;13(9):2513-2526.</li>
    <li>Cox J, Mann M. MaxQuant enables high peptide identification rates, individualized p.p.b.-range mass accuracies and proteome-wide protein quantification. <em>Nat Biotechnol</em> 2008;26(12):1367-1372.</li>
    <li>Hein MY, Hubner NC, Poser I, et al. A human interactome in three quantitative dimensions organized by stoichiometries and abundances. <em>Cell</em> 2015;163(3):712-723. (stoichiometry)</li>
    <li>Schwanhäusser B, Busse D, Li N, et al. Global quantification of mammalian gene expression control. <em>Nature</em> 2011;473(7347):337-342. (iBAQ)</li>
    <li>Tusher VG, Tibshirani R, Chu G. Significance analysis of microarrays applied to the ionizing radiation response. <em>Proc Natl Acad Sci USA</em> 2001;98(9):5116-5121.</li>
  </ol>

  <h3>License</h3>
  <p>msVolcano is free for noncommercial use under the <a href="https://polyformproject.org/licenses/noncommercial/1.0.0">PolyForm Noncommercial License 1.0.0</a>. Research, teaching, personal study and nonprofit use are covered. If you copy, modify or redistribute the code, you must keep the <a href="https://github.com/uksurd88/msvolcano/blob/main/NOTICE"><code>NOTICE</code></a> file, which names the author and asks users to cite the paper. This is source available software, not OSI approved open source.</p>
  <p>Commercial use needs a separate commercial license. Contact the author through the profile below.</p>

  <h3>History and author</h3>
  <p>The original web app ran on a server at BIOTEC, TU Dresden. That server was retired and the link in the paper stopped working. Version 2 is a rewrite, not a copy: the statistics were reimplemented in TypeScript and checked against reference values from SciPy. Behaviour follows the description in the paper.</p>
  <p>Author: Sukhdeep Singh. <a href="https://scholar.google.com/citations?user=y62wfS8AAAAJ&hl=en">Google Scholar profile</a>. Source code and issue tracker: <a href="https://github.com/uksurd88/msvolcano">github.com/uksurd88/msvolcano</a>.</p>
</article>`;
