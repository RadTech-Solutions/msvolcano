// About page content. Plain HTML string, rendered next to the tool.
export const aboutHtml = `
<article class="about">
  <h2>About msVolcano</h2>

  <h3>Two versions</h3>
  <p><strong>Version 1</strong> is the method published in <em>Proteomics</em> in 2016, reimplemented to run in the browser: MaxQuant <code>proteinGroups.txt</code> with LFQ intensities, a Student or Welch t-test, missing values filled in from a shifted normal distribution (bait-only proteins included), the hyperbolic cutoff, and stoichiometry. Its settings are locked to the published behaviour. Use it to reproduce or cite published results. Open it directly with <a href="?v=1#tool">?v=1</a>.</p>
  <p><strong>Version 2</strong> is <strong>experimental and not validated</strong>. It adds more input formats, lists bait-only proteins separately, offers other fill-in methods and FDR based cutoffs, and adds quality checks. None of this was peer reviewed. Results can differ from version 1, so compare before relying on them. Open it with <a href="?v=2#tool">?v=2</a>.</p>
  <p>Both are the same code base. The steps below describe version 2. Where version 1 differs, it does what the 2016 tool did.</p>

  <h3>What it does</h3>
  <p>msVolcano turns the output of a label free affinity purification or affinity enrichment mass spectrometry experiment (AP/MS, AE/MS) into a volcano plot and a list of candidate interactors. You give it a protein table (MaxQuant with LFQ, FragPipe, DIA-NN and others), choose which sample columns are the bait replicates and which are the controls, and adjust a hyperbolic cutoff until the true interactors separate from the background binders.</p>

  <h3>Why it was built</h3>
  <p>In an interactomics experiment, replicates of the affinity enriched bait are compared with negative controls. Proteins that bind non specifically sit around zero on the plot, and enriched interactors move to the right. Choosing the threshold that separates the two is the critical step and usually needs some manual tuning. At the time (2016) this analysis meant specialist desktop software and some scripting, which was a hurdle for bench scientists and a burden for mass spectrometry core facilities. msVolcano put every downstream step behind one simple interface that needs no bioinformatics knowledge. It was developed in the Stewart lab at the Biotechnology Center (BIOTEC) of TU Dresden and published in <em>Proteomics</em>.</p>

  <h3>What happens to your data</h3>
  <ol>
    <li>Rows marked as reverse (decoy) hits or potential contaminants are removed.</li>
    <li>Proteins with no value in any bait replicate are removed.</li>
    <li>Intensities are read as raw or as log2. msVolcano decides from the typical value, tells you what it decided under the plot, and lets you set it in step 3. In raw intensities a zero means "not quantified" and counts as missing. In log2 data an empty cell is missing.</li>
    <li>Proteins seen in at least 2 bait replicates and not detected in any control replicate are set aside in their own list (see "New in version 2" below). Proteins with fewer than 2 measured values in both groups are left out, because they would be tested on filled-in values alone. The 2016 tool filled both kinds in and tested them.</li>
    <li>Remaining missing values are filled in by drawing from a normal distribution shifted down from the group mean (<code>shift</code> times the standard deviation) and narrowed (<code>spread</code> times the standard deviation), which mimics low abundance proteins that fell below detection. As in 2016, bait and control each get their own distribution. Every filled-in value has its own seeded draw, so the same input always gives the same plot and changing a setting elsewhere does not reshuffle other proteins. Points that rest on a filled-in value are drawn hollow.</li>
    <li>A Student or Welch t-test is run for every protein, bait against control.</li>
    <li>The volcano plot shows the difference of the means (x) against the negative log10 p value (y). A protein is called significant when it lies to the right of the hyperbolic curve <code>y = curvature / (x - minFoldChange)</code>, which tightens the p value requirement for weakly enriched proteins and relaxes it for strongly enriched ones. This curve is a tuning device, not a calibrated false discovery rate.</li>
  </ol>
  <p>Your file is read by your browser and processed on your computer. Nothing is uploaded to any server.</p>

  <h3>Stoichiometry</h3>
  <p>Optionally, msVolcano estimates how abundant each enriched protein is relative to the bait, as described in the 2016 paper. The intensity above the control is divided by the number of theoretical tryptic peptides (7 to 30 amino acids) of the protein, then divided by the same quantity for the bait. The peptide tables for nine organisms were rebuilt from UniProt Swiss-Prot for this rewrite and are used by both versions. They count tryptic peptides with the standard rule that trypsin does not cleave before proline, whereas the original tool cleaved there too, so absolute values can differ from the original. Compare numbers only within one analysis.</p>

  <h3>Input</h3>
  <p>A protein table with one intensity column per sample, tab or comma separated. Detected automatically:</p>
  <ul>
    <li>MaxQuant <code>proteinGroups.txt</code> with <strong>LFQ intensity</strong> columns.</li>
    <li>FragPipe <code>combined_protein.tsv</code> (MaxLFQ intensity columns).</li>
    <li>DIA-NN <code>report.pg_matrix.tsv</code>, Spectronaut protein group pivot, and Proteome Discoverer protein exports. These three readers are <strong>beta</strong>: written from the software documentation and not yet checked on many real files. Check the groups in step 2 and report problems on GitHub.</li>
    <li>Perseus matrices (the columns marked Main) and any plain table of numbers.</li>
  </ul>
  <p>Gene names and protein IDs are used for labels when present. Use the example dataset (simulated, not real data) to try the tool.</p>

  <h3>New in version 2 <span class="badge">Not in the 2016 paper</span></h3>
  <p>These additions follow what the field has learned since 2016. They change what you see, so they are described here. Each can be switched off or reverted to the published behaviour.</p>
  <ul>
    <li><strong>Present only in bait.</strong> A protein found in at least 2 bait replicates and not detected in any control is a strong candidate, but a t-test on filled-in control values ranks it using numbers that were invented, and the result depends on the shift setting. It can also be missing from the controls by chance, especially with few controls. Version 2 lists these proteins on their own tab, ranked by how many bait replicates saw them and by signal, with no p value. Imputation accuracy is driven mainly by how much of the missingness is "not at random" (<a href="https://doi.org/10.1038/s41598-021-81279-4">Jin et al. 2021</a>; <a href="https://doi.org/10.1021/acs.jproteome.5b00981">Lazar et al. 2016</a>). Choose "Fill in and test (as published)" in step 3, Advanced, to return to the 2016 behaviour.</li>
    <li><strong>Visible imputation.</strong> Points that include a filled-in value are hollow, the interactors table shows how many values were filled in, and a "Do not fill in" mode tests only the measured values. A sensitivity check re-runs the analysis at nearby shift values and tells you how much of the interactor list moves.</li>
    <li><strong>Quality checks.</strong> Replicate agreement from measured values only, a sample map (PCA), the share of proteins with a value in each sample, a loading balance check, and a bait recovery check (is the bait among the most enriched proteins?). These are warnings to look at, not statistical tests.</li>
  </ul>

  <h3>Experimental features</h3>
  <p>Two labels are used. "New in v2" marks additions that follow current good practice and are described above. <span class="badge">Experimental</span> marks statistical methods that are less settled. Both are new in version 2. They are <strong>not part of the 2016 publication</strong>, have not been peer reviewed or benchmarked, and may change or disappear. They are offered as is, open for anyone to try, test and improve. Please check results from them against your own judgement and an established tool before drawing conclusions, and report problems on GitHub.</p>
  <ul>
    <li><strong>s0 score with permutation FDR (Perseus style).</strong> Scores each protein as the difference divided by its standard error plus a fudge factor s0, shuffles the sample labels to estimate how many false calls to expect, and picks the threshold that keeps the false discovery rate at your chosen level. Idea from Tusher et al. 2001, as used in Perseus (Tyanova et al. 2016). Only proteins enriched in the bait are called.</li>
    <li><strong>Benjamini-Hochberg q value.</strong> Multiple testing correction of the t-test p values (Benjamini and Hochberg 1995), with a q value cutoff you choose.</li>
    <li><strong>Alternative imputation.</strong> MinDet and MinProb replace missing values using a low quantile of each column, following the comparison by Lazar et al. 2016 for values missing because they fell below detection.</li>
  </ul>
  <p>Ideas not built yet: moderated (limma style) t-tests, SAINT or limma result import, contaminant flagging with CRAPome, complex enrichment with CORUM. Pull requests and issues are welcome.</p>

  <h3>Limits to keep in mind</h3>
  <ul>
    <li>P values computed with filled-in values are optimistic, because filled-in numbers count as independent replicates. The hollow points, the "Filled in" column and the sensitivity check exist to show you where that matters.</li>
    <li>msVolcano does not normalise your data. Use the normalised LFQ or abundance columns from your software, and check the loading balance on the Quality checks tab.</li>
    <li>With two or three replicates per group the statistics are weak whatever the method. The tool says so under the plot.</li>
    <li>The hyperbolic curve is a tuning device, not a false discovery rate. The permutation FDR and Benjamini-Hochberg options are experimental, and the permutation FDR needs enough replicates to have more than a handful of distinct label shuffles.</li>
    <li>The sample map (PCA) uses only proteins measured in every sample, so it shows the background more than the specific interactors. The sensitivity check varies the fill-in shift and the random draw, not every analysis choice.</li>
    <li>Control samples usually have fewer values than bait samples, because the tested proteins are chosen by their bait signal. The coverage check therefore compares each sample with the others of its own group.</li>
    <li>The DIA-NN, Spectronaut and Proteome Discoverer readers are beta.</li>
  </ul>

  <h3>How to cite</h3>
  <p>If you use msVolcano, please cite the paper:</p>
  <p class="cite">Singh S, Hein MY, Stewart AF. msVolcano: A flexible web application for visualizing quantitative proteomics data. <em>Proteomics</em> 2016;16(18):2491. <a href="https://doi.org/10.1002/pmic.201600167">doi:10.1002/pmic.201600167</a> (<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC5096246/">PMC5096246</a>, <a href="https://pubmed.ncbi.nlm.nih.gov/27440201/">PubMed</a>)</p>
  <p>To cite the software itself, use the archived release on Zenodo: Singh S. msVolcano. Zenodo. <a href="https://doi.org/10.5281/zenodo.23262108">doi:10.5281/zenodo.23262108</a>. This DOI always points to the latest release. The repository also has a <code>CITATION.cff</code>, so GitHub offers a "Cite this repository" button.</p>

  <h3>References</h3>
  <ol class="refs">
    <li>Keilhauer EC, Hein MY, Mann M. Accurate protein complex retrieval by affinity enrichment mass spectrometry (AE-MS) rather than affinity purification mass spectrometry (AP-MS). <em>Mol Cell Proteomics</em> 2015;14(1):120-135. (hyperbolic cutoff)</li>
    <li>Cox J, Hein MY, Luber CA, Paron I, Nagaraj N, Mann M. Accurate proteome-wide label-free quantification by delayed normalization and maximal peptide ratio extraction, termed MaxLFQ. <em>Mol Cell Proteomics</em> 2014;13(9):2513-2526.</li>
    <li>Cox J, Mann M. MaxQuant enables high peptide identification rates, individualized p.p.b.-range mass accuracies and proteome-wide protein quantification. <em>Nat Biotechnol</em> 2008;26(12):1367-1372.</li>
    <li>Hein MY, Hubner NC, Poser I, et al. A human interactome in three quantitative dimensions organized by stoichiometries and abundances. <em>Cell</em> 2015;163(3):712-723. (stoichiometry)</li>
    <li>Schwanhäusser B, Busse D, Li N, et al. Global quantification of mammalian gene expression control. <em>Nature</em> 2011;473(7347):337-342. (iBAQ)</li>
    <li>Jin L, Bi Y, Hu C, et al. A comparative study of evaluating missing value imputation methods in label-free proteomics. <em>Sci Rep</em> 2021;11:1760. <a href="https://doi.org/10.1038/s41598-021-81279-4">doi:10.1038/s41598-021-81279-4</a></li>
    <li>Lazar C, Gatto L, Ferro M, Bruley C, Burger T. Accounting for the multiple natures of missing values in label-free quantitative proteomics data sets to compare imputation strategies. <em>J Proteome Res</em> 2016;15(4):1116-1125. <a href="https://doi.org/10.1021/acs.jproteome.5b00981">doi:10.1021/acs.jproteome.5b00981</a></li>
    <li>Tusher VG, Tibshirani R, Chu G. Significance analysis of microarrays applied to the ionizing radiation response. <em>Proc Natl Acad Sci USA</em> 2001;98(9):5116-5121.</li>
  </ol>

  <h3>License</h3>
  <p>msVolcano is free for noncommercial use under the <a href="https://polyformproject.org/licenses/noncommercial/1.0.0">PolyForm Noncommercial License 1.0.0</a>. Research, teaching, personal study and nonprofit use are covered. If you copy, modify or redistribute the code, you must keep the <a href="https://github.com/uksurd88/msvolcano/blob/main/NOTICE"><code>NOTICE</code></a> file, which names the author and asks users to cite the paper. This is source available software, not OSI approved open source.</p>
  <p>Commercial use needs a separate commercial license. Contact the author through the profile below.</p>

  <h3>History and author</h3>
  <p>The original web app ran on a server at BIOTEC, TU Dresden. That server was retired and the link in the paper stopped working. This site is a rewrite, not a copy: the statistics were reimplemented in TypeScript, and the t-test and p value arithmetic was checked against reference values from SciPy. Behaviour follows the description in the paper except where listed under "New in version 2".</p>
  <p>Author: Sukhdeep Singh. <a href="https://scholar.google.com/citations?user=y62wfS8AAAAJ&hl=en">Google Scholar profile</a>. Source code and issue tracker: <a href="https://github.com/uksurd88/msvolcano">github.com/uksurd88/msvolcano</a>.</p>
</article>`;
