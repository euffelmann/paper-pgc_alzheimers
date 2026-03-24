R/ — Utility functions for the AD GWAS meta-analysis pipeline
==============================================================

Scripts are grouped below by purpose. "Job" scripts are thin
command-line wrappers that parse arguments and call the
corresponding function script; they are intended to be submitted
as cluster jobs.


Summary statistics formatting and QC
--------------------------------------

filter_sumstats.R
  Filters meta-analysed summary statistics by two criteria:
  (1) effective sample size (neff) as a fraction of the maximum
  observed neff across all SNPs, and (2) the minimum number of
  contributing cohorts per SNP, derived from the METAL direction
  string. Intended to be run as a command-line script; takes
  input file, neff threshold, and cohort-count threshold as
  arguments.

process_sumstats.R
  Core QC function for per-cohort summary statistics prior to
  meta-analysis. Applies a sequential filter pipeline: removes
  variants with missing values, indels, duplicates, multi-allelic
  sites, low MAF, monomorphic SNPs, nonsense values (p > 1,
  SE <= 0, etc.), failed liftOver positions, and high-MAF
  palindromic SNPs. Also harmonises alleles against an
  ancestry-matched reference panel (HRC/EUR, JPB/EAS, AGR/AFR,
  1KG/AMR, 1KG/SAS), applies an imputation quality filter
  (info >= 0.6 for dosage, >= 0.8 for hard calls), rescales
  betas and SEs by a factor of 2 for proxy-case cohorts
  (Marioni et al. 2018; Liu et al. 2017), and recalculates
  p-values with higher precision where they are reported as
  exactly 0. Logs filter counts to a per-cohort QC summary file.

qc_plots.R
  Generates a standard battery of visual QC plots for a single
  set of summary statistics: QQ plot, Manhattan plot, beta and
  SE histograms, beta and SE boxplots per chromosome, p-value
  and MAF histograms, a P-vs-Z consistency plot, imputation
  info histogram (if available), neff histogram, and
  effect-allele frequency scatter plots against one or two
  ancestry-matched reference panels. Output is a set of PNG
  files written to a specified QC directory.

qc_plots_job.R
  Command-line wrapper for qc_plots.R. Accepts sumstats file,
  output file name, QC output directory, and ancestry as
  positional arguments.


Meta-analysis
-------------

meta_analysis_f.R
  Orchestrates the full GWAS meta-analysis workflow using METAL
  (Willer et al. 2010, Bioinformatics). Runs three stages in
  sequence: (1) within-ancestry meta-analysis of case-control
  cohorts, (2) within-ancestry meta-analysis of proxy-case
  cohorts, and (3) a combined meta-analysis across all
  case-control, proxy, and multi-ancestry strata. Handles
  harmonisation of n_case / n_proxy_case column naming across
  input files, syncs files to and from scratch storage via
  rsync, and delegates the actual METAL call to run_metal.R.

meta_analysis_job.R
  Command-line wrapper for meta_analysis_f.R. Accepts a
  comma-separated list of input sumstats files, output
  directory, output file prefix, and meta-analysis scheme
  (STDERR or SAMPLESIZE) as positional arguments.

run_metal.R
  Low-level function that constructs a METAL script, executes
  it, reformats the METAL output to GWAS catalog column names,
  adds rsIDs from a 1000 Genomes AFR reference file, and
  optionally filters output SNPs by neff threshold and minimum
  number of contributing cohorts. Supports both SE-weighted
  (STDERR) and N-weighted (SAMPLESIZE) schemes.


Coordinate and allele harmonisation
-------------------------------------

snp_modifyBuild_local.R
  Local reimplementation of bigsnpr::snp_modifyBuild(). Lifts
  over SNP positions between genome builds (e.g. hg38 to hg19)
  using the UCSC liftOver binary and a local chain file, with an
  optional round-trip consistency check to flag variants whose
  position does not lift back to the original coordinate.
  Created because the original bigsnpr function was failing to
  download chain files in the cluster environment.

flip.R
  Corrects allele strand after liftOver. liftOver updates
  chromosomal positions but does not adjust alleles for regions
  where the reference strand orientation changes between builds.
  flip() identifies such regions via a range-join against a
  CrossMap-formatted chain file, then complements the alleles
  (A<->T, C<->G) of affected variants. Requires a chain file
  pre-processed with CrossMap viewchain.


Helper functions
----------------

neff_f.R
  Computes the effective sample size as neff = 4 / (1/n_case +
  1/n_control). For proxy-case cohorts the result is divided by
  4 to account for the reduced statistical information content
  of proxy phenotypes (Liu et al. 2017, doi:10.1038/ng.3766).

h2l_h2o.R
  Converts SNP-heritability estimates between the observed
  (0-1) scale and the liability scale, using the prevalence
  correction derived in Lee et al. 2011 (AJHG, Eq. 23).
  Provides two functions: h2o_to_h2l() and h2l_to_h2o().

r2l_r2o.R
  Converts polygenic score R^2 between the observed scale and
  the liability scale (and back), following the method of Lee
  et al. 2012 (Genetic Epidemiology, Eq. 15). Provides two
  functions: prs_r2obs_to_r2liab() and prs_r2liab_to_r2obs().

riskloc.R
  Clusters genome-wide significant variants into independent
  loci. Assigns a ±250 kb (configurable) window around each
  significant variant, merges overlapping windows within
  chromosomes, and returns one row per locus with the index
  variant (lowest p-value). Pins the APOE locus to rs429358
  (19:45411941:C:T) when the top p-value is exactly 0.


Plotting
--------

manhattan_f.R
  Produces a publication-ready Manhattan plot from GWAS catalog
  summary statistics using ggplot2. Optionally annotates
  genome-wide significant loci with gene names (pinning the
  APOE locus to its label regardless of the lead SNP), and
  supports a drop-shadow text style for labels. Thins
  sub-significant SNPs by random sampling to reduce plot size.

miami_f.R
  Produces a Miami plot (two mirrored Manhattan plots) for
  side-by-side comparison of two sets of summary statistics,
  using a shared genomic x-axis and symmetric y-axis limits.
  Optionally annotates significant loci with gene symbols.
  Effective sample sizes are shown in the facet labels.


Cluster job utilities
---------------------

save_job.R
  Helper to write a cluster job shell script (job.sh) from a
  multi-line R string. Strips leading whitespace from each line.
  Intended for generating SLURM or SGE submission scripts
  programmatically from within R.
