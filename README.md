# PGC-ALZ3: Genome-wide association study of Alzheimer's disease

This repository contains all code used to process genotype data, run
the meta-analysis, and perform post-GWAS analyses for the PGC-ALZ3
multi-ancestry GWAS of Alzheimer's disease.

**Paper:** [URL will be added upon publication]  
**Summary statistics:** [URL will be added upon acceptance]

---

## Directory structure

```
.
├── analysis/       Post-GWAS analyses and figure generation
├── data/           Scripts to run the GWAS meta-analysis
├── data_raw/       Scripts to clean raw genotype and summary statistic files
├── R/              Custom R functions sourced by analysis scripts
└── src/            Non-R executables (PLINK, METAL, liftOver, LDSC wrappers)
```

---

## analysis/

Post-GWAS analyses. Each subdirectory corresponds to a distinct analysis.

| Subdirectory   | Description |
|----------------|-------------|
| `eQTL_coloc/`  | Colocalization of AD GWAS signals with eQTL data |
| `figures/`     | All manuscript and supplementary figures |
| `flames/`      | FLAMES gene prioritisation scores for genome-wide significant loci |
| `lava/`        | Local genetic correlation analysis (LAVA) between AD and other traits |
| `ldsc/`        | LDSC SNP heritability (h2) and genetic correlation (rg) results |
| `metasoft/`    | METASOFT RE2 heterogeneity analysis outputs |
| `pgs/`         | Polygenic score evaluation across cohorts and ancestry groups |
| `risk_loci/`   | Genome-wide significant locus definitions and gene annotations |

Scripts that generate these outputs:

| Script | Description |
|--------|-------------|
| `eadb_vs_main.R` | Scatter plots comparing EADB and PGC-ALZ3 p-values |
| `h2l_curve.R` | h2 on the liability scale as a function of disease prevalence |
| `h2l_r2l_bar.R` | h2l (LDSC vs. SBayesRC) and incremental PGS R2l bar plots |
| `lava.R` | Prepares LAVA input, submits jobs, and collects local rg results |
| `ldsc_munge.R` | Munges summary statistics for LDSC (cohort, meta, and other traits) |
| `ldsc_h2.R` | Submits LDSC h2 jobs and collects results |
| `ldsc_rg.R` | Submits LDSC rg jobs and collects results |
| `manhattan.R` | Manhattan plots for all ancestry × stratum combinations |
| `pgs_eval.R` | Evaluates PGS R2l and AUC across cohorts and models |
| `pgs_r2l.R` | Bar plots of incremental PGS R2l |
| `repro.R` | Reproducibility analysis: EADB and PGC-ALZ2 loci in PGC-ALZ3 |
| `sex_miami.R` | Miami plots comparing female vs. male meta-analysis results |

---

## data/

Scripts to prepare and run the GWAS meta-analysis, including summary
statistic formatting, QC, and meta-analysis across ancestries and strata.

| Script | Description |
|--------|-------------|
| `02_format_raw_sumstats.R` | Formats raw per-cohort summary statistics to a common GWAS catalog format |
| `03_process_raw_sumstats.R` | Applies QC filters to formatted summary statistics |

---

## data_raw/

Scripts to clean raw genotype data and prepare cohort-level summary
statistics prior to meta-analysis.

| Script | Description |
|--------|-------------|
| `01_process_raw_genotype.R` | Per-batch genotype QC (DemGene), merging by array, and PCA |

---

## R/

Custom R functions sourced by analysis scripts. See `R/readme.txt` for
descriptions of all functions.

| Script | Description |
|--------|-------------|
| `filter_sumstats.R` | Filter meta-analysed sumstats by neff and cohort count |
| `flip.R` | Correct allele strand after liftOver |
| `h2l_h2o.R` | Convert h2 between observed and liability scales |
| `manhattan_f.R` | Manhattan plot function |
| `meta_analysis_f.R` | Orchestrate METAL meta-analysis across strata |
| `meta_analysis_job.R` | Command-line wrapper for meta_analysis_f.R |
| `miami_f.R` | Miami plot function |
| `neff_f.R` | Compute effective sample size |
| `process_sumstats.R` | Per-cohort summary statistic QC pipeline |
| `qc_plots.R` | QC diagnostic plots for summary statistics |
| `qc_plots_job.R` | Command-line wrapper for qc_plots.R |
| `r2l_r2o.R` | Convert PGS R2 between observed and liability scales |
| `riskloc.R` | Cluster genome-wide significant variants into independent loci |
| `run_metal.R` | Low-level METAL execution and output formatting |
| `save_job.R` | Write a cluster job shell script from an R string |
| `snp_modifyBuild_local.R` | Lift over SNP positions between genome builds |

---

## src/

Non-R executables and cluster job submission scripts. See `src/00_readme.txt`
for download sources and version information.

| File | Description |
|------|-------------|
| `liftOver` | UCSC liftOver executable (macOS and Linux) |
| `plink1.9/` | PLINK v1.9 for genotype QC and meta-analysis preparation |
| `plink2/` | PLINK v2 for multi-allelic variant filtering |
| `METAL-2020-05-05/` | METAL for fixed-effects GWAS meta-analysis |
| `generic-metal/` | Older METAL release retained for reproducibility |
| `mvGWAMA.py` | Multivariate GWAS meta-analysis script |
| `ldscH2.sh` | SLURM wrapper for LDSC h2 estimation |
| `ldscMunge.sh` | SLURM wrapper for LDSC sumstat munging |
| `ldscRg.sh` | SLURM wrapper for LDSC genetic correlation estimation |

---

## Reproducibility

R package versions are recorded in `renv.lock`. To restore the exact
package environment used for these analyses:

```r
install.packages("renv")
renv::restore()
```

All file paths are managed with the `here` package and are relative to
the project root (`paper-pgc_alzheimers.Rproj`). No absolute paths are
used in any script.