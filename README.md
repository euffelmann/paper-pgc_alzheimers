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
├── analysis/   
├── data/       
├── data_raw/   
├── R/          
└── src/        
```

---

## analysis/

Post-GWAS analyses. Each subdirectory corresponds to a distinct analysis.

---

## data/

Scripts to prepare and run the GWAS meta-analysis, including summary
statistic formatting, QC, and meta-analysis across ancestries and phenotype
definition.

---

## data_raw/

Scripts to clean raw genotype data and prepare cohort-level summary
statistics prior to meta-analysis.

---

## R/

Custom R functions sourced by analysis scripts. See `R/readme.txt` for
descriptions of all functions.

---

## src/

Non-R executables and cluster job submission scripts. See `src/00_readme.txt`
for download sources and version information.

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