This folder contains the pipeline used to run SBayesRC for SNP heritability estimation on the liability scale.

The pipeline is documented in sbayesrc_pipeline.txt and consists of three steps:
  1. Building a block LD matrix from UK Biobank EUR-ancestry genotypes using GCTB.
  2. Converting logistic-regression summary statistics to the liability scale and imputing
     missing SNPs relative to the LD reference.
  3. Running SBayesRC with annotation-informed priors, treating the APOE lead variant
     (rs429358) as a fixed effect and SNPs in LD with APOE as random effects.

Paths in the pipeline script are illustrative and must be adjusted for your own system.
The analysis was run on EUR ancestry summary statistics
(main_case_control_eur_neff0.6_nsumstats1.txt) using a UKB EUR LD reference.

Software required:
  - GCTB (https://cnsgenomics.com/software/gctb/) for LD matrix construction and SBayesRC
  - R with data.table for the liability-scale conversion step (get_liability_ma_refAF.R)

For questions, please contact the corresponding authors.
