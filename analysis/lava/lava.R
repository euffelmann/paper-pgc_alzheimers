# =============================================================
#
#  lava.R
#  Local Analysis of [co]Variant Association (LAVA) to estimate
#  local SNP heritability and genetic correlations between
#  Alzheimer's disease (PGC-ALZ3) and other complex traits.
#
#  Table of Contents:
#    00. Set-up
#    01. APOE locus — local h2 in PGC-ALZ3
#    02. Genetic correlations with other GWAS traits
#        (i)   Prepare and reformat per-trait sumstats
#        (ii)  Write LAVA info and locus files
#        (iii) Write sample overlap files (from LDSC rg)
#        (iv)  Run LAVA
#        (v)   Collect results
#
# =============================================================

#### 00. Set-up ####

library(here)
library(data.table)
library(dplyr)
library(openxlsx)
library(LAVA)
source(here("R", "save_job.R"))
source(here("R", "snp_modifyBuild_local.R"))

# g1000 EUR reference — used for rsID lookup when sumstats lack rsIDs
ref <- fread(here("data", "reference_data", "g1000_eur", "g1000_eur.bim")) %>%
  mutate(variant_id = paste(V1, V4, pmin(V5, V6), pmax(V5, V6), sep = ":"))

# Helper: save a cleaned sumstats data frame in LAVA format
save_lava_sst <- function(sst_clean, pheno) {
  fwrite(
    sst_clean,
    file  = here("analysis", "lava", "other_sumstats", "clean_sumstats",
                 paste0(pheno, ".txt.gz")),
    sep   = "\t",
    na    = NA,
    quote = FALSE
  )
}


#### 01. APOE locus — local h2 in PGC-ALZ3 ####

##### (i) Prepare LAVA-format sumstats #####

for (proxy_casec in c("combined", "case_control")) {

  out_file <- here(
    "analysis", "lava", "main_sumstats",
    paste0("main_", proxy_casec, "_eur_neff0.6_nsumstats1.txt.gz")
  )

  if (!file.exists(out_file)) {

    fread(here(
        "data", "sumstats", "meta", "main", proxy_casec,
        paste0("main_", proxy_casec, "_eur_neff0.6_nsumstats1.txt.gz")
      )) %>%
      mutate(Z = beta / standard_error) %>%
      rename(SNP = rsid, A1 = effect_allele, A2 = other_allele, N = neff) %>%
      select(SNP, A1, A2, N, Z) %>%
      fwrite(file = out_file, sep = "\t", na = NA, quote = FALSE)
  }
}

##### (ii) Write info and locus files #####

sst_casec <- fread(here(
  "analysis", "lava", "main_sumstats",
  "main_case_control_eur_neff0.6_nsumstats1.txt.gz"
))

# Info file: duplicate the phenotype to enable self-correlation as a h2 check
info_file <- data.frame(
  phenotype  = c("alzh", "copy"),
  cases      = rep(round(max(sst_casec$N) / 2), 2),
  controls   = rep(round(max(sst_casec$N) / 2), 2),
  prevalence = c(0.05, 0.05),
  filename   = rep(
    here("analysis", "lava", "main_sumstats",
         "main_case_control_eur_neff0.6_nsumstats1.txt.gz"),
    2
  )
)
fwrite(info_file,
  file  = here("analysis", "lava", "main_sumstats",
               "main_case_control_eur_neff0.6_nsumstats1.info.txt"),
  sep   = "\t", na = NA, quote = FALSE
)

# Locus file: ±2 Mb window around the APOE e4 variant (rs429358)
loc_file <- data.frame(
  LOC   = "1",
  CHR   = 19,
  START = 45411941 - 2e6,
  STOP  = 45411941 + 2e6
)
fwrite(loc_file,
  file  = here("analysis", "lava", "main_sumstats",
               "main_case_control_eur_neff0.6_nsumstats1.loc.txt"),
  sep   = "\t", na = NA, quote = FALSE
)

##### (iii) Run LAVA #####

loci  <- read.loci(here("analysis", "lava", "main_sumstats",
                        "main_case_control_eur_neff0.6_nsumstats1.loc.txt"))
input <- process.input(
  here("analysis", "lava", "main_sumstats",
       "main_case_control_eur_neff0.6_nsumstats1.info.txt"),
  sample.overlap.file = NULL,
  ref.prefix          = here("data", "reference_data",
                             "lava-ukb-v1.1_chr17-23", "lava-ukb-v1.1_chr19"),
  unlist(strsplit("alzh;copy", ";"))
)

locus <- process.locus(loci[1, ], input)
run.univ.bivar(locus = locus)$univ
# Expected output:
#   phen   h2.obs h2.latent ascertained p
# 1 alzh 0.104918 0.0890187        TRUE 0
# 2 copy 0.104918 0.0890187        TRUE 0


#### 02. Genetic correlations with other GWAS traits ####

##### (i) Prepare and reformat per-trait sumstats #####

df <- read.xlsx(
    here("analysis", "lava", "other_sumstats", "other_sumstats_filepaths.xlsx"),
    sheet = "2025_05_13"
  ) %>%
  filter(!is.na(raw_file)) %>%
  mutate(pheno_unique = paste(pheno_short, tolower(Author), Year, sep = "_"))

RAW_DIR <- here("analysis", "lava", "other_sumstats", "raw_sumstats")

###### 1. hdl_c_teslovich_2010 ######
pheno <- "hdl_c_teslovich_2010"
# N source: https://www.ebi.ac.uk/gwas/studies/GCST000755
fread(here(RAW_DIR, df$raw_file[df$pheno_unique == pheno])) %>%
  mutate(
    N             = 99900,
    effect_allele = toupper(effect_allele),
    other_allele  = toupper(other_allele),
    Z             = beta / standard_error
  ) %>%
  select(variant_id, effect_allele, other_allele, N, Z, p_value) %>%
  rename(SNP = variant_id, A1 = effect_allele, A2 = other_allele, P = p_value) %>%
  filter(N >= 0.6 * max(.$N)) %>%
  save_lava_sst(pheno)

###### 2. ldl_c_teslovich_2010 ######
pheno <- "ldl_c_teslovich_2010"
# N source: https://www.ebi.ac.uk/gwas/studies/GCST000755
fread(here(RAW_DIR, df$raw_file[df$pheno_unique == pheno])) %>%
  mutate(
    N             = 95454,
    effect_allele = toupper(effect_allele),
    other_allele  = toupper(other_allele),
    Z             = beta / standard_error
  ) %>%
  select(variant_id, effect_allele, other_allele, N, Z, p_value) %>%
  rename(SNP = variant_id, A1 = effect_allele, A2 = other_allele, P = p_value) %>%
  filter(N >= 0.6 * max(.$N)) %>%
  save_lava_sst(pheno)

###### 3. ms_sawcer_2011 ######
pheno <- "ms_sawcer_2011"
fread(here(RAW_DIR, df$raw_file[df$pheno_unique == pheno])) %>%
  mutate(
    neff = 4 / ((1 / 9772) + (1 / 16849)),
    Z    = beta / standard_error
  ) %>%
  select(variant_id, effect_allele, other_allele, neff, Z, p_value) %>%
  rename(SNP = variant_id, A1 = effect_allele, A2 = other_allele, N = neff, P = p_value) %>%
  filter(N >= 0.6 * max(.$N), !is.na(Z)) %>%
  save_lava_sst(pheno)

# Note: ldl_c_willer_2013, hdl_c_willer_2013 excluded — too few SNPs after QC.

###### 4. longevity_pilling_2016 (combined parental age at death) ######
pheno <- "longevity_pilling_2016"
# N source: https://www.ebi.ac.uk/gwas/studies/GCST003394
fread(here(RAW_DIR, df$raw_file[df$pheno_unique == pheno])) %>%
  mutate(N = 45627, Z = beta / standard_error) %>%
  select(variant_id, effect_allele, other_allele, N, Z, p_value) %>%
  rename(SNP = variant_id, A1 = effect_allele, A2 = other_allele, P = p_value) %>%
  filter(N >= 0.6 * max(.$N), !is.na(Z)) %>%
  save_lava_sst(pheno)

###### 5. longevity_extreme_pilling_2016 (parental extreme longevity ≥95 years) ######
pheno <- "longevity_extreme_pilling_2016"
# N source: https://www.ebi.ac.uk/gwas/studies/GCST003395
fread(here("analysis", "lava", "other_sumstats", "raw_sumstats",
           "27015805-GCST003395-EFO_0007796-build37.f.tsv.gz")) %>%
  mutate(
    N = 4 / ((1 / 1339) + (1 / 40934)),
    Z = beta / standard_error
  ) %>%
  select(variant_id, effect_allele, other_allele, N, Z, p_value) %>%
  rename(SNP = variant_id, A1 = effect_allele, A2 = other_allele, P = p_value) %>%
  filter(N >= 0.6 * max(.$N), !is.na(Z)) %>%
  save_lava_sst(pheno)

###### 6. ea_davies_2016 ######
pheno <- "ea_davies_2016"
fread(here(RAW_DIR, df$raw_file[df$pheno_unique == pheno])) %>%
  mutate(N = 112151) %>%
  select(variant_id, effect_allele, other_allele, N, beta, p_value) %>%
  rename(SNP = variant_id, A1 = effect_allele, A2 = other_allele, P = p_value, B = beta) %>%
  filter(N >= 0.6 * max(.$N), !is.na(B)) %>%
  save_lava_sst(pheno)

###### 7. cad_van_der_harst_2017 ######
pheno <- "cad_van_der_harst_2017"
fread(here(RAW_DIR, df$raw_file[df$pheno_unique == "cad_van der harst_2017"][1])) %>%
  mutate(
    neff          = round(4 / ((1 / 122733) + (1 / 424528))),
    effect_allele = toupper(effect_allele),
    other_allele  = toupper(other_allele),
    Z             = beta / standard_error
  ) %>%
  select(variant_id, effect_allele, other_allele, neff, Z, p_value) %>%
  rename(SNP = variant_id, A1 = effect_allele, A2 = other_allele, P = p_value, N = neff) %>%
  filter(N >= 0.6 * max(.$N), !is.na(Z)) %>%
  save_lava_sst(pheno)

# Note: ms_beecham_2013 excluded — too few SNPs after QC.

###### 8. int_savage_2018 ######
pheno <- "int_savage_2018"
fread(here(RAW_DIR, df$raw_file[df$pheno_unique == pheno])) %>%
  mutate(
    effect_allele = toupper(effect_allele),
    other_allele  = toupper(other_allele)
  ) %>%
  select(variant_id, effect_allele, other_allele, n_analyzed, z, p_value) %>%
  rename(SNP = variant_id, A1 = effect_allele, A2 = other_allele,
         P = p_value, N = n_analyzed, Z = z) %>%
  filter(N >= 0.6 * max(.$N), !is.na(Z)) %>%
  save_lava_sst(pheno)

###### 9. cog_lee_2018 ######
pheno <- "cog_lee_2018"
fread(here(RAW_DIR, df$raw_file[df$pheno_unique == pheno])) %>%
  mutate(N = 257841, Z = beta / standard_error) %>%
  select(variant_id, effect_allele, other_allele, N, Z, p_value) %>%
  rename(SNP = variant_id, A1 = effect_allele, A2 = other_allele, P = p_value) %>%
  filter(N >= 0.6 * max(.$N), !is.na(Z)) %>%
  save_lava_sst(pheno)

# Note: trem2_liu_2019 excluded — effect allele column missing in source file.

###### 10. pd_nalls_2019 ######
pheno <- "pd_nalls_2019"
fread(here(RAW_DIR, df$raw_file[df$pheno_unique == pheno])) %>%
  mutate(
    neff = round(4 / ((1 / N_cases) + (1 / N_controls))),
    Z    = beta / standard_error
  ) %>%
  select(rsid, effect_allele, other_allele, neff, Z, p_value) %>%
  rename(SNP = rsid, A1 = effect_allele, A2 = other_allele, P = p_value, N = neff) %>%
  filter(N >= 0.6 * max(.$N), !is.na(Z)) %>%
  save_lava_sst(pheno)

# Note: cortical_thickness_van_der_meer_2020 excluded — no beta column.
# Note: pd_tan_2020 excluded — sample size too small (N = 2755).

###### 11. ldl_c_klimentidis_2020 ######
pheno <- "ldl_c_klimentidis_2020"
fread(here(RAW_DIR, df$raw_file[df$pheno_unique == pheno])) %>%
  mutate(
    variant_id = paste(chromosome, base_pair_location,
                       pmin(effect_allele, other_allele),
                       pmax(effect_allele, other_allele), sep = ":"),
    N = 431167,
    Z = beta / standard_error
  ) %>%
  left_join(ref, by = "variant_id") %>%
  arrange(chromosome, base_pair_location) %>%
  select(V2, effect_allele, other_allele, N, Z, p_value) %>%
  rename(SNP = V2, A1 = effect_allele, A2 = other_allele, P = p_value) %>%
  filter(N >= 0.6 * max(.$N), !is.na(Z), !is.na(SNP)) %>%
  save_lava_sst(pheno)

###### 12. triglyceride_barton_2021 ######
pheno <- "triglyceride_barton_2021"
fread(here(RAW_DIR, df$raw_file[df$pheno_unique == pheno])) %>%
  mutate(
    variant_id = paste(chromosome, base_pair_location,
                       pmin(ALLELE1, ALLELE0), pmax(ALLELE1, ALLELE0), sep = ":"),
    N = 437532,
    Z = beta / standard_error,
    P = as.numeric(p_value)
  ) %>%
  left_join(ref, by = "variant_id") %>%
  arrange(chromosome, base_pair_location) %>%
  select(V2, ALLELE1, ALLELE0, N, Z, P) %>%
  rename(SNP = V2, A1 = ALLELE1, A2 = ALLELE0) %>%
  filter(N >= 0.6 * max(.$N), !is.na(Z), !is.na(SNP)) %>%
  save_lava_sst(pheno)

###### 13. als_van_rheenen_2021 ######
pheno <- "als_van_rheenen_2021"
fread(here(RAW_DIR, df$raw_file[df$pheno_unique == "als_van rheenen_2021"])) %>%
  mutate(
    effect_allele = toupper(effect_allele),
    other_allele  = toupper(other_allele),
    Z             = beta / standard_error
  ) %>%
  select(rsid, effect_allele, other_allele, N_effective, Z, p_value) %>%
  rename(SNP = rsid, A1 = effect_allele, A2 = other_allele, N = N_effective, P = p_value) %>%
  filter(N >= 0.6 * max(.$N), !is.na(Z)) %>%
  save_lava_sst(pheno)

###### 14. apo_a1_richardson_2022 ######
pheno <- "apo_a1_richardson_2022"
fread(here(RAW_DIR, df$raw_file[df$pheno_unique == pheno])) %>%
  mutate(N = 115082, Z = beta / standard_error) %>%
  select(variant_id, effect_allele, other_allele, N, Z, p_value) %>%
  rename(SNP = variant_id, A1 = effect_allele, A2 = other_allele, P = p_value) %>%
  filter(N >= 0.6 * max(.$N), !is.na(Z)) %>%
  save_lava_sst(pheno)

###### 15. apo_b_richardson_2022 ######
pheno <- "apo_b_richardson_2022"
fread(here(RAW_DIR, df$raw_file[df$pheno_unique == pheno])) %>%
  mutate(N = 115082, Z = beta / standard_error) %>%
  select(variant_id, effect_allele, other_allele, N, Z, p_value) %>%
  rename(SNP = variant_id, A1 = effect_allele, A2 = other_allele, P = p_value) %>%
  filter(N >= 0.6 * max(.$N), !is.na(Z)) %>%
  save_lava_sst(pheno)

###### 16. total_c_richardson_2022 ######
pheno <- "total_c_richardson_2022"
fread(here(RAW_DIR, df$raw_file[df$pheno_unique == pheno])) %>%
  mutate(N = 115082, Z = beta / standard_error) %>%
  select(variant_id, effect_allele, other_allele, N, Z, p_value) %>%
  rename(SNP = variant_id, A1 = effect_allele, A2 = other_allele, P = p_value) %>%
  filter(N >= 0.6 * max(.$N), !is.na(Z)) %>%
  save_lava_sst(pheno)

###### 17. p_tau_jansen_2022 ######
# Requires liftOver from hg38 to hg19
pheno <- "p_tau_jansen_2022"
fread(here(RAW_DIR, df$raw_file[df$pheno_unique == pheno])) %>%
  rename(chr = chromosome, pos = base_pair_location) %>%
  snp_modifyBuild_local(liftOver = here("src", "liftOver"),
                        from = "hg38", to = "hg19", check_reverse = TRUE) %>%
  mutate(
    effect_allele = toupper(effect_allele),
    other_allele  = toupper(other_allele),
    variant_id    = paste(chr, pos,
                          pmin(effect_allele, other_allele),
                          pmax(effect_allele, other_allele), sep = ":")
  ) %>%
  left_join(ref, by = "variant_id") %>%
  select(V2, effect_allele, other_allele, N, Zscore, p_value) %>%
  rename(SNP = V2, A1 = effect_allele, A2 = other_allele, P = p_value, Z = Zscore) %>%
  filter(N >= 0.6 * max(.$N), !is.na(Z), !is.na(SNP)) %>%
  save_lava_sst(pheno)

###### 18. cad_aragam_2022 ######
pheno <- "cad_aragam_2022"
fread(here(RAW_DIR, df$raw_file[df$pheno_unique == pheno])) %>%
  mutate(
    effect_allele = toupper(effect_allele),
    other_allele  = toupper(other_allele),
    variant_id    = paste(chromosome, base_pair_location,
                          pmin(effect_allele, other_allele),
                          pmax(effect_allele, other_allele), sep = ":"),
    controls = n - cases,
    neff      = round(4 / ((1 / cases) + (1 / controls))),
    Z         = beta / standard_error
  ) %>%
  left_join(ref, by = "variant_id") %>%
  select(V2, effect_allele, other_allele, neff, Z, p_value) %>%
  rename(SNP = V2, A1 = effect_allele, A2 = other_allele, P = p_value, N = neff) %>%
  filter(N >= 0.6 * max(.$N), !is.na(Z), !is.na(SNP)) %>%
  save_lava_sst(pheno)

###### 19. hdl_c_graham_2021 ######
# Requires liftOver from hg38 to hg19
pheno <- "hdl_c_graham_2021"
fread(here(RAW_DIR, df$raw_file[df$pheno_unique == pheno])) %>%
  rename(chr = chromosome, pos = base_pair_location) %>%
  snp_modifyBuild_local(liftOver = here("src", "liftOver"),
                        from = "hg38", to = "hg19", check_reverse = TRUE) %>%
  mutate(
    effect_allele = toupper(effect_allele),
    other_allele  = toupper(other_allele),
    variant_id    = paste(chr, pos,
                          pmin(effect_allele, other_allele),
                          pmax(effect_allele, other_allele), sep = ":"),
    p_value       = as.numeric(p_value),
    Z             = beta / standard_error
  ) %>%
  left_join(ref, by = "variant_id") %>%
  select(V2, effect_allele, other_allele, n, Z, p_value) %>%
  rename(SNP = V2, A1 = effect_allele, A2 = other_allele, P = p_value, N = n) %>%
  filter(N >= 0.6 * max(.$N), !is.na(Z), !is.na(SNP)) %>%
  save_lava_sst(pheno)

###### 20. ldl_c_graham_2021 ######
# Requires liftOver from hg38 to hg19
pheno <- "ldl_c_graham_2021"
fread(here(RAW_DIR, df$raw_file[df$pheno_unique == pheno])) %>%
  filter(n >= 0.6 * max(.$n)) %>%
  mutate(p_value = as.numeric(p_value)) %>%
  rename(chr = chromosome, pos = base_pair_location) %>%
  snp_modifyBuild_local(liftOver = here("src", "liftOver"),
                        from = "hg38", to = "hg19", check_reverse = TRUE) %>%
  mutate(
    effect_allele = toupper(effect_allele),
    other_allele  = toupper(other_allele),
    variant_id    = paste(chr, pos,
                          pmin(effect_allele, other_allele),
                          pmax(effect_allele, other_allele), sep = ":"),
    Z = beta / standard_error
  ) %>%
  left_join(ref, by = "variant_id") %>%
  select(V2, effect_allele, other_allele, n, Z, p_value) %>%
  rename(SNP = V2, A1 = effect_allele, A2 = other_allele, P = p_value, N = n) %>%
  filter(!is.na(Z), !is.na(SNP)) %>%
  save_lava_sst(pheno)

###### 21. triglyceride_graham_2021 ######
# Requires liftOver from hg38 to hg19
pheno <- "triglyceride_graham_2021"
fread(here(RAW_DIR, df$raw_file[df$pheno_unique == pheno])) %>%
  filter(n >= 0.6 * max(.$n)) %>%
  rename(chr = chromosome, pos = base_pair_location) %>%
  snp_modifyBuild_local(liftOver = here("src", "liftOver"),
                        from = "hg38", to = "hg19", check_reverse = TRUE) %>%
  mutate(
    effect_allele = toupper(effect_allele),
    other_allele  = toupper(other_allele),
    variant_id    = paste(chr, pos,
                          pmin(effect_allele, other_allele),
                          pmax(effect_allele, other_allele), sep = ":"),
    Z = beta / standard_error,
    P = as.numeric(p_value)
  ) %>%
  left_join(ref, by = "variant_id") %>%
  select(V2, effect_allele, other_allele, n, Z, P) %>%
  rename(SNP = V2, A1 = effect_allele, A2 = other_allele, N = n) %>%
  filter(!is.na(Z), !is.na(SNP)) %>%
  save_lava_sst(pheno)

###### 22. total_c_graham_2021 ######
# Requires liftOver from hg38 to hg19
pheno <- "total_c_graham_2021"
fread(here(RAW_DIR, df$raw_file[df$pheno_unique == pheno])) %>%
  filter(n >= 0.6 * max(.$n)) %>%
  rename(chr = chromosome, pos = base_pair_location) %>%
  snp_modifyBuild_local(liftOver = here("src", "liftOver"),
                        from = "hg38", to = "hg19", check_reverse = TRUE) %>%
  mutate(
    effect_allele = toupper(effect_allele),
    other_allele  = toupper(other_allele),
    variant_id    = paste(chr, pos,
                          pmin(effect_allele, other_allele),
                          pmax(effect_allele, other_allele), sep = ":"),
    p_value = as.numeric(p_value),
    Z       = beta / standard_error
  ) %>%
  left_join(ref, by = "variant_id") %>%
  select(V2, effect_allele, other_allele, n, Z, p_value) %>%
  rename(SNP = V2, A1 = effect_allele, A2 = other_allele, P = p_value, N = n) %>%
  filter(!is.na(Z), !is.na(SNP)) %>%
  save_lava_sst(pheno)


##### (ii) Write LAVA info and locus files #####

cohort_summary <- read.xlsx(here("analysis", "cohort_summary", "cohort_summary.xlsx")) %>%
  filter(ancestry == "eur")

df <- read.xlsx(
    here("analysis", "lava", "other_sumstats", "other_sumstats_filepaths.xlsx"),
    sheet = "2025_05_13"
  ) %>%
  filter(!is.na(raw_file), clean_file == "pass") %>%
  mutate(
    pheno_unique = paste(pheno_short, tolower(Author), Year, sep = "_"),
    pheno_unique = gsub(" ", "_", pheno_unique)
  )

###### a. Info files ######
for (i in seq_len(nrow(df))) {

  pheno     <- df$pheno_unique[i]
  n_case    <- df$n_case[i]
  n_control <- df$n_control[i]

  for (proxy_casec in c("combined", "case_control")) {

    alz_neff <- if (proxy_casec == "combined") {
      sum(cohort_summary$neff)
    } else {
      sum(cohort_summary$neff[cohort_summary$proxy_casecontrol == "case_control"])
    }

    data.frame(
      phenotype = c(pheno, paste0("alzh_", proxy_casec)),
      cases     = c(n_case,             round(alz_neff / 2)),
      controls  = c(n_control,          round(alz_neff / 2)),
      filename  = c(
        here("analysis", "lava", "other_sumstats", "clean_sumstats",
             paste0(pheno, ".txt.gz")),
        here("analysis", "lava", "main_sumstats",
             paste0("main_", proxy_casec, "_eur_neff0.6_nsumstats1.txt.gz"))
      )
    ) %>%
      fwrite(
        file  = here("analysis", "lava", "other_sumstats", "clean_sumstats",
                     paste0(pheno, "_", proxy_casec, ".info.txt")),
        sep   = "\t", na = NA, quote = FALSE
      )
  }
}

###### b. Locus files ######
alz_loci <- readRDS(here("analysis", "risk_loci", "risk_loci_table.rds"))[["unique_loci_main_combined"]] %>%
  select(locus, chromosome, start, end)

for (i in seq_len(nrow(df))) {

  trait_loci <- unlist(stringr::str_split(df$Loci[i], pattern = ";"))

  alz_loci %>%
    filter(locus %in% trait_loci) %>%
    transmute(LOC = locus, CHR = chromosome, START = start, STOP = end) %>%
    fwrite(
      file  = here("analysis", "lava", "other_sumstats", "clean_sumstats",
                   paste0(df$pheno_unique[i], ".loc.txt")),
      sep   = "\t", na = NA, quote = FALSE
    )
}


##### (iii) Write sample overlap files (from LDSC rg) #####

for (pheno in df$pheno_unique) {
  for (proxy_casec in c("combined", "case_control")) {

    alz_pheno <- paste0("alzh_", proxy_casec)
    gcov_int  <- fread(here(
        "analysis", "ldsc", "ldsc_rg", "other",
        paste0(pheno, "_main_", proxy_casec, "_eur_yeaapoe_rg.table")
      ))$gcov_int

    # Write symmetric 2×2 sample overlap matrix
    phenotypes <- c(pheno, alz_pheno)
    outfile    <- here("analysis", "lava", "other_sumstats", "clean_sumstats",
                       paste0(pheno, "_", proxy_casec, ".sample.overlap.txt"))
    con <- file(outfile, "w")
    cat(phenotypes,          sep = " ", file = con); cat("\n", file = con)
    cat(pheno,    1,         gcov_int,  sep = " ", file = con); cat("\n", file = con)
    cat(alz_pheno, gcov_int, 1,         sep = " ", file = con); cat("\n", file = con)
    close(con)
  }
}


##### (iv) Run LAVA #####

for (pheno in df$pheno_unique) {

  results_missing <- !file.exists(here(
    "analysis", "lava", "other_sumstats", "results",
    paste0(pheno, c("_combined.univ.txt", "_combined.biva.txt",
                    "_case_control.univ.txt", "_case_control.biva.txt"))
  ))

  if (!any(results_missing)) {
    message("Skipping (all outputs exist): ", pheno)
    next
  }

  message("Running LAVA for: ", pheno)

  loci  <- read.loci(here("analysis", "lava", "other_sumstats", "clean_sumstats",
                          paste0(pheno, ".loc.txt")))
  n.loc <- nrow(loci)

  for (proxy_casec in c("combined", "case_control")) {

    input <- process.input(
      here("analysis", "lava", "other_sumstats", "clean_sumstats",
           paste0(pheno, "_", proxy_casec, ".info.txt")),
      sample.overlap.file = here("analysis", "lava", "other_sumstats", "clean_sumstats",
                                 paste0(pheno, "_", proxy_casec, ".sample.overlap.txt")),
      here("data", "reference_data", "g1000_eur", "g1000_eur")
    )

    message("  Processing ", n.loc, " loci (", proxy_casec, ")")
    progress <- ceiling(quantile(seq_len(n.loc), seq(0.05, 1, 0.05)))

    univ_results <- vector("list", n.loc)
    biva_results <- vector("list", n.loc)

    for (i in seq_len(n.loc)) {

      if (i %in% progress) message("  ..", names(progress[which(progress == i)]))

      locus <- process.locus(loci[i, ], input)

      if (!is.null(locus)) {
        loc_info <- data.frame(
          locus  = locus$id,
          chr    = locus$chr,
          start  = locus$start,
          stop   = locus$stop,
          n.snps = locus$n.snps,
          n.pcs  = locus$K
        )
        loc_out            <- run.univ.bivar(locus, univ.thresh = 1e-4)
        univ_results[[i]]  <- cbind(loc_info, loc_out$univ)
        if (!is.null(loc_out$bivar))
          biva_results[[i]] <- cbind(loc_info, loc_out$bivar)
      }
    }

    fwrite(rbindlist(univ_results),
      file  = here("analysis", "lava", "other_sumstats", "results",
                   paste0(pheno, "_", proxy_casec, ".univ.txt")),
      sep   = "\t", na = NA, quote = FALSE
    )
    fwrite(rbindlist(biva_results),
      file  = here("analysis", "lava", "other_sumstats", "results",
                   paste0(pheno, "_", proxy_casec, ".biva.txt")),
      sep   = "\t", na = NA, quote = FALSE
    )
  }
}


##### (v) Collect results #####

results_df <- lapply(df$pheno_unique, function(pheno) {
  lapply(c("combined", "case_control"), function(proxy_casec) {
    fread(here("analysis", "lava", "other_sumstats", "results",
               paste0(pheno, "_", proxy_casec, ".biva.txt")))
  }) %>% rbindlist()
}) %>% rbindlist()

fwrite(
  results_df,
  file  = here("analysis", "lava", "other_sumstats", "results.txt"),
  sep   = "\t", na = NA, quote = FALSE
)
