# =============================================================
#
#  ldsc_rg.R
#  Submit LDSC genetic correlation jobs and collect results.
#
#  Table of Contents:
#    00. Set-up
#    01. Meta-analysis rg: proxy vs. case-control (main)
#    02. Meta-analysis rg: female vs. male
#    03. Other-trait rg: each trait vs. PGC-ALZ3
#    04. Collect other-trait rg results
#
# =============================================================

#### 00. Set-up ####

library(data.table)
library(dplyr)
library(here)
library(openxlsx)

MUNGE_DIR  <- here("analysis", "ldsc", "ldsc_munge")
RG_DIR     <- here("analysis", "ldsc", "ldsc_rg")
LD_PREFIX  <- here("data", "reference_data", "UKBB.ALL.ldscore", "UKBB.")

ANCESTRIES <- c("eur", "afr", "eas", "sas", "amr")
STRATA     <- c("proxy", "case_control", "combined")

# Helper: map ancestry code to UKBB LD score file suffix
anc_to_ldsc <- function(anc) ifelse(anc == "sas", "CSA", toupper(anc))

# Helper: submit a single ldsc rg job if both input files exist
submit_rg <- function(sst_f, sst_r, pheno_out, out_dir, ld_suffix) {
  if (!file.exists(sst_f) || !file.exists(sst_r)) return(invisible(NULL))
  system(paste0(
    "sbatch ", here("src", "ldscRg.sh"), " ",
    "-f ", sst_f, " ",
    "-r ", sst_r, " ",
    "-e ", pheno_out, " ",
    "-o ", out_dir, " ",
    "-l ", LD_PREFIX, ld_suffix, ".rsid"
  ))
}


#### 01. Meta-analysis rg: proxy vs. case-control (main) ####

for (anc in ANCESTRIES) {
  submit_rg(
    sst_f     = here(MUNGE_DIR, "meta", "main", "proxy",
                     paste0("main_proxy_", anc, "_yeaapoe_munge.sumstats.gz")),
    sst_r     = here(MUNGE_DIR, "meta", "main", "case_control",
                     paste0("main_case_control_", anc, "_yeaapoe_munge.sumstats.gz")),
    pheno_out = paste0("main_proxy_main_casec_", anc, "_yeaapoe"),
    out_dir   = here(RG_DIR, "meta"),
    ld_suffix = anc_to_ldsc(anc)
  )
}


#### 02. Meta-analysis rg: female vs. male ####

for (anc in ANCESTRIES) {
  for (proxy_casec in STRATA) {
    submit_rg(
      sst_f     = here(MUNGE_DIR, "meta", "fema", proxy_casec,
                       paste0("fema_", proxy_casec, "_", anc, "_yeaapoe_munge.sumstats.gz")),
      sst_r     = here(MUNGE_DIR, "meta", "male", proxy_casec,
                       paste0("male_", proxy_casec, "_", anc, "_yeaapoe_munge.sumstats.gz")),
      pheno_out = paste0("fema_", proxy_casec, "_male_", proxy_casec, "_", anc, "_yeaapoe"),
      out_dir   = here(RG_DIR, "meta"),
      ld_suffix = anc_to_ldsc(anc)
    )
  }
}


#### 03. Other-trait rg: each trait vs. PGC-ALZ3 ####

phenos <- read.xlsx(
    here("analysis", "lava", "other_sumstats", "other_sumstats_filepaths.xlsx"),
    sheet = "2025_05_13"
  ) %>%
  filter(!is.na(raw_file), clean_file == "pass") %>%
  mutate(
    pheno_unique = gsub(" ", "_", paste(pheno_short, tolower(Author), Year, sep = "_"))
  ) %>%
  pull(pheno_unique)

phenos <- c(phenos, "longevity_deelen_2019", "amyloid_beta_jansen_2022")

for (pheno in phenos) {
  for (proxy_casec in c("case_control", "combined")) {
    submit_rg(
      sst_f     = here(MUNGE_DIR, "other", paste0(pheno, "_munge.sumstats.gz")),
      sst_r     = here(MUNGE_DIR, "meta", "main", proxy_casec,
                       paste0("main_", proxy_casec, "_eur_yeaapoe_munge.sumstats.gz")),
      pheno_out = paste0(pheno, "_main_", proxy_casec, "_eur_yeaapoe"),
      out_dir   = here(RG_DIR, "other"),
      ld_suffix = "EUR"
    )
  }
}


#### 04. Collect other-trait rg results ####

results_df <- lapply(phenos, function(pheno) {
  lapply(c("combined", "case_control"), function(proxy_casec) {
    fread(here(RG_DIR, "other",
               paste0(pheno, "_main_", proxy_casec, "_eur_yeaapoe_rg.table"))) %>%
      mutate(pheno1 = pheno, pheno2 = paste0("alzh_", proxy_casec))
  }) %>% bind_rows()
}) %>% bind_rows()

fwrite(results_df,
  file  = here(RG_DIR, "other.txt"),
  quote = FALSE, sep = "\t", na = NA
)
