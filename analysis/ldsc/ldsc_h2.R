# =============================================================
#
#  ldsc_h2.R
#  Submit LDSC heritability jobs and collect results.
#
#  Table of Contents:
#    00. Set-up
#    01. Submit h2 jobs for meta-analysis sumstats
#    02. Collect results into a summary table
#    03. Submit h2 jobs with --chisq-max 999999 (sensitivity)
#
# =============================================================

#### 00. Set-up ####

library(data.table)
library(dplyr)
library(here)
library(openxlsx)

K          <- 0.05   # population prevalence of AD
H2_DIR     <- here("analysis", "ldsc", "ldsc_h2")
MUNGE_DIR  <- here("analysis", "ldsc", "ldsc_munge")
LD_PREFIX  <- here("data", "reference_data", "UKBB.ALL.ldscore", "UKBB.")

ANALYSES   <- c("main", "fema", "male", "noag", "apoe", "agex")
ANCESTRIES <- c("eur", "afr", "eas", "sas", "amr")
STRATA     <- c("proxy", "case_control", "combined")

# -------------------------------------------------------------
# Helper: map ancestry code to UKBB LD score file suffix
# (SAS uses "CSA" in the UKBB LD score file naming)
# -------------------------------------------------------------
anc_to_ldsc <- function(anc) ifelse(anc == "sas", "CSA", toupper(anc))

# -------------------------------------------------------------
# Helper: extract a single numeric value from an LDSC log file
# using a grep pattern and awk field index. Parentheses are
# stripped for SE fields reported as "(0.123)".
# -------------------------------------------------------------
extract_log_field <- function(log_file, pattern, field) {
  raw <- system(
    paste0("grep '", pattern, "' ", log_file, " | awk '{print $", field, "}'"),
    intern = TRUE
  )
  as.numeric(gsub("[()]", "", raw))
}


#### 01. Submit h2 jobs for meta-analysis sumstats ####

for (analysis in ANALYSES) {
  for (anc in ANCESTRIES) {
    for (proxy_casec in STRATA) {

      pheno      <- paste(analysis, proxy_casec, anc, "yeaapoe", sep = "_")
      sst_munged <- here(MUNGE_DIR, "meta", analysis, proxy_casec,
                         paste0(pheno, "_munge.sumstats.gz"))
      log_out    <- here(H2_DIR, "meta", paste0(pheno, "_h2.log"))

      if (file.exists(sst_munged) && !file.exists(log_out)) {
        system(paste0(
          "sbatch ", here("src", "ldscH2.sh"), " ",
          "-f ", sst_munged, " ",
          "-e ", pheno, " ",
          "-o ", here(H2_DIR, "meta"), " ",
          # samp-prev = 0.5 because neff is used as input to ldsc
          "-h '--pop-prev ", K, " --samp-prev 0.5' ",
          "-l ", LD_PREFIX, anc_to_ldsc(anc), ".rsid"
        ))
      }
    }
  }
}


#### 02. Collect results into a summary table ####

params <- readWorkbook(
    here("analysis", "cohort_summary", "cohort_summary.xlsx"),
    sheet = "main"
  ) %>%
  rename(proxy_casec = proxy_casecontrol, anc = ancestry)

h2_df <- lapply(ANALYSES, function(analysis) {
  lapply(ANCESTRIES, function(anc) {
    lapply(STRATA, function(proxy_casec) {

      log_file <- here(H2_DIR, "meta",
                       paste0(analysis, "_", proxy_casec, "_", anc, "_yeaapoe_h2.log"))

      if (!file.exists(log_file)) return(NULL)

      neff <- if (proxy_casec == "combined") {
        sum(params$neff[params$anc == anc])
      } else {
        sum(params$neff[params$anc == anc & params$proxy_casec == proxy_casec])
      }

      data.frame(
        proxy_casec    = proxy_casec,
        analysis       = analysis,
        anc            = anc,
        apoe           = "yeaapoe",
        K              = K,
        neff           = neff,
        h2             = extract_log_field(log_file, "Total Liability scale h2", 5),
        h2_se          = extract_log_field(log_file, "Total Liability scale h2", 6),
        lambda_gc      = extract_log_field(log_file, "Lambda GC",                3),
        mean_chi2      = extract_log_field(log_file, "Mean Chi\\^2",             3),
        intercept      = extract_log_field(log_file, "Intercept",                2),
        intercept_se   = extract_log_field(log_file, "Intercept",                3),
        ratio          = extract_log_field(log_file, "Ratio",                    2),
        ratio_se       = extract_log_field(log_file, "Ratio",                    3)
      )
    }) %>% bind_rows()
  }) %>% bind_rows()
}) %>% bind_rows()

h2_df %>%
  mutate(
    analysis    = factor(analysis,    levels = c("main", "noag", "agex", "fema", "male")),
    proxy_casec = factor(proxy_casec, levels = c("combined", "case_control", "proxy")),
    anc         = factor(anc,         levels = ANCESTRIES),
    p           = 2 * pnorm(abs(h2 / h2_se), lower.tail = FALSE)
  ) %>%
  relocate(p, .after = h2_se) %>%
  arrange(analysis, proxy_casec, anc, apoe) %>%
  fwrite(here("analysis", "ldsc", "ldsc_h2.txt"),
         quote = FALSE, sep = "\t", na = NA)


#### 03. Submit h2 jobs with --chisq-max 999999 (sensitivity) ####

# Runs only the main/combined/EUR stratum to assess the effect of the
# default --chisq-max filter on h2 estimates.

for (analysis in "main") {
  for (anc in "eur") {
    for (proxy_casec in "combined") {

      pheno      <- paste(analysis, proxy_casec, anc, "yeaapoe", sep = "_")
      sst_munged <- here(MUNGE_DIR, "meta", analysis, proxy_casec,
                         paste0(pheno, "_munge.sumstats.gz"))
      log_out    <- here(H2_DIR, "chi_sq_max", paste0(pheno, "_h2.log"))

      if (file.exists(sst_munged) && !file.exists(log_out)) {
        system(paste0(
          "sbatch ", here("src", "ldscH2.sh"), " ",
          "-f ", sst_munged, " ",
          "-e ", pheno, " ",
          "-o ", here(H2_DIR, "chi_sq_max"), " ",
          "-h '--pop-prev ", K, " --samp-prev 0.5 --chisq-max 999999' ",
          "-l ", LD_PREFIX, anc_to_ldsc(anc), ".rsid"
        ))
      }
    }
  }
}
