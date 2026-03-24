# =============================================================
#
#  ldsc_munge.R
#  Prepare and munge summary statistics for LDSC analysis.
#  Munging is submitted as SLURM jobs via ldscMunge.sh.
#
#  For each sumstats file, a temporary version is first written
#  with n_case/n_control removed (these columns cause ldsc to
#  ignore --N-col) and p-values clamped to ≥ 1e-310. For the
#  APOE-excluded version (_nayapoe), a ±500 kb window around
#  the APOE e4 variant (chr19:45411941, GRCh37) is also masked.
#
#  Table of Contents:
#    00. Set-up
#    01. Munge per-cohort sumstats
#    02. Munge meta-analysis sumstats
#    03. Munge other-trait sumstats (for LAVA / LDSC rg)
#
# =============================================================

#### 00. Set-up ####

library(data.table)
library(dplyr)
library(here)
library(openxlsx)

MUNGE_DIR  <- here("analysis", "ldsc", "ldsc_munge")
TEMP_DIR   <- here(MUNGE_DIR, "temp_files")
APOE_CHR   <- 19
APOE_POS   <- 45411941
APOE_WIN   <- 5e5    # ±500 kb mask for _nayapoe version
P_MIN      <- 1e-310 # minimum representable p-value for ldsc
MUNGE_COLS <- "--N-col neff --a1 effect_allele --a2 other_allele --snp rsid"

# -------------------------------------------------------------
# Helper: write a temporary LDSC-ready sumstats file.
# Strips n_case/n_control, clamps p-values, and optionally
# masks the APOE locus (for _nayapoe analyses).
# -------------------------------------------------------------
write_temp_sst <- function(sst, outfile, mask_apoe = FALSE) {
  temp <- sst %>%
    select(-any_of(c("n_case", "n_control"))) %>%
    mutate(p_value = pmax(as.numeric(p_value), P_MIN))

  if (mask_apoe) {
    temp <- temp %>%
      filter(!(chromosome == APOE_CHR &
               base_pair_location > APOE_POS - APOE_WIN &
               base_pair_location < APOE_POS + APOE_WIN))
  }

  fwrite(temp,
    file  = file.path(TEMP_DIR, paste0(outfile, ".txt.gz")),
    sep   = "\t", na = "NA", quote = FALSE
  )
}


#### 01. Munge per-cohort sumstats ####

params <- readWorkbook(
    here("analysis", "cohort_summary", "cohort_summary.xlsx"),
    sheet = "main"
  ) %>%
  mutate(
    subdir        = ifelse(proxy_casecontrol == "proxy", "proxy_meta", proxy_casecontrol),
    sumstats_file = here(
      "data", "sumstats", subdir,
      paste0("clean_", cohort, "_", proxy_casecontrol, "_", analysis, "_", ancestry, ".txt.gz")
    )
  )

for (i in seq_len(nrow(params))) {
  for (apoe in c("_yeaapoe", "_nayapoe")) {

    outfile <- paste0(
      params$cohort[i], "_", params$proxy_casecontrol[i], "_",
      params$analysis[i], "_", params$ancestry[i], apoe
    )
    out_munged <- here(MUNGE_DIR, "coho", paste0(outfile, "_munge.sumstats.gz"))

    if (!file.exists(out_munged)) {
      message("Munging: ", outfile)
      write_temp_sst(
        sst        = fread(params$sumstats_file[i]),
        outfile    = outfile,
        mask_apoe  = (apoe == "_nayapoe")
      )
      system(paste0(
        "sbatch ", here("src", "ldscMunge.sh"), " ",
        "-f ", file.path(TEMP_DIR, paste0(outfile, ".txt.gz")), " ",
        "-e ", outfile, " ",
        "-o ", here(MUNGE_DIR, "coho"), " ",
        "-m '", MUNGE_COLS, "'"
      ))
    }
  }
}


#### 02. Munge meta-analysis sumstats ####

for (analysis in c("main", "fema", "male", "noag", "apoe", "agex")) {
  for (subdir in c("case_control", "proxy", "combined")) {
    for (ancestry in c("eur", "afr", "eas", "sas", "amr")) {

      # Only _yeaapoe version for meta sumstats
      outfile       <- paste0(analysis, "_", subdir, "_", ancestry, "_yeaapoe")
      sumstats_file <- here(
        "data", "sumstats", "meta", analysis, subdir,
        paste0(analysis, "_", subdir, "_", ancestry, ".txt.gz")
      )
      out_munged <- here(MUNGE_DIR, "meta", analysis, subdir,
                         paste0(outfile, "_munge.sumstats.gz"))

      if (file.exists(sumstats_file) && !file.exists(out_munged)) {
        message("Munging: ", outfile)
        write_temp_sst(sst = fread(sumstats_file), outfile = outfile, mask_apoe = FALSE)
        system(paste0(
          "sbatch ", here("src", "ldscMunge.sh"), " ",
          "-f ", file.path(TEMP_DIR, paste0(outfile, ".txt.gz")), " ",
          "-e ", outfile, " ",
          "-o ", here(MUNGE_DIR, "meta", analysis, subdir), " ",
          "-m '", MUNGE_COLS, "'"
        ))
      }
    }
  }
}


#### 03. Munge other-trait sumstats (for LAVA / LDSC rg) ####

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
  system(paste0(
    "sbatch ", here("src", "ldscMunge.sh"), " ",
    "-f ", here("analysis", "lava", "other_sumstats", "clean_sumstats",
                paste0(pheno, ".txt.gz")), " ",
    "-e ", pheno, " ",
    "-o ", here(MUNGE_DIR, "other"), " ",
    "-m '--N-col N --a1 A1 --a2 A2 --snp SNP'"
  ))
}
