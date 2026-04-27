# =============================================================
#
#  sex_miami.R
#  Miami plots comparing male vs. female AD GWAS summary
#  statistics for the combined and case-control strata.
#
# =============================================================

#### 00. Set-up ####

library(here)
library(data.table)
library(dplyr)
library(ggplot2)
source(here("R", "miami_f.R"))

# Plot dimensions (inches)
FIG_W <- 10
FIG_H <- 8

# -------------------------------------------------------------
# Helper: load and parse a FLAMES scores file, extracting locus
# coordinates from the filename column
# -------------------------------------------------------------
load_flames <- function(path) {
  fread(path) %>%
    tidyr::extract(
      filename,
      into   = c("chromosome", "start", "end"),
      regex  = "locus_(\\d+):(\\d+)-(\\d+)",
      convert = TRUE
    ) %>%
    select(chromosome, start, end, symbol)
}

# -------------------------------------------------------------
# Helper: load sumstats, FLAMES annotations, and loci for one
# sex × stratum combination, then join loci to gene symbols
# -------------------------------------------------------------
load_stratum <- function(sex, stratum) {

  prefix <- paste0(sex, "_", stratum, "_all_neff0.6_nsumstats1")

  sst <- fread(here(
    "data", "sumstats", "meta", sex, stratum,
    paste0(prefix, ".txt.gz")
  ))

  flames <- load_flames(here(
    "analysis", "flames", sex, stratum, prefix, "FLAMES_scores.pred"
  ))

  loci <- fread(here(
    "analysis", "risk_loci", sex,
    paste0(prefix, "_loci.txt")
  )) %>%
    left_join(flames, by = c("chromosome", "start", "end")) %>%
    select(index_variant_id, symbol)

  list(sst = sst, loci = loci)
}

# -------------------------------------------------------------
# Helper: save a plot as both PDF and PNG
# -------------------------------------------------------------
save_plot <- function(p, path_stem) {
  ggsave(paste0(path_stem, ".png"), plot = p,
         device = "png",       width = FIG_W, height = FIG_H)
  ggsave(paste0(path_stem, ".pdf"), plot = p,
         device = cairo_pdf,   width = FIG_W, height = FIG_H)
}


#### 01. Combined stratum ####

male_combined <- load_stratum("male", "combined")
fema_combined <- load_stratum("fema", "combined")

p_combined <- miami_f(
  sst_1          = fema_combined$sst,
  sst_2          = male_combined$sst,
  label_1        = "Female",
  label_2        = "Male",
  loci_1         = fema_combined$loci,
  loci_2         = male_combined$loci,
  annotate_loci  = TRUE,
  max_log_p_value = 30
)

save_plot(p_combined, here("analysis", "figures", "sex_combined_miami"))


#### 02. Case-control stratum ####

male_casec <- load_stratum("male", "case_control")
fema_casec <- load_stratum("fema", "case_control")

p_casec <- miami_f(
  sst_1          = fema_casec$sst,
  sst_2          = male_casec$sst,
  label_1        = "Female",
  label_2        = "Male",
  loci_1         = fema_casec$loci,
  loci_2         = male_casec$loci,
  annotate_loci  = FALSE,
  max_log_p_value = 30
)

save_plot(p_casec, here("analysis", "figures", "sex_casec_miami"))
