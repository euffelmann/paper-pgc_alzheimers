#### 0. setup ####
library(here)
library(data.table)
library(dplyr)
library(ggplot2)
source(here("R/miami_f.R"))

#### 1. plot combined ####
combined <- NULL
combined[["male"]] <- fread(here(
  "data/sumstats/meta/male/combined",
  "male_combined_all_neff0.6_nsumstats1.txt.gz"
))

combined[["fema"]] <- fread(here(
  "data/sumstats/meta/fema/combined/",
  "fema_combined_all_neff0.6_nsumstats1.txt.gz"
)) 

flames_male <- fread(here("analysis/flames/male/combined/male_combined_all_neff0.6_nsumstats1/FLAMES_scores.pred")) %>%
  tidyr::extract(filename, into = c("chromosome", "start", "end"),
          regex = "locus_(\\d+):(\\d+)-(\\d+)",
          convert = TRUE) %>%
  select(chromosome, start, end, symbol)
flames_fema <- fread(here("analysis/flames/fema/combined/fema_combined_all_neff0.6_nsumstats1/FLAMES_scores.pred")) %>%
  tidyr::extract(filename, into = c("chromosome", "start", "end"),
                 regex = "locus_(\\d+):(\\d+)-(\\d+)",
                 convert = TRUE) %>%
  select(chromosome, start, end, symbol)

loci_male <- fread(here("analysis/risk_loci/male/male_combined_all_neff0.6_nsumstats1_loci.txt")) %>%
  left_join(flames_male, by = c("chromosome", "start", "end")) %>%
  select(index_variant_id, symbol)
loci_fema <- fread(here("analysis/risk_loci/fema/fema_combined_all_neff0.6_nsumstats1_loci.txt")) %>%
  left_join(flames_fema, by = c("chromosome", "start", "end")) %>%
  select(index_variant_id, symbol)

p_combined <- miami_f(
    sst_1 = combined[["fema"]],
    sst_2 = combined[["male"]],
    label_1 = "Female",
    label_2 = "Male",
    loci_1 = loci_fema,
    loci_2 = loci_male,
    annotate_loci = T,
    max_log_p_value = 30
  )

ggsave(
  plot = p_combined,
  filename = here("analysis/figures/sex_combined_miami.png"),
  device = "png",
  width = 10,
  height = 8
)


#### 2. plot case_control ####
casec <- NULL
casec[["male"]] <- fread(here(
  "data/sumstats/meta/male/case_control",
  "male_case_control_all_neff0.6_nsumstats1.txt.gz"
))

casec[["fema"]] <- fread(here(
  "data/sumstats/meta/fema/case_control/",
  "fema_case_control_all_neff0.6_nsumstats1.txt.gz"
)) 

flames_male <- fread(here("analysis/flames/male/case_control/male_case_control_all_neff0.6_nsumstats1/FLAMES_scores.pred")) %>%
  tidyr::extract(filename, into = c("chromosome", "start", "end"),
                 regex = "locus_(\\d+):(\\d+)-(\\d+)",
                 convert = TRUE) %>%
  select(chromosome, start, end, symbol)
flames_fema <- fread(here("analysis/flames/fema/case_control/fema_case_control_all_neff0.6_nsumstats1/FLAMES_scores.pred")) %>%
  tidyr::extract(filename, into = c("chromosome", "start", "end"),
                 regex = "locus_(\\d+):(\\d+)-(\\d+)",
                 convert = TRUE) %>%
  select(chromosome, start, end, symbol)

loci_male <- fread(here("analysis/risk_loci/male/male_case_control_all_neff0.6_nsumstats1_loci.txt")) %>%
  left_join(flames_male, by = c("chromosome", "start", "end")) %>%
  select(index_variant_id, symbol)
loci_fema <- fread(here("analysis/risk_loci/fema/fema_case_control_all_neff0.6_nsumstats1_loci.txt")) %>%
  left_join(flames_fema, by = c("chromosome", "start", "end")) %>%
  select(index_variant_id, symbol)

p_casec <- miami_f(
  sst_1 = casec[["fema"]],
  sst_2 = casec[["male"]],
  label_1 = "Female",
  label_2 = "Male",
  loci_1 = loci_fema,
  loci_2 = loci_male,
  annotate_loci = F,
  max_log_p_value = 30
)

ggsave(
  plot = p_casec,
  filename = here("analysis/figures/sex_casec_miami.png"),
  device = "png",
  width = 10,
  height = 8
)

