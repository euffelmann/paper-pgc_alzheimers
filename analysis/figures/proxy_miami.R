#### 0. setup ####
library(here)
library(data.table)
library(dplyr)
library(ggplot2)
source(here("R/miami_f.R"))


#### 1. plot ####
sst <- NULL
sst[["proxy"]] <- fread(here(
  "data/sumstats/meta/main/proxy",
  "main_proxy_all_neff0.6_nsumstats1.txt.gz"
))

sst[["case_control"]] <- fread(here(
  "data/sumstats/meta/main/case_control",
  "main_case_control_all_neff0.6_nsumstats1.txt.gz"
))

p <- miami_f(
  sst_1 = sst[["case_control"]],
  sst_2 = sst[["proxy"]],
  label_1 = "Case-Control",
  label_2 = "Proxy",
  annotate_loci = F,
  max_log_p_value = 100
)

ggsave(
  plot = p,
  filename = here("analysis/figures/proxy_vs_casec_all_miami.png"),
  device = "png",
  width = 10,
  height = 8
)