#### 0. set-up ####
library(here)
library(data.table)
library(dplyr)
library(purrr)
source(here("R/save_job.R"))

#### 1. select sumstats to be analyzed ####
params <- data.frame(
  sumstats_file = c(
    list.files(here("data_raw", "sumstats", "case_control"), pattern = "\\.gz$", full.names = T),
    list.files(here("data_raw", "sumstats", "proxy"), pattern = "\\.gz$", full.names = T),
    list.files(here("data_raw", "sumstats", "combined"), pattern = "\\.gz$", full.names = T)
  ),
  outfile = c(
    sub("\\.txt\\.gz$", "", list.files(here("data_raw", "sumstats", "case_control"), pattern = "\\.gz$")),
    sub("\\.txt\\.gz$", "", list.files(here("data_raw", "sumstats", "proxy"), pattern = "\\.gz$")),
    sub("\\.txt\\.gz$", "", list.files(here("data_raw", "sumstats", "combined"), pattern = "\\.gz$"))
  ),
  qc_dir = here("data_raw", "sumstats", "qc_plots"),
  ancestry = c(
    sub(".*_(eur|eas|afr|amr|sas).*", "\\1", list.files(here("data_raw", "sumstats", "case_control"), pattern = "\\.gz$")),
    sub(".*_(eur|eas|afr|amr|sas).*", "\\1", list.files(here("data_raw", "sumstats", "proxy"), pattern = "\\.gz$")),
    sub(".*_(eur|eas|afr|amr|sas).*", "\\1", list.files(here("data_raw", "sumstats", "combined"), pattern = "\\.gz$"))
  )) %>%
  mutate(
    cohort = stringr::str_extract(outfile, ".*(?=_case_control|_proxy|_combined)"),
    proxy_casecontrol = ifelse(grepl("proxy", outfile), "proxy", "case_control"),
    proxy_casecontrol = ifelse(grepl("combined", outfile), "combined", proxy_casecontrol)
  ) %>%
  rowwise() %>%
  mutate(
    ## run if any of these files are missing
    run = any(c(
      ## check if cleaned sumstats already exists
      !file.exists(here(paste0("data/sumstats/", proxy_casecontrol, "/clean_", outfile, ".txt.gz"))),
      ## check if qc plots of the raw sumstats exist
      !file.exists(here(paste0("data_raw/sumstats/qc_plots/", outfile, "_qq_plot.png"))),
      !file.exists(here(paste0("data_raw/sumstats/qc_plots/", outfile, "_manhattan_plot.png"))),
      !file.exists(here(paste0("data_raw/sumstats/qc_plots/", outfile, "_beta_histogram_plot.png"))),
      !file.exists(here(paste0("data_raw/sumstats/qc_plots/", outfile, "_beta_box_per_chr_plot.png"))),
      !file.exists(here(paste0("data_raw/sumstats/qc_plots/", outfile, "_se_box_per_chr_plot.png"))),
      !file.exists(here(paste0("data_raw/sumstats/qc_plots/", outfile, "_p_histogram_plot.png"))),
      !file.exists(here(paste0("data_raw/sumstats/qc_plots/", outfile, "_maf_histogram_plot.png"))),
      !file.exists(here(paste0("data_raw/sumstats/qc_plots/", outfile, "_neff_histogram_plot.png"))),
      !file.exists(here(paste0("data_raw/sumstats/qc_plots/", outfile, "_pz_plot.png"))),
      ## check if qc plots of the clean sumstats exist
      !file.exists(here(paste0("data/sumstats/qc_plots/clean_", outfile, "_qq_plot.png"))),
      !file.exists(here(paste0("data/sumstats/qc_plots/clean_", outfile, "_manhattan_plot.png"))),
      !file.exists(here(paste0("data/sumstats/qc_plots/clean_", outfile, "_beta_histogram_plot.png"))),
      !file.exists(here(paste0("data/sumstats/qc_plots/clean_", outfile, "_beta_box_per_chr_plot.png"))),
      !file.exists(here(paste0("data/sumstats/qc_plots/clean_", outfile, "_se_box_per_chr_plot.png"))),
      !file.exists(here(paste0("data/sumstats/qc_plots/clean_", outfile, "_p_histogram_plot.png"))),
      !file.exists(here(paste0("data/sumstats/qc_plots/clean_", outfile, "_maf_histogram_plot.png"))),
      !file.exists(here(paste0("data/sumstats/qc_plots/clean_", outfile, "_neff_histogram_plot.png"))),
      !file.exists(here(paste0("data/sumstats/qc_plots/clean_", outfile, "_pz_plot.png")))      
      ))
  ) %>%
  ungroup() %>%
  filter(run == TRUE)

setwd(here("data_raw"))

for (cohort in unique(params$cohort)) {
  
  save_job(
    job = paste0(
      "#!/bin/bash

      #SBATCH --nodes=1
      #SBATCH -t 10:00:00
      #SBATCH --job-name=03_", cohort, " 
      #SBATCH --output=logs/03_process_raw_sumstats_", cohort, "_%A_out.txt
      #SBATCH --error=logs/03_process_raw_sumstats_", cohort, "_%A_error.txt
      #SBATCH --partition=rome
      #SBATCH --mem=48G
      
      ## load modules
      module load 2022
      module load R/4.2.1-foss-2022a
      
      Rscript 03_process_raw_sumstats.R ", cohort, " 
      "
   ),
   directory = here("data_raw/"))
  
  ## submit job
  system(paste0("sbatch ", here("data_raw/job.sh")))
  
  ## remove job script
  system(paste0("rm ", here("data_raw/job.sh")))
  
}