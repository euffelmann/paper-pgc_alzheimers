#### 0. setup ####
library(here)
library(data.table)
library(dplyr)
source(here("R/save_job.R"))

#### 1. set parameters ####
params <- data.frame(
  sumstats_file = c(
    list.files(here("data/sumstats/meta/male/combined"), pattern = "nsumstats[0-9]+\\.txt\\.gz$", full.names = T),
    list.files(here("data/sumstats/meta/male/case_control"), pattern = "nsumstats[0-9]+\\.txt\\.gz$", full.names = T),
    list.files(here("data/sumstats/meta/male/proxy"), pattern = "nsumstats[0-9]+\\.txt\\.gz$", full.names = T)
  ),
  outfile = c(
    sub("\\.txt\\.gz$", "", list.files(here("data/sumstats/meta/male/combined"), pattern = "nsumstats[0-9]+\\.txt\\.gz$")),
    sub("\\.txt\\.gz$", "", list.files(here("data/sumstats/meta/male/case_control"), pattern = "nsumstats[0-9]+\\.txt\\.gz$")),
    sub("\\.txt\\.gz$", "", list.files(here("data/sumstats/meta/male/proxy"), pattern = "nsumstats[0-9]+\\.txt\\.gz$"))
  ),
  qc_dir = ".",
  ancestry = c(
    sub(".*_(eur|eas|afr|amr|sas|all).*", "\\1", list.files(here("data/sumstats/meta/male/combined"), pattern = "nsumstats[0-9]+\\.txt\\.gz$")),
    sub(".*_(eur|eas|afr|amr|sas|all).*", "\\1", list.files(here("data/sumstats/meta/male/case_control"), pattern = "nsumstats[0-9]+\\.txt\\.gz$")),
    sub(".*_(eur|eas|afr|amr|sas|all).*", "\\1", list.files(here("data/sumstats/meta/male/proxy"), pattern = "nsumstats[0-9]+\\.txt\\.gz$"))
  ))

#### 2. job file ####
setwd(here("data/sumstats/meta/male"))
for (i in 1:nrow(params)) {
  
  save_job(
    job = paste0(
      "#!/bin/bash

      #SBATCH --nodes=1
      #SBATCH -t 02:00:00
      #SBATCH --job-name=qc_", params$outfile[i], "  
      #SBATCH --output=logs/qc_", params$outfile[i], "_%A_out.txt
      #SBATCH --error=logs/qc_", params$outfile[i], "_%A_error.txt
      #SBATCH --partition=rome
      #SBATCH --mem=32G
      
      ## load modules
      module load 2022
      module load R/4.2.1-foss-2022a
      
      ## copy files to $TMPDIR
      cp -r ", here("R/qc_plots.R") , " $TMPDIR/
      cp -r ", here("R/qc_plots_job.R") , " $TMPDIR/
      mkdir $TMPDIR/data
      cp -r ", here("data/reference_data") , " $TMPDIR/data/
      cp -r ", params$sumstats_file[i], " $TMPDIR/ 
      cp -r ", here("data/sumstats/meta/male/qc_plots/*png"), " $TMPDIR/ 
      
      ## save current dir
      qc_dir=\"", here("data/sumstats/meta/male/qc_plots"), "\"
      
      cd $TMPDIR
      
      ## run script
      Rscript qc_plots_job.R ", paste0(params$outfile[i], ".txt.gz"), " ", params$outfile[i], " ", params$qc_dir[i], " ", params$ancestry[i], " 
      
      ## copy all plots back to home directory
      cp *png $qc_dir
      
      ## remove everything again
      rm -r $TMPDIR/*
      
      "
    ),
    directory = here("data/sumstats/meta/male/"))
  
  ## submit job
  system(paste0("sbatch ", here("data/sumstats/meta/male/job.sh")))
  
  ## remove job script
  system(paste0("rm ", here("data/sumstats/meta/male/job.sh")))
  
}