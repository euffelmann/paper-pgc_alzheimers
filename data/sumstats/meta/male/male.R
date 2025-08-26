#### 0. setup ####
library(here)
library(dplyr)
library(data.table)
source(here("R/save_job.R"))

#### 1. define input ####

## sumstats to include in meta-analysis
params_cc <- data.frame(
  sumstats_file = list.files(
    here("data", "sumstats", "case_control"),
    pattern = "\\.gz$",
    full.names = T
  ),
  ancestry = c(
    sub(".*_(eur|eas|afr|amr|sas).*", "\\1", list.files(here("data", "sumstats", "case_control"), pattern = "\\.gz$"))) 
) %>%
  filter(grepl("male", sumstats_file))

params_proxy <- data.frame(
  sumstats_file = list.files(
    here("data", "sumstats", "proxy"),
    pattern = "\\.gz$",
    full.names = T
  ),
  ancestry = c(
    sub(".*_(eur|eas|afr|amr|sas).*", "\\1", list.files(here("data", "sumstats", "proxy"), pattern = "\\.gz$"))) 
) %>%
  filter(grepl("father", sumstats_file) & grepl("main", sumstats_file))

params_mvp_fathers <- data.frame(
  sumstats_file = list.files(
    here("data", "sumstats", "proxy"),
    pattern = "\\.gz$",
    full.names = T
  ),
  ancestry = c(
    sub(".*_(eur|eas|afr|amr|sas).*", "\\1", list.files(here("data", "sumstats", "proxy"), pattern = "\\.gz$"))) 
) %>%
  filter(
    grepl("mvplo_proxy_father", sumstats_file),
    !grepl("inref", sumstats_file)
  )

params <- rbind(params_cc, params_proxy, params_mvp_fathers)
sumstats <- paste(params$sumstats_file, collapse = ",")
outfile <- "male"
out_dir <- here("data/sumstats/meta/male")
scheme <- "STDERR"

setwd(here("data/sumstats/meta/male"))

#### 2. submit job to snellius ####
save_job(
  job = paste0(
    "#!/bin/bash

      #SBATCH --nodes=1
      #SBATCH -t 03:00:00
      #SBATCH --job-name=meta_", outfile, "  
      #SBATCH --output=logs/meta_", outfile, "_%A_out.txt
      #SBATCH --error=logs/meta_", outfile, "_%A_error.txt
      #SBATCH --partition=rome
      #SBATCH --mem=32G
      
      ## load modules
      module load 2023
      module load R/4.3.2-gfbf-2023a
      
      ## copy files to $TMPDIR
      echo 'current directory: $TMPDIR'
      cp -r ", here("R/*") , " $TMPDIR/
      cp -r ", out_dir , " $TMPDIR/
      cp -r ", here("src/METAL-2020-05-05") , " $TMPDIR/
      cp -r ", here("data/reference_data/g1000_afr/g1000_afr.frq") , " $TMPDIR/
      
      cd $TMPDIR
      
      Rscript meta_analysis_job.R ", sumstats, " ", out_dir, " ", outfile, " ", scheme, " 
      "
  ),
  directory = here("data/sumstats/meta/male/"))

## submit job
system(paste0("sbatch ", here("data/sumstats/meta/male/job.sh")))

## remove job script
system(paste0("rm ", here("data/sumstats/meta/male/job.sh")))


#### 3. filter sumstats ####
setwd(here("data/sumstats/meta/male/"))
for (casec in c("combined", "case_control", "proxy")) {
  for (anc in c("all", "eur", "afr", "eas", "sas", "amr")) {
    
    name <- paste0("male_", casec, "_", anc)
    input <- paste0(casec, "/", name)
    neff_th <- 0.6
    n_sumstats_th <- 1
    
    if (!file.exists(here("data/sumstats/meta/male", paste0(input, "_neff", neff_th, "_nsumstats", n_sumstats_th, ".txt.gz")))) {
      
      save_job(
        job = paste0(
          "#!/bin/bash
    
      #SBATCH --nodes=1
      #SBATCH -t 00:05:00
      #SBATCH --job-name=filter_", name, "  
      #SBATCH --output=logs/filter_", name, "_%A_out.txt
      #SBATCH --error=logs/filter_", name, "_%A_error.txt
      #SBATCH --partition=rome
      #SBATCH --mem=28G
      
      ## load modules
      module load 2023
      module load R/4.3.2-gfbf-2023a
      
      ## copy files to $TMPDIR
      echo 'current directory: $TMPDIR'
      cp -r ", here("R/*") , " $TMPDIR/

      cd $TMPDIR
      
      Rscript filter_sumstats.R ", here("data/sumstats/meta/male", input), " ", neff_th, " ", n_sumstats_th, " 
      "
        ),
        directory = here("data/sumstats/meta/male/"))
      
      ## submit job
      system(paste0("sbatch ", here("data/sumstats/meta/male/job.sh")))
      
      ## remove job script
      system(paste0("rm ", here("data/sumstats/meta/male/job.sh")))
      
    } 
  }
}