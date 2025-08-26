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
  filter(
    # filter for main analysis or those that didn't have age available
    grepl("main", sumstats_file) | 
      (
        grepl("noag", sumstats_file) & 
          (
            grepl("mvplo", sumstats_file) |
              grepl("biovu", sumstats_file) |
              grepl("grace", sumstats_file) | 
              grepl("xstsa", sumstats_file) | 
              grepl("twing", sumstats_file) | 
              grepl("gothe", sumstats_file)
          )
      )
  )
params_proxy <- data.frame(
  sumstats_file = list.files(
    here("data", "sumstats", "proxy_meta"),
    pattern = "\\.gz$",
    full.names = T
  ),
  ancestry = c(
    sub(".*_(eur|eas|afr|amr|sas).*", "\\1", list.files(here("data", "sumstats", "proxy_meta"), pattern = "\\.gz$"))) 
) %>%
  filter(
    grepl("main", sumstats_file),
    !grepl("inref", sumstats_file)
  )

params_mvp_mothers <- data.frame(
  sumstats_file = list.files(
    here("data", "sumstats", "proxy"),
    pattern = "\\.gz$",
    full.names = T
  ),
  ancestry = c(
    sub(".*_(eur|eas|afr|amr|sas).*", "\\1", list.files(here("data", "sumstats", "proxy"), pattern = "\\.gz$"))) 
) %>%
  filter(
    grepl("mvplo_proxy_mother", sumstats_file),
    !grepl("inref", sumstats_file)
  )

params <- rbind(params_cc, params_proxy, params_mvp_mothers)
sumstats <- paste(params$sumstats_file, collapse = ",")
outfile <- "main"
out_dir <- here("data/sumstats/meta/main")
scheme <- "STDERR"

setwd(here("data/sumstats/meta/main"))

#### 2. run meta-anlaysis ####
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
      cp -r ", here("src", "METAL-2020-05-05") , " $TMPDIR/
      cp -r ", here("data/reference_data/g1000_afr/g1000_afr.frq") , " $TMPDIR/
      
      cd $TMPDIR
      
      Rscript meta_analysis_job.R ", sumstats, " ", out_dir, " ", outfile, " ", scheme, " 
      "
  ),
  directory = here("data/sumstats/meta/main/"))

## submit job
system(paste0("sbatch ", here("data/sumstats/meta/main/job.sh")))

## remove job script
system(paste0("rm ", here("data/sumstats/meta/main/job.sh")))


#### 3. filter sumstats ####
setwd(here("data/sumstats/meta/main/"))
for (casec in c("combined", "case_control", "proxy")) {
  for (anc in c("all", "eur", "afr", "eas", "sas", "amr")) {
    
    name <- paste0("main_", casec, "_", anc)
    input <- paste0(casec, "/", name)
    neff_th <- 0.6
    n_sumstats_th <- 1
    
    if (!file.exists(here("data/sumstats/meta/main", paste0(input, "_neff", neff_th, "_nsumstats", n_sumstats_th, ".txt.gz")))) {
      
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
      
      Rscript filter_sumstats.R ", here("data/sumstats/meta/main", input), " ", neff_th, " ", n_sumstats_th, " 
      
      ## remove everything  again
      rm -r $TMPDIR/*
      
      "
        ),
        directory = here("data/sumstats/meta/main/"))
      
      ## submit job
      system(paste0("sbatch ", here("data/sumstats/meta/main/job.sh")))
      
      ## remove job script
      system(paste0("rm ", here("data/sumstats/meta/main/job.sh")))
      
    } 
  }
}