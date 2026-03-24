#### 0. set-up ####
library(here)
library(dplyr)
library(data.table)
library(openxlsx)
source(here("R/save_job.R"))


#### 1. risk loci tables ####
setwd(here("analysis/risk_loci"))
for (analysis in c("main", "male", "fema")) {
  
  system(paste0("mkdir ", here("analysis/risk_loci", analysis)))
  
  for (casec in c("combined", "case_control", "proxy")) {
    
    input_files <- list.files(
      path = here("data/sumstats/meta", analysis, casec), 
      pattern = "nsumstats[0-9]+\\.txt\\.gz$"
    )
    
    for (input in input_files) {
      
      input <- gsub(x = input, replacement = "", pattern = ".txt.gz")
      clump_kb <- 250 
      
      if (!file.exists(here("analysis/risk_loci", analysis, paste0(input, "_loci.txt")))) {
        
        save_job(
          job = paste0(
            "#!/bin/bash
    
      #SBATCH --nodes=1
      #SBATCH -t 00:03:00
      #SBATCH --job-name=loci_", input, "  
      #SBATCH --output=logs/loci_", input, "_%A_out.txt
      #SBATCH --error=logs/loci_", input, "_%A_error.txt
      #SBATCH --partition=rome
      #SBATCH --mem=28G
      
      ## load modules
      module load 2023
      module load R/4.3.2-gfbf-2023a
      
      ## copy files to $TMPDIR
      echo 'current directory: $TMPDIR'
      cp -r ", here("analysis/risk_loci/define_loci.R") , " $TMPDIR/
      cp -r ", here("R") , " $TMPDIR/
      cd $TMPDIR
      
      Rscript define_loci.R ", here("data/sumstats/meta", analysis, casec, input), ".txt.gz ", here("analysis/risk_loci", analysis, input), "_loci.txt ", clump_kb, " 
      "
          ),
          directory = here("analysis/risk_loci/"))
        
        ## submit job
        system(paste0("sbatch ", here("analysis/risk_loci/job.sh")))
        
        ## remove job script
        system(paste0("rm ", here("analysis/risk_loci/job.sh")))
        
      }
    }
  }
}


#### 2. risk loci tables (loo) ####
cohorts <- read.xlsx(here("analysis/cohort_summary/cohort_summary.xlsx"), sheet = "main")$cohort
cohorts <- unique(cohorts)
setwd(here("analysis/risk_loci"))

for (cohort in cohorts) {
  
  system(paste0("mkdir ", here("analysis/risk_loci/loo", cohort)))
  
  for (casec in c("combined", "case_control", "proxy")) {
    
    input_files <- list.files(
      path = here("data/sumstats/meta/loo", cohort, casec), 
      pattern = "nsumstats[0-9]+\\.txt\\.gz$"
    )
    
    for (input in input_files) {
      
      input <- gsub(x = input, replacement = "", pattern = ".txt.gz")
      clump_kb <- 250 
      
      if (!file.exists(here("analysis/risk_loci/loo", cohort, paste0(input, "_loci.txt")))) {
        
        save_job(
          job = paste0(
            "#!/bin/bash
    
      #SBATCH --nodes=1
      #SBATCH -t 00:03:00
      #SBATCH --job-name=loo_loci_", input, "  
      #SBATCH --output=logs/loo_loci_", input, "_%A_out.txt
      #SBATCH --error=logs/loo_loci_", input, "_%A_error.txt
      #SBATCH --partition=rome
      #SBATCH --mem=28G
      
      ## load modules
      module load 2023
      module load R/4.3.2-gfbf-2023a
      
      ## copy files to $TMPDIR
      echo 'current directory: $TMPDIR'
      cp -r ", here("analysis/risk_loci/define_loci.R") , " $TMPDIR/
      cp -r ", here("R") , " $TMPDIR/
      cd $TMPDIR
      
      Rscript define_loci.R ", here("data/sumstats/meta/loo", cohort, casec, input), ".txt.gz ", here("analysis/risk_loci/loo", cohort, input), "_loci.txt ", clump_kb, " 
      "
          ),
          directory = here("analysis/risk_loci/"))
        
        ## submit job
        system(paste0("sbatch ", here("analysis/risk_loci/job.sh")))
        
        ## remove job script
        system(paste0("rm ", here("analysis/risk_loci/job.sh")))
        
      }
    }
  }
}

  

