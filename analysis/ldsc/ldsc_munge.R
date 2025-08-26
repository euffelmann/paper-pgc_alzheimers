#### 0. set-up ####
library(data.table)
library(dplyr)
library(here)
library(openxlsx)
setwd(here("analysis/ldsc/ldsc_munge"))

#### 1. munge individual cohorts #####
params <- readWorkbook(here("analysis/cohort_summary/cohort_summary.xlsx"), sheet = "main") %>%
  mutate(
    subdir = ifelse(proxy_casecontrol == "proxy", "proxy_meta", proxy_casecontrol),
    ancestry_ldsc = ifelse(ancestry == "sas", "csa", ancestry)
  ) %>%
  mutate(
    sumstats_file = paste0(
      here("data/sumstats/"),
      subdir,
      "/clean_",
      cohort,
      "_",
      proxy_casecontrol,
      "_",
      analysis,
      "_",
      ancestry,
      ".txt.gz"
    )
  )

for (i in 1:nrow(params)) {
  for (apoe in c("_yeaapoe", "_nayapoe")) { 
    
    cohort <- params$cohort[i]
    subdir <- params$subdir[i]
    analysis <- params$analysis[i]
    ancestry <- params$ancestry[i]
    ancestry_ldsc <- toupper(params$ancestry_ldsc[i])
    sumstats_file <- params$sumstats_file[i]
    outfile <- gsub(here(paste0("data/sumstats/", subdir, "/")), "", sumstats_file)
    outfile <- gsub(".txt.gz", "", outfile)
    outfile <- gsub("clean_", "", outfile)
    outfile <- paste0(outfile, apoe)
    
    if (!file.exists(paste0(here("analysis/ldsc/ldsc_munge/coho/"), outfile, "_munge.sumstats.gz"))) {
      
      ## make temporary files with n_case and n_control columns removed, 
      ## because they interfere with ldsc (i.e. their presence makes ldsc
      ## ignore the --N-col flag). make a version with the APOE locus removed
      print(paste0("making temporary sumstats files for ... ", outfile))
      sst <- fread(sumstats_file)
      if (apoe == "_yeaapoe") {
        sst %>%
          select(-c(n_case, n_control)) %>%
          mutate(p_value = as.numeric(p_value)) %>%
          mutate(p_value = ifelse(p_value < 1e-310, 1e-310, p_value)) %>%
          fwrite(
            file = paste0(here("analysis/ldsc/ldsc_munge/temp_files/"), outfile, ".txt.gz"),
            sep = "\t",
            na = "NA",
            quote = FALSE
          )
      } else {
        sst %>%
          select(-c(n_case, n_control)) %>%
          mutate(p_value = as.numeric(p_value)) %>%
          mutate(p_value = ifelse(p_value < 1e-310, 1e-310, p_value)) %>%
          # remove 1Mb window around APOE-e4 SNP on GRCh37 / hg19 coordinates
          filter(!(chromosome == 19 & base_pair_location > (45411941 - 5e5) & base_pair_location < (45411941 + 5e5))) %>%
          fwrite(
            file = paste0(here("analysis/ldsc/ldsc_munge/temp_files/"), outfile, ".txt.gz"),
            sep = "\t",
            na = "NA",
            quote = FALSE
          )
      }
      
      system(
        paste0(
          "sbatch ", here("src/ldscMunge.sh"), " ",
          "-f ", sumstats_file, " ",
          "-e ", outfile, " ",
          "-o ", here("analysis/ldsc/ldsc_munge/coho"), " ",
          "-m '--N-col neff --a1 effect_allele --a2 other_allele --snp rsid' "
        )
      ) 
    }
  }
}


#### 2. munge meta sumstats #####
for (analysis in c("main", "fema", "male", "noag", "apoe", "agex")) {
  for (subdir in c("case_control", "proxy", "combined")) {
    for (ancestry in c("eur", "afr", "eas", "sas", "amr")) {
      for (apoe in c("_yeaapoe")) {
        
        sumstats_file <- here(
          paste0( "data/sumstats/meta/", analysis, "/", subdir, "/", analysis,  "_", subdir, "_", ancestry, ".txt.gz")
        )
        
        # check if sumstats file exists
        if (file.exists(sumstats_file)) {
          
          outfile <- paste0(analysis, "_", subdir, "_", ancestry, apoe)
          
          if (!file.exists(paste0(here("analysis/ldsc/ldsc_munge/meta", analysis, subdir, ""), outfile, "_munge.sumstats.gz"))) {
            
            ## make temporary files with n_case and n_control columns removed, 
            ## because they interfere with ldsc (i.e. their presence makes ldsc
            ## ignore the --N-col flag). make a version with the APOE locus removed
            print(paste0("making temporary sumstats files for ... ", outfile))
            sst <- fread(sumstats_file)
            if (apoe == "_yeaapoe") {
              sst %>%
                select(-c(n_case, n_control)) %>%
                mutate(p_value = as.numeric(p_value)) %>%
                mutate(p_value = ifelse(p_value < 1e-310, 1e-310, p_value)) %>%
                fwrite(
                  file = paste0(here("analysis/ldsc/ldsc_munge/temp_files/"), outfile, ".txt.gz"),
                  sep = "\t",
                  na = "NA",
                  quote = FALSE
                )
            } else {
              sst %>%
                select(-c(n_case, n_control)) %>%
                mutate(p_value = as.numeric(p_value)) %>%
                mutate(p_value = ifelse(p_value < 1e-310, 1e-310, p_value)) %>%
                # remove 1Mb window around APOE-e4 SNP on GRCh37 / hg19 coordinates
                filter(!(chromosome == 19 & base_pair_location > (45411941 - 5e5) & base_pair_location < (45411941 + 5e5))) %>%
                fwrite(
                  file = paste0(here("analysis/ldsc/ldsc_munge/temp_files/"), outfile, ".txt.gz"),
                  sep = "\t",
                  na = "NA",
                  quote = FALSE
                )
            }
            
            system(
              paste0(
                "sbatch ", here("src/ldscMunge.sh"), " ",
                "-f ", paste0(here("analysis/ldsc/ldsc_munge/temp_files/"), outfile, ".txt.gz"), " ",
                "-e ", outfile, " ",
                "-o ", here("analysis/ldsc/ldsc_munge/meta", analysis, subdir), " ",
                "-m '--N-col neff --a1 effect_allele --a2 other_allele --snp rsid' "
              )
            ) 
          }
        }
      }
    }
  }
}



#### 3. munge other sumstats #####
phenos <- read.xlsx(here("analysis/lava/other_sumstats/other_sumstats_filepaths.xlsx"), sheet = "2025_05_13") %>%
  filter(!is.na(raw_file), clean_file == "pass") %>%
  mutate(
    pheno_unique = paste(pheno_short, tolower(Author), Year, sep = "_"),
    pheno_unique = gsub(" ", "_", pheno_unique)
  ) %>%
  pull(pheno_unique)

phenos <- c(phenos, "longevity_deelen_2019", "amyloid_beta_jansen_2022")

for (pheno in phenos) {
    
  system(
    paste0(
      "sbatch ", here("src/ldscMunge.sh"), " ",
      "-f ", paste0(here("analysis/lava/other_sumstats/clean_sumstats/"), pheno, ".txt.gz"), " ",
      "-e ", pheno, " ",
      "-o ", here("analysis/ldsc/ldsc_munge/other"), " ",
      "-m '--N-col N --a1 A1 --a2 A2 --snp SNP' "
    )
  ) 
  
}