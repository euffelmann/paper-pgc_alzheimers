#### 0. set-up ####
library(data.table)
library(dplyr)
library(here)
plink1.9 <- here("src/plink1.9/plink1.9")


#### 1. plink --score for Swedish cohorts ####
for (cohort in c(
  "gothe",
  "demge",
  "xstsa",
  "twing"
)) {
  
  for (apoe in c("yea", "nay")) {
    
    for (anc in c("eur", "eas", "afr")) {
      
      for (q in c("_q0.9", "_q0.95")) { 
        
        if (anc != "afr") {q <- ""}
        
        ## if the sumstats ancestry is European, use loo sumstats
        if (anc == "eur") {
          post_betas <- here(
            "analysis/pgs/sbayesrc/post_betas", 
            paste0(cohort, "_combined_", anc, "_neff0.6_nsumstats1_lia_rc_multichain_idhg38_", apoe, "apoe.snpRes")
          )
          
          ## if the sumstats ancestry is not European, use full sumstats
        } else if (anc != "eur") {
          post_betas <- here(
            "analysis/pgs/sbayesrc/post_betas", 
            paste0("main_combined_", anc, "_neff0.6_nsumstats1_lia_rc_multichain", q, "_idhg38_", apoe, "apoe.snpRes")
          )
        }
        
        outfile <- here("analysis/pgs/sbayesrc/plink_scores", paste0(cohort, "_", anc, q, "_", apoe, "apoe"))
        testing_sample <- here("analysis/pgs/testing_samples", cohort)
        bim_file <- here("analysis/pgs/sbayesrc/input", paste0(cohort, ".bim"))
        
        
        if (!file.exists(paste0(outfile, ".profile"))) {
          
          print(paste("running plink --score for", cohort, apoe, anc, q))
          
          system(
            paste0(
              plink1.9, " ",
              "--bfile ", testing_sample, " ",
              "--bim ", bim_file, " ",
              "--score ", post_betas, " 21 5 8 header sum center ",
              "--memory 25000 ", 
              "--out ", outfile
            )
          )  
          
        }
      }
    }
  }
}

#### 2. plink --score for UKB ####

for (cohort in c(
  "xukbb_eur_casec",
  "xukbb_afr_casec",
  "xukbb_sas_casec"
)) {
  
  for (apoe in c("yea", "nay")) {
    
    for (anc in c("eur", "eas", "afr")) {
      
      for (q in c("_q0.9", "_q0.95")) { 
        
        if (anc != "afr") {q <- ""}
        
        ## if the UKB ancestry matches the sumstats ancestry, use loo sumstats
        if (grepl(anc, cohort)) {
          post_betas <- here(
            "analysis/pgs/sbayesrc/post_betas", 
            paste0("xukbb_combined_", anc, "_neff0.6_nsumstats1_lia_rc_multichain", q, "_idhg37_", apoe, "apoe.snpRes")
          )
        }
        
        ## if the UKB ancestry does not match the sumstats ancestry, use full sumstats
        if (!grepl(anc, cohort)) { 
          post_betas <- here(
            "analysis/pgs/sbayesrc/post_betas", 
            paste0("main_combined_", anc, "_neff0.6_nsumstats1_lia_rc_multichain", q, "_idhg37_", apoe, "apoe.snpRes")
          )
        }
        
        outfile <- here("analysis/pgs/sbayesrc/plink_scores", paste0(cohort, "_", anc, q, "_", apoe, "apoe"))
        testing_sample <- here("analysis/pgs/testing_samples", cohort)
        bim_file <- here("analysis/pgs/sbayesrc/input", paste0(cohort, ".bim"))
        
        
        if (!file.exists(paste0(outfile, ".profile"))) {
          
          print(paste("running plink --score for", cohort, apoe, anc, q))
          
          system(
            paste0(
              plink1.9, " ",
              "--bfile ", testing_sample, " ",
              "--bim ", bim_file, " ",
              "--score ", post_betas, " 20 5 8 header sum center ",
              "--memory 25000 ", 
              "--out ", outfile
            )
          )  
          
        }
      }
    }
  }
}