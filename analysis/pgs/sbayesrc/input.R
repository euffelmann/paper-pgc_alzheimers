#### 0. set-up ####
library(data.table)
library(dplyr)
library(here)
library(stringr)


#### 1. adjust .bim files ####
for (cohort in c(
  "gothe",
  "demge",
  "xstsa",
  "twing",
  "xukbb_eur_casec",
  "xukbb_afr_casec",
  "xukbb_sas_casec"
)) {
  
  if (!file.exists(here("analysis/pgs/sbayesrc/input", paste0(cohort, ".bim")))) {
    
    print(paste0("making .bim file for ", cohort))
    
    fread(here("analysis/pgs/testing_samples", paste0(cohort, ".bim"))) %>%
      mutate(V2 = paste(V1, V4, pmin(V5, V6), pmax(V5, V6), sep = ":")) %>%
      fwrite(
        file = here("analysis/pgs/sbayesrc/input", paste0(cohort, ".bim")),
        quote = F,
        sep = "\t",
        na = NA,
        col.names = F
      ) 
    
  }
}
