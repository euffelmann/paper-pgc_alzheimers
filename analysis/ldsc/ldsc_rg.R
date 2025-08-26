#### 0. setup ####
library(data.table)
library(dplyr)
library(here)
setwd(here("analysis/ldsc/ldsc_rg"))

#### 1. meta-analysis sumstats ####

##### i. ldsc rg: proxy vs. casec (main) #####

for (anc in c("eur", "afr", "eas", "sas", "amr")) {
  
  ## ldsc ancestries
  if (anc == "sas") {
    anc_ldsc <- "CSA"
  } else {
    anc_ldsc <- toupper(anc)
  }
  
  sst_proxy <- here(
    "analysis/ldsc/ldsc_munge/meta/main/proxy", 
    paste0("main_proxy_", anc, "_yeaapoe_munge.sumstats.gz")
  )
  
  sst_casec <- here(
    "analysis/ldsc/ldsc_munge/meta/main/case_control", 
    paste0("main_case_control_", anc, "_yeaapoe_munge.sumstats.gz")
  )
  
  if (file.exists(sst_proxy) & file.exists(sst_casec)) {
    
    system(paste0(
      "sbatch ", here("src/ldscRg.sh"), " ",
      "-f ", sst_proxy, " ",
      "-r ", sst_casec, " ",
      "-e main_proxy_main_casec_", anc, "_yeaapoe ",
      "-o ", here("analysis/ldsc/ldsc_rg/meta"), " ",
      "-l ", here("data/reference_data/UKBB.ALL.ldscore/UKBB."), anc_ldsc, ".rsid "
    ))
 
  }
}


##### ii. ldsc rg: fema vs. male #####
for (anc in c("eur", "afr", "eas", "sas", "amr")) {
  for (proxy_casec in c("proxy", "case_control", "combined")) {
    
    ## ldsc ancestries
    if (anc == "sas") {
      anc_ldsc <- "CSA"
    } else {
      anc_ldsc <- toupper(anc)
    }
    
    sst_fema <- here(
      "analysis/ldsc/ldsc_munge/meta/fema", proxy_casec, 
      paste0("fema_", proxy_casec, "_", anc, "_yeaapoe_munge.sumstats.gz")
    )
    
    sst_male <- here(
      "analysis/ldsc/ldsc_munge/meta/male", proxy_casec, 
      paste0("male_", proxy_casec, "_", anc, "_yeaapoe_munge.sumstats.gz")
    )
    
    if (file.exists(sst_fema) & file.exists(sst_male)) {
      
      system(paste0(
        "sbatch ", here("src/ldscRg.sh"), " ",
        "-f ", sst_fema, " ",
        "-r ", sst_male, " ",
        "-e fema_", proxy_casec, "_male_", proxy_casec, "_", anc, "_yeaapoe ",
        "-o ", here("analysis/ldsc/ldsc_rg/meta"), " ",
        "-l ", here("data/reference_data/UKBB.ALL.ldscore/UKBB."), anc_ldsc, ".rsid "
      ))
      
    }
  }
}


#### 2. other sumstats ####
phenos <- read.xlsx(here("analysis/lava/other_sumstats/other_sumstats_filepaths.xlsx"), sheet = "2025_05_13") %>%
  filter(!is.na(raw_file), clean_file == "pass") %>%
  mutate(
    pheno_unique = paste(pheno_short, tolower(Author), Year, sep = "_"),
    pheno_unique = gsub(" ", "_", pheno_unique)
  ) %>%
  pull(pheno_unique)

phenos <- c(phenos, "longevity_deelen_2019", "amyloid_beta_jansen_2022")

for (pheno in phenos) {
  
  for (proxy_casec in c("case_control", "combined")) {
    
    sst_other <- here(
      "analysis/ldsc/ldsc_munge/other", paste0(pheno, "_munge.sumstats.gz")
      )
    
    sst_alz <- here(
      "analysis/ldsc/ldsc_munge/meta/main", proxy_casec, 
      paste0("main_", proxy_casec, "_eur_yeaapoe_munge.sumstats.gz")
    )
    
    if (file.exists(sst_other) & file.exists(sst_alz)) {
      
      system(paste0(
        "sbatch ", here("src/ldscRg.sh"), " ",
        "-f ", sst_other, " ",
        "-r ", sst_alz, " ",
        "-e ", pheno, "_main_", proxy_casec, "_eur_yeaapoe ",
        "-o ", here("analysis/ldsc/ldsc_rg/other"), " ",
        "-l ", here("data/reference_data/UKBB.ALL.ldscore/UKBB."), "EUR.rsid "
      ))
      
    }
    
  }
}


## make results table
results_df <- NULL

for (pheno in phenos) {
  
  for (proxy_casec in c("combined", "case_control")) {
    
    ldsc <-
      fread(here(
        "analysis/ldsc/ldsc_rg/other",
        paste0(pheno, "_main_", proxy_casec, "_eur_yeaapoe_rg.table")
      ))
    
    ldsc[1, "pheno1"] <- pheno
    ldsc[1, "pheno2"] <- paste0("alzh_", proxy_casec)
     
    results_df <- rbind(results_df, ldsc)
    
  }
}

fwrite(
  x = results_df,
  file = here("analysis/ldsc/ldsc_rg/other.txt"),
  quote = F,
  sep = "\t",
  na = NA
)