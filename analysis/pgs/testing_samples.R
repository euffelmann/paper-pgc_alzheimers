#### 0. set-up ####
library(dplyr)
library(data.table)
library(here)
plink1.9 <- here("src/plink1.9/plink1.9")

#### 1. remove duplicate SNP IDs and move PLINK files ####
for (cohort in c("demgene_OmniExpress", "Gothenburg", "STSA", "TwinGene")) {
  
  if (cohort == "demgene_OmniExpress") {abbr <- "demge"}
  if (cohort == "Gothenburg") {abbr <- "gothe"}
  if (cohort == "STSA") {abbr <- "xstsa"}
  if (cohort == "TwinGene") {abbr <- "twing"}
  
  bim <- fread(here("random/tempshare/2025_01_16_doug_testing_samples", paste0(cohort, "mafINFO.bim"))) %>%
    mutate(variant_id = paste(V1, V4, pmin(V5, V6), pmax(V5, V6), sep = ":"))
  
  uniq_snps <- bim %>%
    distinct(variant_id, .keep_all = TRUE) %>%
    pull(V2)
    
  dups <- bim %>%
    filter(!V2 %in% uniq_snps) %>%
    pull(V2)
  
  writeLines(dups, here("random/tempshare/2025_01_16_doug_testing_samples/qc", paste0(abbr, "_dups.txt")))
  
  system(paste0(
    plink1.9, " ",
    "--bfile ", here("random/tempshare/2025_01_16_doug_testing_samples", paste0(cohort, "mafINFO")), " ",
    "--exclude ", here("random/tempshare/2025_01_16_doug_testing_samples/qc", paste0(abbr, "_dups.txt")), " ",
    "--make-bed ",
    "--memory 25000 ",
    "--out ", here("analysis/pgs/testing_samples", abbr)
  ))

}


#### 2. copy phenotype files ####
system(paste0(
  "rsync -av ",
  here("random/tempshare/2025_01_16_doug_testing_samples", paste0("GothenburgPheno.txt")),
  " ",
  here("analysis/pgs/testing_samples", paste0("gothe.cov"))
))

system(paste0(
  "rsync -av ",
  here("random/tempshare/2025_01_16_doug_testing_samples", paste0("demgene_OmniExpress_PhenoCovars.txt")),
  " ",
  here("analysis/pgs/testing_samples", paste0("demge.cov"))
))
fread(here("analysis/pgs/testing_samples", paste0("demge.cov"))) %>%
  rename(Phenotype = AD) %>%
  fwrite(file = here("analysis/pgs/testing_samples", paste0("demge.cov")), quote = F, sep = "\t", na = NA)
  
system(paste0(
  "rsync -av ",
  here("random/tempshare/2025_01_16_doug_testing_samples", paste0("TwinGenePheno.txt")),
  " ",
  here("analysis/pgs/testing_samples", paste0("twing.cov"))
))

system(paste0(
  "rsync -av ",
  here("random/tempshare/2025_01_16_doug_testing_samples", paste0("STSAPheno.txt")),
  " ",
  here("analysis/pgs/testing_samples", paste0("xstsa.cov"))
))



#### 3. extract UKB individuals ####
set.seed(20250213)
for (ancestry in c("eur", "afr", "sas")) {
  
  geno_files <- paste0(
    "/projects/0/ctgukbio/datasets/ukbio/applicationID1640/qc/final/genotypes/release2b/",
    toupper(ancestry),
    "/unrel_",
    toupper(ancestry),
    "/unrel_",
    toupper(ancestry)
  )
  
  ## extract IDs to keep
  if (ancestry == "eur") {
    
    unrel_eur <- fread(paste0(geno_files, ".fam")) %>%
      select(V2)
    
    cov_eur <- fread(here(
      "analysis/pgs/testing_samples",
      paste0("xukbb_", ancestry, "_casec.cov")
    )) %>% 
      inner_join(unrel_eur, by = c("IID" = "V2"))
    
    n_cases <- sum(cov_eur$AD == 1)
    
    cases <-
      cov_eur %>%
      filter(AD == 1) %>%
      slice_sample(n = n_cases)
    
    controls <-
      cov_eur %>%
      filter(AD == 0) %>%
      slice_sample(n = n_cases)
      
    df <- rbind(cases, controls)
    if (sum(df$AD == 1) != n_cases | sum(df$AD == 0) != n_cases) {
      
      stop("sampling cases and controls did not work")
      
    }
    
    fwrite(
      x = df[, c("FID", "IID")],
      file = here(
        "analysis/pgs/testing_samples",
        paste0("temp_xukbb_", ancestry, "_id.txt")
      ),
      quote = F,
      sep = "\t",
      col.names = F
    )
    
  } else {
    
    system(paste0(
      "awk 'NR > 1 { print $1, $2 }' ", 
      here("analysis/pgs/testing_samples", paste0("xukbb_", ancestry, "_casec.cov")), 
      " > ",
      here("analysis/pgs/testing_samples", paste0("temp_xukbb_", ancestry, "_id.txt"))))
    
  }
  
  ## make new plink files with individuals of interest
  system(paste0(
    plink1.9, " ",
    "--bfile ", geno_files, " ",
    "--keep ", here("analysis/pgs/testing_samples", paste0("temp_xukbb_", ancestry, "_id.txt")), " ",
    "--make-bed ",
    "--out ", here("analysis/pgs/testing_samples", paste0("xukbb_", ancestry, "_casec"))
  ))
  
}



#### 4. count APOE2 allele ####
## demge
cohort <- "demge"

apoe2 <- "19:44908822:C:T T"
writeLines(apoe2, here("analysis/pgs/testing_samples/recode_apoe2.txt"))

apoe2 <- "19:44908822:C:T"
writeLines(apoe2, here("analysis/pgs/testing_samples/extract_apoe2.txt"))


system(paste0(
  plink1.9, " ",
  "--bfile ", here("analysis/pgs/testing_samples", cohort), " ",
  "--extract ", here("analysis/pgs/testing_samples/extract_apoe2.txt"), " ",
  "--recode A ",
  "--memory 25000 ",
  "--recode-allele ", here("analysis/pgs/testing_samples/recode_apoe2.txt"), " ",
  "--out ", here("analysis/pgs/testing_samples", paste0(cohort, "_apoe2"))
))

## the other three swedish cohorts
for (cohort in c("xstsa", "twing", "gothe")) {
  
  apoe2 <- "chr19:44908822:C:T T"
  writeLines(apoe2, here("analysis/pgs/testing_samples/recode_apoe2.txt"))
  
  apoe2 <- "chr19:44908822:C:T"
  writeLines(apoe2, here("analysis/pgs/testing_samples/extract_apoe2.txt"))
  
  system(paste0(
    plink1.9, " ",
    "--bfile ", here("analysis/pgs/testing_samples", cohort), " ",
    "--extract ", here("analysis/pgs/testing_samples/extract_apoe2.txt"), " ",
    "--recode A ",
    "--memory 25000 ",
    "--recode-allele ", here("analysis/pgs/testing_samples/recode_apoe2.txt"), " ",
    "--out ", here("analysis/pgs/testing_samples", paste0(cohort, "_apoe2"))
  ))
  
}

## for xukbb
for (cohort in c("xukbb_eur_casec", "xukbb_afr_casec", "xukbb_sas_casec")) {
  
  apoe2 <- "19:45412079:C_T T"
  writeLines(apoe2, here("analysis/pgs/testing_samples/recode_apoe2.txt"))
  
  apoe2 <- "19:45412079:C_T"
  writeLines(apoe2, here("analysis/pgs/testing_samples/extract_apoe2.txt"))
  
  system(paste0(
    plink1.9, " ",
    "--bfile ", here("analysis/pgs/testing_samples", cohort), " ",
    "--extract ", here("analysis/pgs/testing_samples/extract_apoe2.txt"), " ",
    "--recode A ",
    "--memory 25000 ",
    "--recode-allele ", here("analysis/pgs/testing_samples/recode_apoe2.txt"), " ",
    "--out ", here("analysis/pgs/testing_samples", paste0(cohort, "_apoe2"))
  ))
  
}

