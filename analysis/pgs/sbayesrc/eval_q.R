#### 0. set-up ####
library(data.table)
library(dplyr)
library(here)
library(pROC)
source(here("R/r2l_r2o.R"))

#### 1. load data ####
pgs_0.90 <- fread(here("analysis/pgs/sbayesrc/plink_scores", "xukbb_afr_casec_afr_q0.9_yeaapoe.profile")) %>%
  dplyr::select(FID, IID, SCORESUM) %>%
  rename_with(~"PGS_q0.9", .cols = "SCORESUM")

pgs_0.95 <- fread(here("analysis/pgs/sbayesrc/plink_scores", "xukbb_afr_casec_afr_q0.95_yeaapoe.profile")) %>%
  dplyr::select(FID, IID, SCORESUM) %>%
  rename_with(~"PGS_q0.95", .cols = "SCORESUM")

cohort <- "xukbb_afr_casec"
pheno <- fread(here("analysis/pgs/testing_samples", paste0("xukbb_afr_casec", ".cov"))) %>%
  rename(
    Phenotype = ifelse(grepl("xukbb", cohort), "AD", "Phenotype"),
    PC1 = ifelse(grepl("xukbb", cohort), "pop_pc1", "PC1"),
    PC2 = ifelse(grepl("xukbb", cohort), "pop_pc2", "PC2"),
    PC3 = ifelse(grepl("xukbb", cohort), "pop_pc3", "PC3"),
    PC4 = ifelse(grepl("xukbb", cohort), "pop_pc4", "PC4"),
    PC5 = ifelse(grepl("xukbb", cohort), "pop_pc5", "PC5"),
    PC6 = ifelse(grepl("xukbb", cohort), "pop_pc6", "PC6"),
    PC7 = ifelse(grepl("xukbb", cohort), "pop_pc7", "PC7"),
    PC8 = ifelse(grepl("xukbb", cohort), "pop_pc8", "PC8"),
    PC9 = ifelse(grepl("xukbb", cohort), "pop_pc9", "PC9"),
    PC10 = ifelse(grepl("xukbb", cohort), "pop_pc10", "PC10")) %>%
  filter(Phenotype %in% c(0, 1))

m <- pheno %>%
  inner_join(pgs_0.90, by = c("FID", "IID")) %>%
  inner_join(pgs_0.95, by = c("FID", "IID"))

## q0.9 performs better 
summary(lm(Phenotype ~ PGS_q0.9, data = m))$r.sq # 0.07596734
summary(lm(Phenotype ~ PGS_q0.95, data = m))$r.sq # 0.0559696