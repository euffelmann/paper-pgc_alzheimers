#### 0. set-up ####
library(data.table)
library(dplyr)
library(here)
library(pROC)
source(here("R/r2l_r2o.R"))


#### 1. SBayesRC: evaluate EUR PGS (base model: PC[1-10]) ####
pgs_eval <- data.frame(
  cohort = as.character(),
  method = as.character(),
  model = as.character(),
  r2o = as.numeric(),
  r2l = as.numeric(),
  auc = as.numeric(),
  auc_analytic = as.numeric()
)

pgs_models <- c(
  "PC[1-10]",
  "APOE",
  "PGS",
  "APOE + PGS",
  "APOE + PGS + sex"
)

pgs_model <- NULL
pgs_model[[pgs_models[1]]] <- "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
pgs_model[[pgs_models[2]]] <- "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + APOE4 + APOE2 "
pgs_model[[pgs_models[3]]] <- "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + pgs_yeaapoe"
pgs_model[[pgs_models[4]]] <- "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + APOE4 + APOE2 + pgs_nayapoe"
pgs_model[[pgs_models[5]]] <- "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + APOE4 + APOE2 + pgs_nayapoe + sex"

for (cohort in c("xstsa", "gothe", "twing", "demge", "xukbb_eur_casec")) { 
  
  ### load pgs scores 
  ## (yeaapoe = including APOE region; nayapoe = excluding APOE region)
  profile_yeaapoe <- fread(here("analysis/pgs/sbayesrc/plink_scores", paste0(cohort, "_eur_yeaapoe.profile"))) %>%
    dplyr::select(FID, IID, SCORESUM) %>%
    rename(pgs_yeaapoe = SCORESUM)
  profile_nayapoe <- fread(here("analysis/pgs/sbayesrc/plink_scores", paste0(cohort, "_eur_nayapoe.profile"))) %>%
    dplyr::select(FID, IID, SCORESUM) %>%
    rename(pgs_nayapoe = SCORESUM)
  
  ## pheno file
  pheno <- fread(here("analysis/pgs/testing_samples", paste0(cohort, ".cov"))) %>%
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
  
  apoe2 <- fread(here("analysis/pgs/testing_samples", paste0(cohort, "_apoe2.raw")))[,c(1,2,7)]
  colnames(apoe2)[3] <- "APOE2"
  
  pheno <- pheno %>%
    inner_join(apoe2, by = c("FID", "IID")) %>%
    filter(!is.na(APOE2))
  
  if (cohort == "demge") {pheno <- pheno %>% filter(age >= 40)}
  
  ## merge files
  m <- pheno %>%
    inner_join(profile_yeaapoe, by = c("FID", "IID")) %>%
    inner_join(profile_nayapoe, by = c("FID", "IID"))
  
  for (i in 1:length(pgs_models)) {
    
    formula_r2 <- as.formula(paste("Phenotype ~", pgs_model[[pgs_models[i]]]))
    r2o <- summary(lm(formula_r2, data = m))$adj.r.sq
    r2l <- prs_r2obs_to_r2liab(
      K = 0.05,
      P = mean(m$Phenotype),
      prs_r2obs = r2o
    )
    mod1 <- glm(formula_r2, data = m, family = "binomial")
    auc <- as.numeric(roc(response = m$Phenotype, predictor = predict(mod1), quiet = TRUE)$auc)

    
    ## r2o to auc
    n1    <- sum(m$Phenotype == 1); n2 <- (sum(m$Phenotype == 0))
    alpha <- (n1 + n2)^2 / (n1 * n2)
    d     <- (sqrt(alpha) * sqrt(r2o)) / sqrt(1 - r2o)
    auc_analytic   <- pnorm(d/sqrt(2), 0, 1)
    
    
    pgs_eval[nrow(pgs_eval) + 1, "cohort"] <- cohort
    pgs_eval[nrow(pgs_eval), "method"] <- "sbayesrc"
    pgs_eval[nrow(pgs_eval), "model"] <- pgs_models[i]
    pgs_eval[nrow(pgs_eval), "r2o"] <- r2o
    pgs_eval[nrow(pgs_eval), "r2l"] <- r2l
    pgs_eval[nrow(pgs_eval), "auc"] <- auc 
    pgs_eval[nrow(pgs_eval), "auc_analytic"] <- auc_analytic 
    
  }
}

pgs_eval_lst <- NULL
pgs_eval_lst[["sbayesrc"]] <- pgs_eval


#### 2. SBayesRC-mult: evaluate EUR+EAS+AFR PGS (base model: PC[1-10]) ####
## split each cohort into tuning (1/3) and testing (2/3)
set.seed(20250124)
pgs_eval <- data.frame(
  cohort = as.character(),
  method = as.character(),
  model = as.character(),
  r2o = as.numeric(),
  r2l = as.numeric(),
  auc = as.numeric(),
  auc_analytic = as.numeric()
)

pgs_models <- c(
  "PC[1-10]",
  "PGS_eur",
  "PGS_eas",
  "PGS_afr",
  "PGS_all"
)

pgs_model <- NULL
pgs_model[[pgs_models[1]]] <- "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
pgs_model[[pgs_models[2]]] <- "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PGS_eur"
pgs_model[[pgs_models[3]]] <- "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PGS_eas"
pgs_model[[pgs_models[4]]] <- "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PGS_afr"
pgs_model[[pgs_models[5]]] <- "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PGS_all"

for (cohort in c("xstsa", "gothe", "twing", "demge", "xukbb_eur_casec", "xukbb_sas_casec", "xukbb_afr_casec")) { 
  
  ## pgs files
  profiles <- NULL
  for (anc in c("eur", "eas", "afr")) {
    
    if (anc == "afr") {q <- "_q0.9"} else {q <- ""}
    
    profiles[[anc]] <- fread(here("analysis/pgs/sbayesrc/plink_scores", paste0(cohort, "_", anc, q, "_yeaapoe.profile"))) %>%
      dplyr::select(FID, IID, SCORESUM) %>%
      rename_with(~paste0("PGS_", anc), .cols = "SCORESUM")
    
  }
  
  ## pheno file
  pheno <- fread(here("analysis/pgs/testing_samples", paste0(cohort, ".cov"))) %>%
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
  
  if (cohort == "demge") {pheno <- pheno %>% filter(age >= 40)}
  
  ## merge files
  m <- pheno %>%
    inner_join(profiles[["eur"]], by = c("FID", "IID")) %>%
    inner_join(profiles[["eas"]], by = c("FID", "IID")) %>%
    inner_join(profiles[["afr"]], by = c("FID", "IID"))
  
  ## split into tuning and testing set
  index <- sample(1:nrow(m), size = round((1/3) * nrow(m)), replace = F)
  tune_set <- m[index,]
  test_set <- m[!index,]
  
  mod_tune <- lm(Phenotype ~ PGS_eur + PGS_eas + PGS_afr, data = tune_set)
  test_set$PGS_all <- predict(mod_tune, test_set)
  
  for (i in 1:length(pgs_models)) {
    
    formula_r2 <- as.formula(paste("Phenotype ~", pgs_model[[pgs_models[i]]]))
    r2o <- summary(lm(formula_r2, data = test_set))$adj.r.sq
    if (r2o < 0) {r2o <- 0}
    r2l <- prs_r2obs_to_r2liab(
      K = 0.05,
      P = mean(test_set$Phenotype),
      prs_r2obs = r2o
    )
    mod1 <- glm(formula_r2, data = test_set, family = "binomial")
    auc <- as.numeric(roc(response = test_set$Phenotype, predictor = predict(mod1), quiet = TRUE)$auc)
    
    
    ## r2o to auc
    n1    <- sum(test_set$Phenotype == 1); n2 <- (sum(test_set$Phenotype == 0))
    alpha <- (n1 + n2)^2 / (n1 * n2)
    d     <- (sqrt(alpha) * sqrt(r2o)) / sqrt(1 - r2o)
    auc_analytic   <- pnorm(d/sqrt(2), 0, 1)
    
    
    pgs_eval[nrow(pgs_eval) + 1, "cohort"] <- cohort
    pgs_eval[nrow(pgs_eval), "method"] <- "sbayesrc"
    pgs_eval[nrow(pgs_eval), "model"] <- pgs_models[i]
    pgs_eval[nrow(pgs_eval), "r2o"] <- r2o
    pgs_eval[nrow(pgs_eval), "r2l"] <- r2l
    pgs_eval[nrow(pgs_eval), "auc"] <- auc 
    pgs_eval[nrow(pgs_eval), "auc_analytic"] <- auc_analytic 
    
  }
}

pgs_eval_bbj <- readRDS(here("analysis/pgs/sbayesrc/taka/BBJ_yeaapoe.20250428.rds"))$pgs_eval %>%
  filter(model %in% c("PC[1-10]", "PC[1-10] + PGS_eur", "PC[1-10] + PGS_eas", "PC[1-10] + PGS_afr", "PC[1-10] + PGS_all")) %>% 
  mutate(model = case_when(
    model == "PC[1-10] + PGS_eur" ~ "PGS_eur",
    model == "PC[1-10] + PGS_eas" ~ "PGS_eas",
    model == "PC[1-10] + PGS_afr" ~ "PGS_afr",
    model == "PC[1-10] + PGS_all" ~ "PGS_all",
    .default = model
    ))

pgs_eval <- rbind(pgs_eval, pgs_eval_bbj)

pgs_eval_lst[["sbayesrc_mult"]] <- pgs_eval
  
saveRDS(pgs_eval_lst, file = here("analysis/pgs/pgs_eval.rds"))
