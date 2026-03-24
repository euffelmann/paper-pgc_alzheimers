# =============================================================
#
#  pgs_eval.R
#  Evaluate PGS performance (R2 on liability scale, AUC) across
#  cohorts and models. Results are saved as a named list RDS for
#  use in downstream plotting scripts.
#
#  Table of Contents:
#    00. Set-up
#    01. SBayesRC — EUR PGS (APOE-included and -excluded)
#    02. SBayesRC-mult — EUR + EAS + AFR PGS
#
# =============================================================

#### 00. Set-up ####

library(data.table)
library(dplyr)
library(here)
library(pROC)
source(here("R", "r2l_r2o.R"))

K   <- 0.05   # population prevalence of AD
PCS <- paste(paste0("PC", 1:10), collapse = " + ")

# UKB cohorts use different column names for phenotype and PCs.
# This lookup maps UKB names to the standard names used elsewhere.
UKB_RENAME <- c(
  AD       = "Phenotype",
  pop_pc1  = "PC1",  pop_pc2  = "PC2",  pop_pc3  = "PC3",
  pop_pc4  = "PC4",  pop_pc5  = "PC5",  pop_pc6  = "PC6",
  pop_pc7  = "PC7",  pop_pc8  = "PC8",  pop_pc9  = "PC9",
  pop_pc10 = "PC10"
)

# -------------------------------------------------------------
# Helper: load and harmonise the phenotype/covariate file for
# a given cohort. Renames UKB columns, filters to cases and
# controls, and applies the DemGene age ≥ 40 criterion.
# -------------------------------------------------------------
load_pheno <- function(cohort) {
  pheno <- fread(here("analysis", "pgs", "testing_samples",
                      paste0(cohort, ".cov")))
  if (grepl("xukbb", cohort)) pheno <- rename(pheno, any_of(UKB_RENAME))
  pheno %>%
    filter(Phenotype %in% c(0, 1)) %>%
    { if (cohort == "demge") filter(., age >= 40) else . }
}

# -------------------------------------------------------------
# Helper: load APOE2 dosage and join to phenotype data frame
# -------------------------------------------------------------
load_apoe2 <- function(cohort, pheno) {
  apoe2 <- fread(here("analysis", "pgs", "testing_samples",
                      paste0(cohort, "_apoe2.raw")))[, c(1, 2, 7)]
  colnames(apoe2)[3] <- "APOE2"
  inner_join(pheno, apoe2, by = c("FID", "IID")) %>%
    filter(!is.na(APOE2))
}

# -------------------------------------------------------------
# Helper: load a single PLINK profile score file
# -------------------------------------------------------------
load_profile <- function(path, score_col) {
  fread(path) %>%
    select(FID, IID, SCORESUM) %>%
    rename({{ score_col }} := SCORESUM)
}

# -------------------------------------------------------------
# Helper: fit linear and logistic models, and return a one-row
# data frame of evaluation metrics for a single model formula
# -------------------------------------------------------------
eval_model <- function(data, formula_str, cohort, method, model_name) {
  formula_r2 <- as.formula(paste("Phenotype ~", formula_str))
  r2o        <- max(summary(lm(formula_r2, data = data))$adj.r.sq, 0)
  r2l        <- prs_r2obs_to_r2liab(K = K, P = mean(data$Phenotype), prs_r2obs = r2o)
  glm_fit    <- glm(formula_r2, data = data, family = "binomial")
  auc        <- as.numeric(roc(data$Phenotype, predict(glm_fit), quiet = TRUE)$auc)

  n1    <- sum(data$Phenotype == 1)
  n0    <- sum(data$Phenotype == 0)
  alpha <- (n1 + n0)^2 / (n1 * n0)
  d     <- sqrt(alpha) * sqrt(r2o) / sqrt(1 - r2o)

  data.frame(
    cohort       = cohort,
    method       = method,
    model        = model_name,
    r2o          = r2o,
    r2l          = r2l,
    auc          = auc,
    auc_analytic = pnorm(d / sqrt(2))
  )
}

pgs_eval_lst <- list()


#### 01. SBayesRC — EUR PGS (APOE-included and -excluded) ####

# yeaapoe = PGS including APOE region
# nayapoe = PGS excluding APOE region (used when APOE terms are in the model)
pgs_models_sbayesrc <- c(
  "PC[1-10]",
  "APOE",
  "PGS",
  "APOE + PGS",
  "APOE + PGS + sex"
)

pgs_formulas_sbayesrc <- c(
  "PC[1-10]"         = PCS,
  "APOE"             = paste(PCS, "APOE4 + APOE2",                      sep = " + "),
  "PGS"              = paste(PCS, "pgs_yeaapoe",                         sep = " + "),
  "APOE + PGS"       = paste(PCS, "APOE4 + APOE2 + pgs_nayapoe",        sep = " + "),
  "APOE + PGS + sex" = paste(PCS, "APOE4 + APOE2 + pgs_nayapoe + sex",  sep = " + ")
)

cohorts_sbayesrc <- c("xstsa", "gothe", "twing", "demge", "xukbb_eur_casec")

pgs_eval_lst[["sbayesrc"]] <- lapply(cohorts_sbayesrc, function(cohort) {

  m <- load_pheno(cohort) %>%
    load_apoe2(cohort, pheno = .) %>%
    inner_join(
      load_profile(
        here("analysis", "pgs", "sbayesrc", "plink_scores",
             paste0(cohort, "_eur_yeaapoe.profile")),
        "pgs_yeaapoe"
      ),
      by = c("FID", "IID")
    ) %>%
    inner_join(
      load_profile(
        here("analysis", "pgs", "sbayesrc", "plink_scores",
             paste0(cohort, "_eur_nayapoe.profile")),
        "pgs_nayapoe"
      ),
      by = c("FID", "IID")
    )

  lapply(pgs_models_sbayesrc, function(mod) {
    eval_model(m, pgs_formulas_sbayesrc[[mod]], cohort, "sbayesrc", mod)
  }) %>% bind_rows()

}) %>% bind_rows()


#### 02. SBayesRC-mult — EUR + EAS + AFR PGS ####

# A combined multi-ancestry PGS (PGS_all) is derived by fitting ancestry-
# specific scores on a tuning set (1/3 of each cohort) and predicting on
# the held-out test set (2/3).
set.seed(20250124)

pgs_models_mult <- c(
  "PC[1-10]",
  "PGS_eur",
  "PGS_eas",
  "PGS_afr",
  "PGS_all"
)

pgs_formulas_mult <- c(
  "PC[1-10]" = PCS,
  "PGS_eur"  = paste(PCS, "PGS_eur", sep = " + "),
  "PGS_eas"  = paste(PCS, "PGS_eas", sep = " + "),
  "PGS_afr"  = paste(PCS, "PGS_afr", sep = " + "),
  "PGS_all"  = paste(PCS, "PGS_all", sep = " + ")
)

cohorts_mult <- c("xstsa", "gothe", "twing", "demge",
                  "xukbb_eur_casec", "xukbb_sas_casec", "xukbb_afr_casec")

# AFR scores use a quality-filtered subset (_q0.9); other ancestries do not
profile_suffix <- function(anc) if (anc == "afr") "_q0.9" else ""

pgs_eval_mult <- lapply(cohorts_mult, function(cohort) {

  # Load per-ancestry profiles and merge with phenotype data
  profiles <- lapply(c("eur", "eas", "afr"), function(anc) {
    load_profile(
      here("analysis", "pgs", "sbayesrc", "plink_scores",
           paste0(cohort, "_", anc, profile_suffix(anc), "_yeaapoe.profile")),
      paste0("PGS_", anc)
    )
  })

  m <- load_pheno(cohort) %>%
    inner_join(profiles[[1]], by = c("FID", "IID")) %>%
    inner_join(profiles[[2]], by = c("FID", "IID")) %>%
    inner_join(profiles[[3]], by = c("FID", "IID"))

  # Tune-test split: fit PGS_all weights on 1/3, evaluate on 2/3
  tune_idx  <- sample(seq_len(nrow(m)), size = round(nrow(m) / 3))
  tune_set  <- m[tune_idx, ]
  test_set  <- m[-tune_idx, ]
  test_set$PGS_all <- predict(lm(Phenotype ~ PGS_eur + PGS_eas + PGS_afr,
                                  data = tune_set),
                               newdata = test_set)

  lapply(pgs_models_mult, function(mod) {
    eval_model(test_set, pgs_formulas_mult[[mod]], cohort, "sbayesrc_mult", mod)
  }) %>% bind_rows()

}) %>% bind_rows()

# Append BioBank Japan results (provided externally, harmonise model labels)
pgs_eval_bbj <- readRDS(
    here("analysis", "pgs", "sbayesrc", "taka", "BBJ_yeaapoe.20250428.rds")
  )$pgs_eval %>%
  filter(model %in% c(
    "PC[1-10]", "PC[1-10] + PGS_eur", "PC[1-10] + PGS_eas",
    "PC[1-10] + PGS_afr", "PC[1-10] + PGS_all"
  )) %>%
  mutate(model = recode(model,
    "PC[1-10] + PGS_eur" = "PGS_eur",
    "PC[1-10] + PGS_eas" = "PGS_eas",
    "PC[1-10] + PGS_afr" = "PGS_afr",
    "PC[1-10] + PGS_all" = "PGS_all"
  ))

pgs_eval_lst[["sbayesrc_mult"]] <- bind_rows(pgs_eval_mult, pgs_eval_bbj)

saveRDS(pgs_eval_lst, file = here("analysis", "pgs", "pgs_eval.rds"))
