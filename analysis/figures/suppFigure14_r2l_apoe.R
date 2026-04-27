# =============================================================
#
#  suppFigure14_r2l_apoe.R
#  Evaluate polygenic score (PGS) performance with and without
#  the APOE region included. Compares incremental R2 on the
#  liability scale for:
#    - APOE genotype alone (e2 + e4 alleles)
#    - PGS excluding APOE region
#    - Full PGS (including APOE region)
#    - Combined APOE + PGS (excluding APOE)
#
#  This analysis quantifies the relative contributions of APOE
#  and genome-wide polygenic effects to AD risk prediction.
#
# =============================================================

#### 0. Set-up ####

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
library(pROC)

#### Helper functions: Convert R2 between observed and liability scales ####

# Function to convert observed-scale R2 to liability-scale R2
# Based on Lee et al. 2012 Genetic Epidemiology
prs_r2obs_to_r2liab <- function(K, P, prs_r2obs) {
  ## Lee et al. 2012 Genet Epidemiology
  t <- -qnorm(K, mean = 0, sd = 1)  # Disease threshold on liability scale
  z <- dnorm(t)                      # Height of normal distribution at threshold
  i1 <- z / K                        # Mean liability of cases (Falconer & Mackay)
  k1 <- i1 * (i1 - t)               # Reduction in variance in cases
  i0 <- -z / (1 - K)                # Mean liability of controls
  k0 <- i0 * (i0 - t)               # Reduction in variance in controls

  theta <- i1 * (P - K) / (1 - K) * (i1 * (P - K) / (1 - K) - t)  # theta parameter
  cv <- K * (1 - K) / z^2 * K * (1 - K) / (P * (1 - P))           # C parameter
  R2 <- prs_r2obs * cv / (1 + prs_r2obs * theta * cv)             # Liability-scale R2

  return(R2)
}

# Function to convert liability-scale R2 to observed-scale R2
prs_r2liab_to_r2obs <- function(K, P, prs_r2liab) {
  ## Lee et al. 2012 Genet Epidemiology
  t <- -qnorm(K, mean = 0, sd = 1)
  z <- dnorm(t)
  i1 <- z / K
  k1 <- i1 * (i1 - t)
  i0 <- -z / (1 - K)
  k0 <- i0 * (i0 - t)

  theta <- i1 * (P - K) / (1 - K) * (i1 * (P - K) / (1 - K) - t)
  cv <- K * (1 - K) / z^2 * K * (1 - K) / (P * (1 - P))

  prs_r2obs <- prs_r2liab / {cv - prs_r2liab * theta * cv}  # Observed-scale R2

  return(prs_r2obs)
}


#### 1. Download input data from DNAnexus ####

# Download phenotype/covariate files
system("dx download -r Dougs_data:/PGCALZ3/PGCALZ3_PRS/2026_04_10_pgsAnalysis_emil/testing_samples -o .")

# Download PGS score files (with and without APOE region)
system("dx download -r Dougs_data:/PGCALZ3/PGCALZ3_PRS/2026_04_10_pgsAnalysis_emil/plink_scores -o .")

# Download existing PGS evaluation results
system("dx download -r Dougs_data:/PGCALZ3/PGCALZ3_PRS/2026_04_10_pgsAnalysis_emil/pgs_eval.rds -o .")


#### 2. Evaluate PGS excluding APOE region ####
# Initialize results data frame
pgs_eval <- data.frame(
  cohort = as.character(),
  method = as.character(),
  model = as.character(),
  r2o = as.numeric(),           # Observed-scale R2
  r2l = as.numeric(),           # Liability-scale R2
  auc = as.numeric(),           # Area under ROC curve
  auc_analytic = as.numeric()   # Analytic AUC from R2
)

# Define model to test: PGS excluding APOE region
pgs_models <- c("PGS_exclAPOE")

# Regression formula: PCs + PGS (without APOE region)
pgs_model <- NULL
pgs_model[[pgs_models[1]]] <- "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + pgs_nayapoe"


# Loop through each cohort and evaluate PGS performance
for (cohort in c("xstsa", "gothe", "twing", "demge", "xukbb_eur_casec")) { 


  ### Load PGS scores
  # yeaapoe = PGS including APOE region
  # nayapoe = PGS excluding APOE region
  profile_yeaapoe <- fread(paste0("plink_scores/", cohort, "_eur_yeaapoe.profile")) %>%
    dplyr::select(FID, IID, SCORESUM) %>%
    rename(pgs_yeaapoe = SCORESUM)

  profile_nayapoe <- fread(paste0("plink_scores/", cohort, "_eur_nayapoe.profile")) %>%
    dplyr::select(FID, IID, SCORESUM) %>%
    rename(pgs_nayapoe = SCORESUM)

  ## Load phenotype and covariate file
  pheno <- fread(paste0("testing_samples/", cohort, ".cov")) %>%
    # Standardize column names across cohorts (UK Biobank uses different naming)
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
    filter(Phenotype %in% c(0, 1))  # Keep only cases (1) and controls (0)

  # Load APOE-e2 genotype dosage
  apoe2 <- fread(paste0("testing_samples/", cohort, "_apoe2.raw"))[, c(1, 2, 7)]
  colnames(apoe2)[3] <- "APOE2"

  # Merge APOE genotype with phenotype data
  pheno <- pheno %>%
    inner_join(apoe2, by = c("FID", "IID")) %>%
    filter(!is.na(APOE2))

  # For DemGene, exclude individuals < 40 years (too young for AD diagnosis)
  if (cohort == "demge") {pheno <- pheno %>% filter(age >= 40)}

  ## Merge phenotype data with PGS scores
  m <- pheno %>%
    inner_join(profile_yeaapoe, by = c("FID", "IID")) %>%
    inner_join(profile_nayapoe, by = c("FID", "IID"))

  # Evaluate PGS model
  for (i in 1:length(pgs_models)) {

    # Fit linear model to calculate adjusted R2
    formula_r2 <- as.formula(paste("Phenotype ~", pgs_model[[pgs_models[i]]]))
    r2o <- summary(lm(formula_r2, data = m))$adj.r.sq

    # Convert R2 from observed to liability scale
    r2l <- prs_r2obs_to_r2liab(
      K = 0.05,                 # Disease prevalence (5%)
      P = mean(m$Phenotype),    # Sample prevalence
      prs_r2obs = r2o
    )

    # Fit logistic regression and calculate AUC
    mod1 <- glm(formula_r2, data = m, family = "binomial")
    auc <- as.numeric(roc(response = m$Phenotype, predictor = predict(mod1), quiet = TRUE)$auc)

    ## Calculate analytic AUC from R2 (approximation)
    n1    <- sum(m$Phenotype == 1); n2 <- (sum(m$Phenotype == 0))
    alpha <- (n1 + n2)^2 / (n1 * n2)
    d     <- (sqrt(alpha) * sqrt(r2o)) / sqrt(1 - r2o)
    auc_analytic   <- pnorm(d / sqrt(2), 0, 1)

    # Store results
    pgs_eval[nrow(pgs_eval) + 1, "cohort"] <- cohort
    pgs_eval[nrow(pgs_eval), "method"] <- "sbayesrc"
    pgs_eval[nrow(pgs_eval), "model"] <- pgs_models[i]
    pgs_eval[nrow(pgs_eval), "r2o"] <- r2o
    pgs_eval[nrow(pgs_eval), "r2l"] <- r2l
    pgs_eval[nrow(pgs_eval), "auc"] <- auc
    pgs_eval[nrow(pgs_eval), "auc_analytic"] <- auc_analytic

  }
}

# Save results for PGS excluding APOE
pgs_eval_exclAPOE <- pgs_eval

# Load existing PGS evaluation results (full PGS with APOE)
temp <- readRDS("pgs_eval.rds")[["sbayesrc"]] %>%
  filter(method == "sbayesrc_new")

#### 3. Create bar plot comparing APOE and PGS contributions ####
# Combine results: full PGS + PGS excluding APOE
pgs_eval <- rbind(
  temp,
  pgs_eval_exclAPOE
)

# Define model order for plotting
pgs_models <- c(
  "PC[1-10]",           # Baseline
  "PGS_exclAPOE",       # PGS without APOE
  "APOE",               # APOE genotype only
  "PGS",                # Full PGS (with APOE)
  "APOE + PGS",         # APOE + PGS (excl. APOE)
  "APOE + PGS + sex"    # Full model with sex
)

# Calculate incremental R2 relative to baseline (PCs only)
pgs_eval <- pgs_eval %>%
  filter(model != "APOE + PGS + sex") %>%  # Exclude sex model
  group_by(cohort) %>%
  mutate(
    r2l_base = r2l[model == "PC[1-10]"],    # Baseline R2
    r2l_delta = r2l - r2l_base               # Incremental R2
  ) %>%
  ungroup()

dodge <- 0.7

# Create grouped bar plot
p1 <- pgs_eval %>%
  filter(model != "PC[1-10]") %>%  # Remove baseline from plot
  mutate(
    model = factor(model, levels = pgs_models),
    # Create readable cohort labels
    cohort_label = case_when(
      cohort == "xstsa" ~ "STSA",
      cohort == "demge" ~ "DemGene",
      cohort == "gothe" ~ "Gothenburg",
      cohort == "twing" ~ "TwinGene",
      grepl("xukbb", cohort) ~ "UK Biobank",
    )) %>%
  ggplot(aes(y = r2l_delta, x = reorder(cohort_label, -r2l_delta), fill = model)) +
  stat_summary(
    fun = "mean",
    geom = "bar",
    position = position_dodge(width = dodge),
    width = 0.7,
    colour = "black",
    alpha = 0.7
  ) +
  ylim(0, 0.19) +
  scale_fill_manual(
    values = c("#f1eef6", "#bdc9e1", "#74a9cf", "#0570b0"),
    labels = c(
      "APOE" = "APOE e2+e4",
      "PGS" = expression(PGS[full]),
      "PGS_exclAPOE" = expression(PGS["excl. APOE"]),
      "APOE + PGS" = bquote("APOE e2+e4 + " * PGS["excl. APOE"])
    ),
    guide = guide_legend(title = "Predictor")
  ) +
  guides(fill = guide_legend(title = "Predictor")) +
  xlab("Cohort") +
  ylab(expression("Incremental" ~ italic(R)[liability]^2)) +
  theme_classic()

#### 4. Save figure and upload to DNAnexus ####

ggsave(
  plot = p1,
  filename = paste0("supFig_r2l_apoe.png"),
  device = "png",
  width = 8,
  height = 4,
  type = "cairo"
)

system("dx upload supFig_r2l_apoe.png --path Dougs_data:/PGCALZ3/PGCALZ3_PRS/2026_04_10_pgsAnalysis_emil/")