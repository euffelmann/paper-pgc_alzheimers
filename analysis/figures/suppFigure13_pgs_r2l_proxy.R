# =============================================================
#
#  suppFigure13_pgs_r2l_proxy.R
#  Compare polygenic score (PGS) performance when trained on
#  proxy phenotype summary statistics versus full summary
#  statistics (proxy + case-control combined).
#
#  Tests whether PGS derived from proxy phenotypes alone
#  (family history, parental dementia) can predict clinical
#  AD diagnosis with comparable accuracy to PGS derived from
#  full case-control + proxy data.
#
# =============================================================

#### 0. Set-up ####
library(data.table)
library(dplyr)
library(here)
library(ggplot2)
library(pROC)
library(forcats)
source(here("R/r2l_r2o.R"))


#### 1. Load and prepare PGS evaluation data ####

# Define prediction models to evaluate
pgs_models <- c(
  "PC[1-10]",         # Baseline: ancestry PCs only
  "APOE",             # APOE genotype only
  "PGS",              # Polygenic score only
  "APOE + PGS",       # Combined APOE + PGS
  "APOE + PGS + sex"  # Full model with sex
)

# Load PGS performance when trained on proxy phenotype data only
pgs_eval_proxy <- readRDS("analysis/pgs/proxy_prediction/proxy_pgs_eval.rds")[["sbayesrc"]] %>%
  filter(model != "APOE + PGS + sex" & model != "APOE + PGS") %>%  # Remove complex models
  group_by(cohort) %>%
  mutate(
    r2l_base = r2l[model == "PC[1-10]"],  # Baseline R2 (PCs only)
    r2l_delta = r2l - r2l_base             # Incremental R2 from adding PGS
  ) %>%
  ungroup() %>%
  filter(model == "PGS") %>%  # Keep PGS-only model
  mutate(sst = "proxy") %>%   # Label as proxy-trained
  select(cohort, r2l_delta, sst)

# Load PGS performance when trained on full data (proxy + case-control)
pgs_eval_full <- readRDS("analysis/pgs/pgs_eval.rds")[["sbayesrc_mult"]] %>%
  group_by(cohort) %>%
  mutate(
    r2l_base = r2l[model == "PC[1-10]"],  # Baseline R2 (PCs only)
    r2l_delta = r2l - r2l_base             # Incremental R2 from adding PGS
  ) %>%
  ungroup() %>%
  # Assign ancestry labels to cohorts
  mutate(ancestry = case_when(
    cohort %in% c("xstsa", "demge", "twing", "xukbb_eur_casec", "gothe") ~ "EUR",
    cohort == "xukbb_afr_casec" ~ "AFR",
    cohort == "BBJ" ~ "EAS",
    cohort == "xukbb_sas_casec" ~ "SAS"
  )) %>%
  # Filter to European ancestry PGS in European cohorts
  filter(model == "PGS_eur", cohort %in% c("xstsa", "gothe", "twing", "demge")) %>%
  mutate(sst = "full") %>%  # Label as full-data trained
  select(cohort, r2l_delta, sst)

# Combine proxy and full PGS results
pgs_eval <- rbind(pgs_eval_proxy, pgs_eval_full)

#### 2. Create bar plot comparing proxy vs full PGS performance ####

dodge <- 0.7

p1 <- pgs_eval %>%
  mutate(
    sst = factor(sst, levels = c("full", "proxy")),
    # Create readable cohort labels
    cohort_label = case_when(
      cohort == "xstsa" ~ "STSA",
      cohort == "demge" ~ "DemGene",
      cohort == "gothe" ~ "Gothenburg",
      cohort == "twing" ~ "TwinGene"
      )) %>%
  ggplot(aes(y = r2l_delta, x = reorder(cohort_label, -r2l_delta), fill = sst)) +
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
    values = c("#ece7f2", "#2b8cbe"),
    labels = c(
      "proxy" = expression(PGS[proxy]),   # PGS trained on proxy data
      "full" = expression(PGS[full])),    # PGS trained on full data
    guide = guide_legend(title = "Predictor")
  ) +
  guides(fill = guide_legend(title = "Predictor")) +
  xlab("Cohort") +
  ylab(expression("Incremental" ~ italic(R)[liability]^2)) +
  theme_classic()

#### 3. Save figure ####

ggsave(
  plot = p1,
  filename = paste0("analysis/figures/sbayesrc_r2l_increment_proxy.png"),
  device = "png",
  width = 8,
  height = 4
)
