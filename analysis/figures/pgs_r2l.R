# =============================================================
#
#  pgs_r2l.R
#  Bar plots of incremental PGS R2 on the liability scale (R2l)
#  across cohorts, PGS models, and ancestries.
#
#  Table of Contents:
#    00. Set-up
#    01. SBayesRC — incremental R2l by predictor (p1)
#    02. SBayesRC-mult — incremental R2l by ancestry (p2)
#
# =============================================================

#### 00. Set-up ####

library(data.table)
library(dplyr)
library(here)
library(ggplot2)
library(pROC)
library(forcats)
source(here("R", "r2l_r2o.R"))

# Shared plot constants
DODGE <- 0.7
FIG_W <- 8
FIG_H <- 4

# Cohort display labels (shared across sections)
COHORT_LABELS <- c(
  "xstsa"           = "STSA",
  "demge"           = "DemGene",
  "gothe"           = "Gothenburg",
  "twing"           = "TwinGene",
  "xukbb_eur_casec" = "UK Biobank",
  "xukbb_afr_casec" = "UK Biobank",
  "xukbb_sas_casec" = "UK Biobank",
  "BBJ"             = "BioBank Japan"
)

# Cohort ancestry labels (shared across sections)
COHORT_ANCESTRY <- c(
  "xstsa"           = "EUR",
  "demge"           = "EUR",
  "twing"           = "EUR",
  "xukbb_eur_casec" = "EUR",
  "gothe"           = "EUR",
  "xukbb_afr_casec" = "AFR",
  "BBJ"             = "EAS",
  "xukbb_sas_casec" = "SAS"
)


#### 01. SBayesRC — incremental R2l by predictor (p1) ####

pgs_models_sbayesrc <- c(
  "PC[1-10]",
  "APOE",
  "PGS",
  "APOE + PGS",
  "APOE + PGS + sex"
)

pgs_eval_sbayesrc <- readRDS(here("analysis", "pgs", "pgs_eval.rds"))[["sbayesrc"]] %>%
  filter(model != "APOE + PGS + sex") %>%
  group_by(cohort) %>%
  mutate(
    r2l_base  = r2l[model == "PC[1-10]"],
    r2l_delta = r2l - r2l_base
  ) %>%
  ungroup()

p1 <- pgs_eval_sbayesrc %>%
  filter(model != "PC[1-10]") %>%
  mutate(
    model        = factor(model, levels = pgs_models_sbayesrc),
    cohort_label = recode(cohort, !!!COHORT_LABELS)
  ) %>%
  ggplot(aes(y = r2l_delta, x = reorder(cohort_label, -r2l_delta), fill = model)) +
  stat_summary(
    fun      = "mean",
    geom     = "bar",
    position = position_dodge(width = DODGE),
    width    = DODGE,
    colour   = "black",
    alpha    = 0.7
  ) +
  ylim(0, 0.19) +
  scale_fill_manual(
    values = c("#ece7f2", "#a6bddb", "#2b8cbe"),
    labels = c(
      "APOE"       = "APOE e2+e4",
      "PGS"        = expression(PGS[full]),
      "APOE + PGS" = bquote("APOE e2+e4 + " * PGS["excl. APOE"])
    ),
    guide = guide_legend(title = "Predictor")
  ) +
  xlab("Cohort") +
  ylab(expression("Incremental" ~ italic(R)[liability]^2)) +
  theme_classic()

ggsave(
  here("analysis", "figures", "sbayesrc_r2l_increment.png"),
  plot = p1, device = "png", width = FIG_W, height = FIG_H
)


#### 02. SBayesRC-mult — incremental R2l by ancestry (p2) ####

pgs_models_mult <- c(
  "PC[1-10]",
  "PGS_eur",
  "PGS_eas",
  "PGS_afr",
  "PGS_all"
)

pgs_eval_mult <- readRDS(here("analysis", "pgs", "pgs_eval.rds"))[["sbayesrc_mult"]] %>%
  group_by(cohort) %>%
  mutate(
    r2l_base  = r2l[model == "PC[1-10]"],
    r2l_delta = r2l - r2l_base
  ) %>%
  ungroup() %>%
  mutate(ancestry = recode(cohort, !!!COHORT_ANCESTRY))

p2 <- pgs_eval_mult %>%
  filter(model %in% c("PGS_eur", "PGS_all")) %>%
  mutate(
    model = case_when(
      model == "PGS_eur" ~ "EUR",
      model == "PGS_all" ~ "MULTI"
    ),
    model        = factor(model, levels = c("EUR", "MULTI")),
    ancestry     = factor(ancestry, levels = c("EUR", "AFR", "EAS", "SAS")),
    cohort_label = recode(cohort, !!!COHORT_LABELS)
  ) %>%
  ggplot(aes(
    y    = r2l_delta,
    x    = fct_reorder(cohort_label, desc(r2l_delta), .fun = "min"),
    fill = model
  )) +
  stat_summary(
    fun      = "mean",
    geom     = "bar",
    position = position_dodge(width = DODGE),
    width    = DODGE,
    colour   = "black",
    alpha    = 0.7
  ) +
  facet_grid(~ancestry, scales = "free", space = "free") +
  ylim(0, 0.18) +
  scale_fill_manual(values = c("EUR" = "#ece7f2", "MULTI" = "#2b8cbe")) +
  guides(fill = guide_legend(title = "PGS Model")) +
  xlab("Cohort") +
  ylab(expression("Incremental" ~ italic(R)[liability]^2)) +
  theme_classic()

ggsave(
  here("analysis", "figures", "sbayesrc_mult_r2l_increment.png"),
  plot = p2, device = "png", width = FIG_W, height = FIG_H
)
