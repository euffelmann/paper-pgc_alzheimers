# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # #     PGS prediction of AD cases based on       # # # # # # # # #
# # # # # # # #     summary statistics of proxy cases only    # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#### 0. setup ####
library(data.table)
library(dplyr)
library(here)
library(ggplot2)
library(pROC)
library(forcats)
source(here("R/r2l_r2o.R"))


#### 1. plot incremental R2 ####
pgs_models <- c(
  "PC[1-10]",
  "APOE",
  "PGS",
  "APOE + PGS",
  "APOE + PGS + sex"
)

pgs_eval_proxy <- readRDS("analysis/pgs/proxy_prediction/proxy_pgs_eval.rds")[["sbayesrc"]] %>%
  filter(model != "APOE + PGS + sex" & model != "APOE + PGS") %>%
  group_by(cohort) %>%
  mutate(
    r2l_base = r2l[model == "PC[1-10]"], 
    r2l_delta = r2l - r2l_base    
  ) %>%
  ungroup() %>%
  filter(model == "PGS") %>%
  mutate(sst = "proxy") %>%
  select(cohort, r2l_delta, sst)

pgs_eval_full <- readRDS("analysis/pgs/pgs_eval.rds")[["sbayesrc_mult"]] %>%
  group_by(cohort) %>%
  mutate(
    r2l_base = r2l[model == "PC[1-10]"], 
    r2l_delta = r2l - r2l_base    
  ) %>%
  ungroup() %>%
  mutate(ancestry = case_when(
    cohort %in% c("xstsa", "demge", "twing", "xukbb_eur_casec", "gothe") ~ "EUR",
    cohort == "xukbb_afr_casec" ~ "AFR",
    cohort == "BBJ" ~ "EAS",
    cohort == "xukbb_sas_casec" ~ "SAS"
  )) %>%
  filter(model == "PGS_eur", cohort %in% c("xstsa", "gothe", "twing", "demge")) %>%
  mutate(sst = "full") %>%
  select(cohort, r2l_delta, sst)

pgs_eval <- rbind(pgs_eval_proxy, pgs_eval_full)

dodge <- 0.7
p1 <- pgs_eval %>%
  mutate(
    sst = factor(sst, levels = c("full", "proxy")),
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
      "proxy" = expression(PGS[proxy]),
      "full" = expression(PGS[full])),
    guide = guide_legend(title = "Predictor")
  ) +
  guides(fill=guide_legend(title="Predictor")) +
  xlab("Cohort") +
  ylab(expression("Incremental" ~ italic(R)[liability]^2)) + 
  theme_classic()

ggsave(
  plot = p1,
  filename = paste0("analysis/figures/sbayesrc_r2l_increment_proxy.png"),
  device = "png",
  width = 8,
  height = 4
)
