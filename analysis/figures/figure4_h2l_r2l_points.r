# =============================================================
#
#  figure4_h2l_r2l_points.R
#  Two-panel figure combining SNP heritability on the liability
#  scale (h2l) and incremental R2l from polygenic scores (PGS).
#
#  Panel a: h2l estimates from LDSC and SBayesRC across three
#  ancestry groups (EUR, AFR, EAS), shown as dodged points with
#  95% CI error bars.
#
#  Panel b: Paired points showing incremental R2l (PGS vs.
#  PC-only baseline) per cohort, split by cohort type (clinical
#  vs. biobank) and coloured by PGS model (EUR-trained vs.
#  multi-ancestry).
#
# =============================================================

#### 0. Set-up ####

library(dplyr)
library(ggplot2)
library(here)
library(data.table)
library(openxlsx)
library(patchwork)

# Shared plot constants
DODGE        <- 0.5
FIG_W_WIDE   <- 8
FIG_W_NARROW <- 6
FIG_H_PANEL  <- 4
FIG_H_DOUBLE <- 6


#### 1. Load and prepare h2l data ####

# LDSC: keep main case-control analysis for EUR/AFR/EAS
ldsc_h2 <- fread(here("analysis", "ldsc", "ldsc_h2.txt")) %>%
  filter(
    proxy_casec == "case_control",
    analysis    == "main",
    anc         %in% c("eur", "afr", "eas")
  ) %>%
  mutate(method = "ldsc") %>%
  select(anc, h2, h2_se, method)

# SBayesRC: neff0.6 datasets correspond to sample size filtering
# used in the main analysis; hsq is the heritability parameter
sbay_h2 <- read.xlsx(
  here("analysis", "sbayesrc", "h2_estimates_8Apr2026.xlsx"),
  sheet = "8Apr2026"
) %>%
  filter(
    Dataset %in% c("main_case_control_eur_neff0.6_nsumstats1",
                   "main_case_control_eas_neff0.6_nsumstats1",
                   "main_case_control_afr_neff0.6_nsumstats1"),
    h2.parameter == "hsq") %>%
  mutate(
    anc    = tolower(Ancestry),
    method = "sbayesrc"
  ) %>%
  rename(h2 = Estimate, h2_se = Posterior.SE) %>%
  select(anc, h2, h2_se, method)

h2 <- rbind(ldsc_h2, sbay_h2)


#### 02. Figure: h2l (LDSC vs. SBayesRC) ####

p1 <- h2 %>%
  mutate(
    anc = toupper(anc),
    anc = factor(anc, levels = c("EUR", "AFR", "EAS"))
  ) %>%
  ggplot(aes(y = h2, x = anc, fill = method)) +
  geom_errorbar(
    aes(ymin = h2 - h2_se * 1.96, ymax = h2 + h2_se * 1.96),
    position = position_dodge(width = DODGE),
    width    = 0.15,
    colour   = "black"
  ) +
  geom_point(
    position = position_dodge(width = DODGE),
    size     = 2.8,
    shape    = 21,
    stroke   = 0.7,
    colour   = "black"
  ) +
  scale_fill_manual(
    values = c("ldsc" = "#8c96c6", "sbayesrc" = "#8856a7"),
    labels = c("ldsc" = "LDSC", "sbayesrc" = "SBayesRC"),
    guide  = guide_legend(title = "Method")
  ) +
  coord_cartesian(ylim = c(0, 0.25)) +
  xlab("Ancestry") +
  ylab(expression(italic(h)[liability]^2)) +
  theme_classic()


#### 03. Prepare R2l by cohort type ####

cohort_type <- c(
  "demge"           = "Clinical samples",
  "gothe"           = "Clinical samples",
  "twing"           = "Biobanks",
  "xukbb"           = "Biobanks",
  "BBJ"             = "Biobanks",
  "xukbb_afr_casec" = "Biobanks",
  "xukbb_eur_casec" = "Biobanks",
  "xukbb_sas_casec" = "Biobanks",
  "xstsa"           = "unknown"
)

pgs_eval_swedish_ukb <- readRDS(here("analysis/pgs/sbayesrc/RAP/pgs_eval.rds"))[["sbayesrc_mult"]] %>%
  filter(method == "sbayesrc_new")

pgs_eval_bbj <- readRDS(here("analysis/pgs/sbayesrc/taka/BBJ_yeaapoe.20260410.rds"))[["pgs_eval"]] %>%
  filter(model %in% c(
    "PC[1-10]",
    "PC[1-10] + PGS_eur",
    "PC[1-10] + PGS_eas",
    "PC[1-10] + PGS_afr",
    "PC[1-10] + PGS_all"
  )) %>%
  mutate(model = stringr::str_remove(model, "PC\\[1-10\\] \\+ "))

pgs_eval <- rbind(pgs_eval_swedish_ukb, pgs_eval_bbj) %>%
  group_by(cohort) %>%
  mutate(
    r2l_base  = r2l[model == "PC[1-10]"],
    r2l_delta = r2l - r2l_base
  ) %>%
  ungroup() %>%
  mutate(ancestry = case_when(
    cohort %in% c("xstsa", "demge", "twing", "xukbb_eur_casec", "gothe") ~ "EUR",
    cohort == "xukbb_afr_casec"                                           ~ "AFR",
    cohort == "BBJ"                                                        ~ "EAS",
    cohort == "xukbb_sas_casec"                                            ~ "SAS"
  ))

pgs_points <- pgs_eval %>%
  mutate(pheno = recode(cohort, !!!cohort_type, .default = NA_character_)) %>%
  filter(
    model    %in% c("PGS_eur", "PGS_all"),
    ancestry != "SAS",
    pheno    != "unknown",
    !is.na(pheno)
  ) %>%
  mutate(
    model    = case_when(model == "PGS_eur" ~ "EUR", model == "PGS_all" ~ "MULTI"),
    model    = factor(model, levels = c("EUR", "MULTI")),
    pheno    = factor(pheno, levels = c("Clinical samples", "Biobanks"))
  ) %>%
  # Numeric x is local to its own facet (resets to 1 at the start of each panel)
  group_by(pheno) %>%
  mutate(
    ancestry = droplevels(factor(ancestry, levels = c("EUR", "AFR", "EAS"))),
    x        = as.numeric(ancestry) + ifelse(model == "EUR", -DODGE / 4, DODGE / 4)
  ) %>%
  ungroup()

## Per-facet tick/label data (because x is numeric)
x_breaks <- pgs_points %>%
  distinct(pheno, ancestry) %>%
  mutate(x_break = as.numeric(ancestry),
         label   = as.character(ancestry))


#### 04. Figure: paired points of R2l ####

p2 <- ggplot(pgs_points, aes(x = x, y = r2l_delta)) +
  geom_line(
    aes(group = cohort),
    colour    = "grey60",
    linewidth = 0.5,
    alpha     = 0.8
  ) +
  geom_point(
    aes(fill = model),
    size   = 2.5,
    shape  = 21,
    stroke = 0.6,
    colour = "black"
  ) +
  facet_grid(
    ~pheno,
    scales = "free_x",
    space  = "free_x"
  ) +
  # Use the same breaks across facets; free_x will only show those that fall
  # inside each panel's data range.
  scale_x_continuous(
    breaks = 1:3,
    labels = c("EUR", "AFR", "EAS"),
    expand = expansion(add = 0.5)
  ) +
  coord_cartesian(ylim = c(0, 0.22)) +
  scale_fill_manual(values = c("EUR" = "#2b8cbe", "MULTI" = "#ece7f2")) +
  guides(fill = guide_legend(title = "PGS Model")) +
  xlab("Cohort") +
  ylab(expression("Incremental" ~ italic(R)[liability]^2)) +
  theme_classic() +
  theme(
    # Facet strip text only, no box
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 9,
                                    margin = margin(t = 5, b = 5))
  )


#### 05. Combine with patchwork ####

p_clinic_bio <- (p1 / p2) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold"))

ggsave(here("analysis", "figures", "fig4_h2l_r2l_points.pdf"),
       plot = p_clinic_bio, width = FIG_W_NARROW, height = FIG_H_DOUBLE,
       device = pdf)