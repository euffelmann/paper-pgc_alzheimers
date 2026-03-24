# =============================================================
#
#  h2l_r2l_bar.R
#  Bar plots comparing SNP heritability (h2l) and incremental
#  PGS R2 on the liability scale (R2l) across ancestries,
#  methods, and cohort types (clinical vs. biobank).
#
#  Table of Contents:
#    00. Set-up
#    01. Load and plot h2l (LDSC vs. SBayesRC)
#    02. Load PGS R2l and compute incremental R2l
#    03. Plot R2l by cohort — main figure (p2)
#    04. Prepare weighted R2l by cohort type (shared base data)
#    05. Plot weighted R2l with min/max error bars (p4 / v1)
#    06. Plot weighted R2l without error bars (p5 / v2)
#    07. Plot weighted R2l with error bars, no y-limit (p6 / v3)
#    08. Combine and save final figures
#
# =============================================================

#### 00. Set-up ####

library(dplyr)
library(ggplot2)
library(here)
library(data.table)
library(openxlsx)

# Shared plot constants
DODGE        <- 0.7
FIG_W_WIDE   <- 8
FIG_W_NARROW <- 6
FIG_H_PANEL  <- 4
FIG_H_DOUBLE <- 8

#### 01. Load and plot h2l (LDSC vs. SBayesRC) ####

## LDSC heritability estimates
ldsc_h2 <- fread(here("analysis", "ldsc", "ldsc_h2.txt")) %>%
  filter(
    proxy_casec == "case_control",
    analysis    == "main",
    anc         %in% c("eur", "afr", "eas")
  ) %>%
  mutate(method = "ldsc") %>%
  select(anc, h2, h2_se, method)

## SBayesRC heritability estimates
sbay_h2 <- read.xlsx(
    here("analysis", "sbayesrc", "h2_estimates_14Apr2025.xlsx"),
    sheet = "14Apr2025"
  ) %>%
  filter(Dataset %in% c(
    "main_case_control_eur_neff0.6_nsumstats1",
    "main_case_control_eas_neff0.6_nsumstats1",
    "main_case_control_afr_neff0.6_nsumstats1"
  )) %>%
  mutate(
    anc    = tolower(Ancestry),
    method = "sbayesrc"
  ) %>%
  rename(h2 = h2.estimate, h2_se = Posterior.SE) %>%
  select(anc, h2, h2_se, method)

h2 <- rbind(ldsc_h2, sbay_h2)

## Plot
p1 <- h2 %>%
  mutate(
    anc = toupper(anc),
    anc = factor(anc, levels = c("EUR", "AFR", "EAS"))
  ) %>%
  ggplot(aes(y = h2, x = anc, fill = method)) +
  geom_bar(
    stat     = "identity",
    position = position_dodge(width = DODGE),
    width    = DODGE,
    colour   = "black",
    alpha    = 0.8
  ) +
  geom_errorbar(
    aes(ymin = h2 - h2_se * 1.96, ymax = h2 + h2_se * 1.96),
    position = position_dodge(width = DODGE),
    width    = 0.1
  ) +
  geom_point(
    position  = position_dodge(width = DODGE),
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c("ldsc" = "#8c96c6", "sbayesrc" = "#8856a7"),
    labels = c("ldsc" = "LDSC", "sbayesrc" = "SBayesRC"),
    guide  = guide_legend(title = "Method")
  ) +
  xlab("Ancestry") +
  ylab(expression(italic(h)[liability]^2)) +
  theme_classic()


#### 02. Load PGS R2l and compute incremental R2l ####

pgs_eval <- readRDS(here("analysis", "pgs", "pgs_eval.rds"))[["sbayesrc_mult"]] %>%
  group_by(cohort) %>%
  mutate(
    # Baseline R2l from a PCs-only model
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


#### 03. Plot R2l by cohort — main figure (p2) ####

## Cohort display labels
cohort_labels <- c(
  "xstsa"          = "STSA",
  "demge"          = "DemGene",
  "gothe"          = "Gothenburg",
  "twing"          = "TwinGene",
  "xukbb_eur_casec" = "UK Biobank",
  "xukbb_afr_casec" = "UK Biobank",
  "xukbb_sas_casec" = "UK Biobank",
  "BBJ"            = "BioBank Japan"
)

p2 <- pgs_eval %>%
  filter(model %in% c("PGS_eur", "PGS_all")) %>%
  mutate(
    model = case_when(
      model == "PGS_eur" ~ "EUR",
      model == "PGS_all" ~ "MULTI"
    ),
    model        = factor(model, levels = c("EUR", "MULTI")),
    ancestry     = factor(ancestry, levels = c("EUR", "AFR", "EAS", "SAS")),
    cohort_label = recode(cohort, !!!cohort_labels)
  ) %>%
  ggplot(aes(
    y = r2l_delta,
    x = forcats::fct_reorder(cohort_label, desc(r2l_delta), .fun = "min"),
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


#### 04. Prepare weighted R2l by cohort type (shared base data) ####

## Cohort-type labels
cohort_type <- c(
  "demge"          = "Clinical samples",
  "gothe"          = "Clinical samples",
  "twing"          = "Biobanks",
  "xukbb"          = "Biobanks",
  "BBJ"            = "Biobanks",
  "xukbb_afr_casec" = "Biobanks",
  "xukbb_eur_casec" = "Biobanks",
  "xukbb_sas_casec" = "Biobanks",
  "xstsa"          = "unknown"
)

## neff overrides for cohorts missing or incorrectly valued in cohort_summary
neff_overrides <- c(
  "xukbb_eur_casec" = 9608,
  "xukbb_afr_casec" = 4 / ((1 / 110)  + (1 / 440)),
  "BBJ"             = 4 / ((1 / 404)  + (1 / 156134))
)

cohort_summary <- read.xlsx(
    here("analysis", "cohort_summary", "cohort_summary.xlsx")
  ) %>%
  select(cohort, neff)

## Shared preparation pipeline for sections 05–07
pgs_eval_pheno <- pgs_eval %>%
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
    ancestry = factor(ancestry, levels = c("EUR", "AFR", "EAS", "SAS")),
    pheno    = factor(pheno, levels = c("Clinical samples", "Biobanks"))
  ) %>%
  left_join(cohort_summary, by = "cohort") %>%
  mutate(neff = ifelse(cohort %in% names(neff_overrides), neff_overrides[cohort], neff))

## Weighted summary per model x pheno x ancestry
pgs_weighted <- pgs_eval_pheno %>%
  group_by(model, pheno, ancestry) %>%
  summarise(
    r2l_delta_weighted = sum(r2l_delta * neff) / sum(neff),
    r2l_delta_min      = min(r2l_delta),
    r2l_delta_max      = max(r2l_delta),
    .groups = "drop"
  )

## Ancestry-specific h2l reference lines (for p4)
reference_lines <- data.frame(
  ancestry   = c("EUR", "AFR", "EAS"),
  yintercept = c(0.19, 0.17, 0.18)
)


#### 05. Plot weighted R2l with min/max error bars and h2l reference lines (p4 / v1) ####

p4 <- pgs_weighted %>%
  ggplot(aes(
    y    = r2l_delta_weighted,
    x    = forcats::fct_reorder(ancestry, desc(r2l_delta_weighted), .fun = "min"),
    fill = model
  )) +
  stat_summary(
    fun      = "identity",
    geom     = "bar",
    position = position_dodge(width = DODGE),
    width    = DODGE,
    colour   = "black",
    alpha    = DODGE
  ) +
  geom_errorbar(
    aes(ymin = r2l_delta_min, ymax = r2l_delta_max),
    position = position_dodge(width = DODGE),
    width    = 0.2
  ) +
  geom_hline(
    data = reference_lines,
    aes(yintercept = yintercept),
    linetype = "dashed",
    color    = "#de2d26"
  ) +
  facet_grid(~pheno, scales = "free", space = "free") +
  ylim(0, 0.2) +
  scale_fill_manual(values = c("EUR" = "#2b8cbe", "MULTI" = "#ece7f2")) +
  guides(fill = guide_legend(title = "PGS Model")) +
  xlab("Cohort") +
  ylab(expression("Incremental" ~ italic(R)[liability]^2)) +
  theme_classic()

ggsave(here("analysis", "figures", "temp", "r2l_v1.png"),
       plot = p4, width = FIG_W_NARROW, height = FIG_H_PANEL, device = "png")
ggsave(here("analysis", "figures", "temp", "r2l_v1.pdf"),
       plot = p4, width = FIG_W_NARROW, height = FIG_H_PANEL, device = cairo_pdf)


#### 06. Plot weighted R2l without error bars (p5 / v2) ####

p5 <- pgs_weighted %>%
  ggplot(aes(
    y    = r2l_delta_weighted,
    x    = forcats::fct_reorder(ancestry, desc(r2l_delta_weighted), .fun = "min"),
    fill = model
  )) +
  stat_summary(
    fun      = "mean",
    geom     = "bar",
    position = position_dodge(width = DODGE),
    width    = DODGE,
    colour   = "black",
    alpha    = DODGE
  ) +
  facet_grid(~pheno, scales = "free", space = "free") +
  scale_fill_manual(values = c("EUR" = "#2b8cbe", "MULTI" = "#ece7f2")) +
  guides(fill = guide_legend(title = "PGS Model")) +
  xlab("Cohort") +
  ylab(expression("Incremental" ~ italic(R)[liability]^2)) +
  theme_classic()

ggsave(here("analysis", "figures", "temp", "r2l_v2.png"),
       plot = p5, width = FIG_W_NARROW, height = FIG_H_PANEL, device = "png")
ggsave(here("analysis", "figures", "temp", "r2l_v2.pdf"),
       plot = p5, width = FIG_W_NARROW, height = FIG_H_PANEL, device = cairo_pdf)


#### 07. Plot weighted R2l with error bars, no y-limit or reference lines (p6 / v3) ####

p6 <- pgs_weighted %>%
  ggplot(aes(
    y    = r2l_delta_weighted,
    x    = forcats::fct_reorder(ancestry, desc(r2l_delta_weighted), .fun = "min"),
    fill = model
  )) +
  stat_summary(
    fun      = "identity",
    geom     = "bar",
    position = position_dodge(width = DODGE),
    width    = DODGE,
    colour   = "black",
    alpha    = DODGE
  ) +
  geom_errorbar(
    aes(ymin = r2l_delta_min, ymax = r2l_delta_max),
    position = position_dodge(width = DODGE),
    width    = 0.2
  ) +
  facet_grid(~pheno, scales = "free", space = "free") +
  scale_fill_manual(values = c("EUR" = "#2b8cbe", "MULTI" = "#ece7f2")) +
  guides(fill = guide_legend(title = "PGS Model")) +
  xlab("Cohort") +
  ylab(expression("Incremental" ~ italic(R)[liability]^2)) +
  theme_classic()

ggsave(here("analysis", "figures", "temp", "r2l_v3.png"),
       plot = p6, width = FIG_W_NARROW, height = FIG_H_PANEL, device = "png")
ggsave(here("analysis", "figures", "temp", "r2l_v3.pdf"),
       plot = p6, width = FIG_W_NARROW, height = FIG_H_PANEL, device = cairo_pdf)


#### 08. Combine and save final figures ####

## Main figure: h2l (p1) and per-cohort R2l (p2)
p_main <- ggpubr::ggarrange(p2, p1, labels = c("a", "b"), ncol = 1, nrow = 2)

ggsave(here("analysis", "figures", "h2l_r2l.png"),
       plot = p_main, width = FIG_W_WIDE, height = FIG_H_DOUBLE, device = "png")
ggsave(here("analysis", "figures", "h2l_r2l.pdf"),
       plot = p_main, width = FIG_W_WIDE, height = FIG_H_DOUBLE, device = cairo_pdf)

## Supplementary figure: h2l (p1) and clinical vs. biobank R2l with reference lines (p4)
p_clinic_bio <- ggpubr::ggarrange(p1, p4, labels = c("a", "b"), ncol = 1, nrow = 2)

ggsave(here("analysis", "figures", "h2l_r2l_clinic_vs_bio.png"),
       plot = p_clinic_bio, width = FIG_W_NARROW, height = FIG_H_DOUBLE, device = "png")
ggsave(here("analysis", "figures", "h2l_r2l_clinic_vs_bio.pdf"),
       plot = p_clinic_bio, width = FIG_W_NARROW, height = FIG_H_DOUBLE, device = cairo_pdf)
