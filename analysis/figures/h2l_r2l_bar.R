#### 0. set-up ####
library(dplyr)
library(ggplot2)
library(here)
library(data.table)
library(openxlsx)

#### 1. load h2l ####
## ldsc
ldsc_h2 <- fread(here("analysis/ldsc/ldsc_h2.txt")) %>%
  filter(
    proxy_casec == "case_control",
    analysis == "main",
    anc %in% c("eur", "afr", "eas")) %>%
  mutate(method = "ldsc") %>%
  select(anc, h2, h2_se, method)

## sbayesrc
sbay_h2 <- read.xlsx(here("analysis/sbayesrc/h2_estimates_14Apr2025.xlsx"), sheet = "14Apr2025") %>%
  filter(Dataset %in% c("main_case_control_eur_neff0.6_nsumstats1",
                        "main_case_control_eas_neff0.6_nsumstats1",
                        "main_case_control_afr_neff0.6_nsumstats1")) %>%
  mutate(
    anc = tolower(Ancestry),
    method = "sbayesrc") %>%
  rename(h2 = h2.estimate, h2_se = Posterior.SE) %>%
  select(anc, h2, h2_se, method)


h2 <- rbind(ldsc_h2, sbay_h2)

dodge <- 0.4
p1 <- h2 %>%
  mutate(
    anc = toupper(anc),
    anc = factor(anc, levels = c("EUR", "AFR", "EAS"))) %>%
  ggplot(aes(y = h2, x = anc, fill = method)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = dodge),
    width = dodge,
    colour = "black",
    alpha = 0.8
  ) +
  geom_errorbar(
    aes(ymin = h2 - h2_se * 1.96, ymax = h2 + h2_se * 1.96),
    position = position_dodge(width = dodge),
    width = 0.1
  ) +
  geom_point(position = position_dodge(width = dodge),
             show.legend = FALSE) +
  scale_fill_manual(
    values = c("#8c96c6", "#8856a7"),
    labels = c(
      "ldsc" = "LDSC",
      "sbayesrc" = "SBayesRC"
    ),
    guide = guide_legend(title = "Method")
  ) +
  xlab("Ancestry") +
  ylab(expression(italic(h)[liability]^2)) +
  theme_classic()

#### 2. pgs r2l ####
pgs_models <- c(
  "PC[1-10]",
  "PGS_eur",
  "PGS_eas",
  "PGS_afr",
  "PGS_all"
)

pgs_eval <- readRDS(here("analysis/pgs/pgs_eval.rds"))[["sbayesrc_mult"]] %>%
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
  ))

dodge <- 0.7
p2 <- pgs_eval %>%
  filter(model %in% c("PGS_eur", "PGS_all")) %>%
  mutate(
    model = case_when(
      model == "PGS_eur" ~ "EUR",
      model == "PGS_all" ~ "MULTI"
    ),
    model = factor(model, levels = c("EUR", "MULTI")),
    ancestry = factor(ancestry, levels = c("EUR", "AFR","EAS", "SAS")),
    cohort_label = case_when(
      cohort == "xstsa" ~ "STSA",
      cohort == "demge" ~ "DemGene",
      cohort == "gothe" ~ "Gothenburg",
      cohort == "twing" ~ "TwinGene",
      grepl("xukbb", cohort) ~ "UK Biobank",
      cohort == "BBJ" ~ "BioBank Japan"
    )) %>%
  ggplot(aes(y = r2l_delta, x = forcats::fct_reorder(cohort_label, desc(r2l_delta), .fun='min'), fill = model)) +
  stat_summary(
    fun = "mean", 
    geom = "bar", 
    position = position_dodge(width = dodge),
    width = 0.7, 
    colour = "black",
    alpha = 0.7
  ) +
  facet_grid(~ancestry, scales = "free", space = "free") +
  ylim(0, 0.18) +
  scale_fill_manual(values = c("#ece7f2", "#2b8cbe")) +
  guides(fill=guide_legend(title="PGS Model")) +
  xlab("Cohort") +
  ylab(expression("Incremental" ~ italic(R)[liability]^2)) + 
  theme_classic()

#### 3. combine figures ####

p3 <- ggpubr::ggarrange(p2, p1, labels = c("a", "b"), ncol = 1, nrow = 2)

ggsave(
  plot = p3,
  filename = paste0("analysis/figures/h2l_r2l.png"),
  device = "png",
  width = 8,
  height = 8
)

ggsave(
  here("analysis/figures/h2l_r2l.pdf"), 
  plot = p3, width = 8, height = 8, device = cairo_pdf
  )


#### 4. h2l_r2l (clinical vs. biobank v.1) ####

cohort_summary <- read.xlsx(here("analysis/cohort_summary/cohort_summary.xlsx")) %>%
  select(cohort, neff)

pgs_eval_pheno <- pgs_eval %>%
  mutate(
    pheno = case_when(
      cohort %in% c("demge", "gothe") ~ "Clinical samples",
      cohort %in% c(
        "twing",
        "xukbb",
        "BBJ",
        "xukbb_afr_casec",
        "xukbb_eur_casec",
        "xukbb_sas_casec"
      ) ~ "Biobanks",
      cohort == "xstsa" ~ "unknown"
    )
  )

dodge <- 0.7


reference_lines <- data.frame(
  ancestry = c("EUR", "AFR", "EAS"),
  yintercept = c(0.19, 0.17, 0.18)
)

p4 <- pgs_eval_pheno %>%
  filter(
    model %in% c("PGS_eur", "PGS_all"),
    ancestry != "SAS",
    pheno != "unknown") %>%
  mutate(
    model = case_when(
      model == "PGS_eur" ~ "EUR",
      model == "PGS_all" ~ "MULTI"
    ),
    model = factor(model, levels = c("EUR", "MULTI")),
    ancestry = factor(ancestry, levels = c("EUR", "AFR","EAS", "SAS")),
    pheno = factor(pheno, levels = c("Clinical samples", "Biobanks"))
  ) %>%
  ## add neff ##
  left_join(cohort_summary, by = "cohort") %>%
  mutate(
    neff = case_when(
      cohort == "xukbb_eur_casec" ~ 9608,
      cohort == "xukbb_afr_casec" ~ 4 / ((1/110) + (1/440)),
      cohort == "BBJ" ~ 4 / ((1/404) + (1/156134)),
      .default = neff
    )
  ) %>%
  group_by(model, pheno, ancestry) %>%
  summarise(
    r2l_delta_weighted = sum(r2l_delta * neff) / sum(neff),
    r2l_delta_min = min(r2l_delta),
    r2l_delta_max = max(r2l_delta),
    .groups = "drop"
  ) %>%
  ggplot(aes(y = r2l_delta_weighted, x = forcats::fct_reorder(ancestry, desc(r2l_delta_weighted), .fun='min'), fill = model)) +
  stat_summary(
    fun = "identity", 
    geom = "bar", 
    position = position_dodge(width = dodge),
    width = dodge, 
    colour = "black",
    alpha = dodge
  ) +
  geom_errorbar(
    aes(
      ymin = r2l_delta_min,
      ymax = r2l_delta_max
    ),
    position = position_dodge(width = dodge),
    width = 0.2
  ) +
  facet_grid(~pheno, scales = "free", space = "free") +
  ylim(0, 0.2) +
  geom_hline(
    data = reference_lines,
    aes(yintercept = yintercept),
    linetype = "dashed",
    color = "#de2d26"
  ) + 
  scale_fill_manual(values = c("#2b8cbe", "#ece7f2")) +
  guides(fill=guide_legend(title="PGS Model")) +
  xlab("Cohort") +
  ylab(expression("Incremental" ~ italic(R)[liability]^2)) + 
  theme_classic()

ggsave(
  here("analysis/figures/temp/r2l_v1.png"), 
  plot = p4, width = 6, height = 4, device = "png"
)

ggsave(
  here("analysis/figures/temp/r2l_v1.pdf"), 
  plot = p4, width = 6, height = 4, device = cairo_pdf
)


#### 5. h2l_r2l (clinical vs. biobank v.2) ####

p5 <- pgs_eval_pheno %>%
  filter(
    model %in% c("PGS_eur", "PGS_all"),
    ancestry != "SAS",
    pheno != "unknown") %>%
  mutate(
    model = case_when(
      model == "PGS_eur" ~ "EUR",
      model == "PGS_all" ~ "MULTI"
    ),
    model = factor(model, levels = c("EUR", "MULTI")),
    ancestry = factor(ancestry, levels = c("EUR", "AFR","EAS", "SAS")),
    pheno = factor(pheno, levels = c("Clinical samples", "Biobanks"))
  ) %>%
  ## add neff ##
  left_join(cohort_summary, by = "cohort") %>%
  mutate(
    neff = case_when(
      cohort == "xukbb_eur_casec" ~ 9608,
      cohort == "xukbb_afr_casec" ~ 4 / ((1/110) + (1/440)),
      cohort == "BBJ" ~ 4 / ((1/404) + (1/156134)),
      .default = neff
    )
  ) %>%
  group_by(model, pheno, ancestry) %>%
  summarise(
    r2l_delta_weighted = sum(r2l_delta * neff) / sum(neff),
    .groups = "drop"
  ) %>%
  ggplot(aes(y = r2l_delta_weighted, x = forcats::fct_reorder(ancestry, desc(r2l_delta_weighted), .fun='min'), fill = model)) +
  stat_summary(
    fun = "mean", 
    geom = "bar", 
    position = position_dodge(width = dodge),
    width = dodge, 
    colour = "black",
    alpha = dodge
  ) +
  facet_grid(~pheno, scales = "free", space = "free") +
  scale_fill_manual(values = c("#2b8cbe", "#ece7f2")) +
  guides(fill=guide_legend(title="PGS Model")) +
  xlab("Cohort") +
  ylab(expression("Incremental" ~ italic(R)[liability]^2)) + 
  theme_classic()

ggsave(
  here("analysis/figures/temp/r2l_v2.png"), 
  plot = p5, width = 6, height = 4, device = "png"
)

ggsave(
  here("analysis/figures/temp/r2l_v2.pdf"), 
  plot = p5, width = 6, height = 4, device = cairo_pdf
)


#### 6. h2l_r2l (clinical vs. biobank v.3) ####

p6 <- pgs_eval_pheno %>%
  filter(
    model %in% c("PGS_eur", "PGS_all"),
    ancestry != "SAS",
    pheno != "unknown") %>%
  mutate(
    model = case_when(
      model == "PGS_eur" ~ "EUR",
      model == "PGS_all" ~ "MULTI"
    ),
    model = factor(model, levels = c("EUR", "MULTI")),
    ancestry = factor(ancestry, levels = c("EUR", "AFR","EAS", "SAS")),
    pheno = factor(pheno, levels = c("Clinical samples", "Biobanks"))
  ) %>%
  ## add neff ##
  left_join(cohort_summary, by = "cohort") %>%
  mutate(
    neff = case_when(
      cohort == "xukbb_eur_casec" ~ 9608,
      cohort == "xukbb_afr_casec" ~ 4 / ((1/110) + (1/440)),
      cohort == "BBJ" ~ 4 / ((1/404) + (1/156134)),
      .default = neff
    )
  ) %>%
  group_by(model, pheno, ancestry) %>%
  summarise(
    r2l_delta_weighted = sum(r2l_delta * neff) / sum(neff),
    r2l_delta_min = min(r2l_delta),
    r2l_delta_max = max(r2l_delta),
    .groups = "drop"
  ) %>%
  ggplot(aes(y = r2l_delta_weighted, x = forcats::fct_reorder(ancestry, desc(r2l_delta_weighted), .fun='min'), fill = model)) +
  stat_summary(
    fun = "identity", 
    geom = "bar", 
    position = position_dodge(width = dodge),
    width = dodge, 
    colour = "black",
    alpha = dodge
  ) +
  geom_errorbar(
    aes(
      ymin = r2l_delta_min,
      ymax = r2l_delta_max
    ),
    position = position_dodge(width = dodge),
    width = 0.2
  ) +
  facet_grid(~pheno, scales = "free", space = "free") +
  scale_fill_manual(values = c("#2b8cbe", "#ece7f2")) +
  guides(fill=guide_legend(title="PGS Model")) +
  xlab("Cohort") +
  ylab(expression("Incremental" ~ italic(R)[liability]^2)) + 
  theme_classic()

ggsave(
  here("analysis/figures/temp/r2l_v3.png"), 
  plot = p6, width = 6, height = 4, device = "png"
)

ggsave(
  here("analysis/figures/temp/r2l_v3.pdf"), 
  plot = p6, width = 6, height = 4, device = cairo_pdf
)


#### 7. h2l_r2l (clinical vs. biobank v.3) final figure ####

p3 <- ggpubr::ggarrange(p1, p4, labels = c("a", "b"), ncol = 1, nrow = 2)

ggsave(
  plot = p3,
  filename = paste0("analysis/figures/h2l_r2l.png"),
  device = "png",
  width = 8,
  height = 8
)

ggsave(
  here("analysis/figures/h2l_r2l_clinic_vs_bio.pdf"), 
  plot = p3, width = 6, height = 8, device = cairo_pdf
)
