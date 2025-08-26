#### 0. set-up ####
library(data.table)
library(dplyr)
library(ggplot2)
library(here)

#### 1. load data ####
other_rg <- fread(here("analysis/ldsc/ldsc_rg/other.txt"))

#### 2. lollipop plot
p1 <- other_rg %>%
  filter(
    pheno2 == "alzh_case_control",
    !pheno1 %in% c(
      "hdl_c_teslovich_2010",
      "ldl_c_teslovich_2010",
      "ldl_c_willer_2013",
      "hdl_c_willer_2013",
      "ldl_c_klimentidis_2020",
      "triglyceride_barton_2021",
      "total_c_richardson_2022",
      "total_c_graham_2021",
      "cad_van_der_harst_2017"
    )) %>%
  mutate(lab = case_when(
    pheno1 == "als_van_rheenen_2021" ~ "Amyotrophic Lateral Sclerosis",
    pheno1 == "pd_nalls_2019" ~ "Parkinson's Disease",
    pheno1 == "apo_a1_richardson_2022" ~ "Apolipoprotein A1",
    pheno1 == "p_tau_jansen_2022" ~ "Phosporylated Tau (CSF)",
    pheno1 == "amyloid_beta_jansen_2022" ~ "Amyloid Beta (CSF)",
    pheno1 == "hdl_c_graham_2021" ~ "HDL cholesterol",
    pheno1 == "ms_sawcer_2011" ~ "Multiple Sclerosis",
    pheno1 == "ldl_c_graham_2021" ~ "LDL cholesterol",
    pheno1 == "apo_b_richardson_2022" ~ "Apolipoprotein B",
    pheno1 == "cad_aragam_2022" ~ "Coronoary Artery Disease",
    pheno1 == "triglyceride_graham_2021" ~ "Triglyceride",
    pheno1 == "longevity_pilling_2016" ~ "Parental longevity (combined parental age at death)",
    pheno1 == "ea_davies_2016" ~ "Educational Attainment",
    pheno1 == "cog_lee_2018" ~ "Cognitive Performance",
    pheno1 == "int_savage_2018" ~ "Intelligence",
    pheno1 == "longevity_deelen_2019" ~ "Longevity (>90th survival percentile)", # https://www.ebi.ac.uk/gwas/studies/GCST008598
  )) %>%
  mutate(lab = factor(lab, levels = lab[order(rg)])) %>%
  ggplot(aes(x = rg, y = lab)) +
  geom_pointrange(aes(xmin = rg - 1.96 * se, xmax = rg + 1.96 * se), size = 0.1, linewidth = 0.7, alpha = 0.8, colour = '#2E4881') +
  geom_text(aes(label = case_when(
    p > 0.05 ~ "",
    p < 0.05 & p > 0.05 / 9 ~ "*",
    p < 0.05 / 9 ~ "**",
    TRUE ~ NA_character_)), size = 5, colour = "black") +
  geom_vline(xintercept = 0, color = "darkgrey", linetype = 2) +
  theme_minimal() +
  theme(
    panel.spacing = unit(0.1, "lines"),
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size=10),
    legend.position="top",
    legend.text = element_text(size=10),
    legend.margin=margin(0,0,0,0),
    legend.box.margin=margin(-10,-10,-10,-10),
    legend.title=element_blank(),
    strip.text.y = element_text(size = 10, angle=25),
    strip.clip = "off"
  ) +
  guides(color = guide_legend(override.aes = list(linetype = c(0, 1) ) ) ) +
  coord_cartesian(clip = "off") +
  guides(shape = guide_legend(override.aes = list(size = 1))) +
  xlab(bquote(italic(r)[g*","*global] ~ "with ALZ")) +
  ylab("")


ggsave(
  plot = p1,
  filename = paste0("analysis/figures/ldsc_rg_other.pdf"),
  device = "pdf",
  width = 7,
  height = 5
)

ggsave(
  plot = p1,
  filename = paste0("analysis/figures/ldsc_rg_other.png"),
  device = "png",
  width = 7,
  height = 5
)
