# =============================================================
#
#  suppFigure9_apoe_effectSize.R
#  Plot APOE-e4 effect sizes across ancestries and phenotype
#  definitions (case-control, proxy, combined) to examine
#  consistency of the strongest AD genetic risk factor.
#
# =============================================================

#### 0. Set-up ####

library(here)
library(dplyr)
library(ggplot2)

#### 1. Load and prepare APOE effect size data ####

apoe <- readRDS(here("analysis/sensitivity_analysis/apoe_effectSize", "apoe_effectSize.rds"))

# Combine results from all phenotype definitions
df <- bind_rows(
  apoe$combined     %>% mutate(phenotype = "combined"),
  apoe$case_control %>% mutate(phenotype = "case_control"),
  apoe$proxy        %>% mutate(phenotype = "proxy")
) %>%
  mutate(
    # Calculate 95% confidence intervals
    ci_lo = beta - 1.96 * standard_error,
    ci_hi = beta + 1.96 * standard_error,
    phenotype = factor(phenotype, levels = c("case_control", "proxy", "combined"))
  ) %>%
  filter(anc != "amr")  # AMR does not contain the APOE-e4 locus

# Order ancestries by mean effect size
anc_order <- df %>%
  group_by(anc) %>%
  summarise(mean_beta = mean(beta)) %>%
  arrange(mean_beta) %>%
  pull(anc)

df <- df %>% mutate(anc = factor(anc, levels = anc_order))

#### 2. Create faceted point-range plot ####

p1 <- ggplot(df, aes(x = phenotype, y = beta, ymin = ci_lo, ymax = ci_hi)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +  # Null effect reference line
  geom_pointrange() +  # Points with 95% CI error bars
  facet_wrap(~ anc, nrow = 1) +  # Separate panel for each ancestry
  labs(
    x = NULL,
    y = "Beta (log-OR)",
    title = "APOE-e4 (19:45411941:C:T) effect sizes by ancestry and phenotype"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

#### 3. Save figure ####

ggsave(
  plot = p1,
  filename = here("analysis/figures", "apoe_effects.png"),
  width = 12,
  height = 4,
  dpi = 600,
  device = png,
  type = "cairo"
)