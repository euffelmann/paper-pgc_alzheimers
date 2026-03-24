library(here)
library(dplyr)
library(ggplot2)

apoe <- readRDS(here("analysis/sensitivity_analysis/apoe_effectSize", "apoe_effectSize.rds"))

df <- bind_rows(
  apoe$combined     %>% mutate(phenotype = "combined"),
  apoe$case_control %>% mutate(phenotype = "case_control"),
  apoe$proxy        %>% mutate(phenotype = "proxy")
) %>%
  mutate(
    ci_lo = beta - 1.96 * standard_error,
    ci_hi = beta + 1.96 * standard_error,
    phenotype = factor(phenotype, levels = c("case_control", "proxy", "combined"))
  ) %>%
  filter(anc != "amr") # because it does not contain the APOE-e4 locus

anc_order <- df %>%
  group_by(anc) %>%
  summarise(mean_beta = mean(beta)) %>%
  arrange(mean_beta) %>%
  pull(anc)

df <- df %>% mutate(anc = factor(anc, levels = anc_order))

p1 <- ggplot(df, aes(x = phenotype, y = beta, ymin = ci_lo, ymax = ci_hi)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange() +
  facet_wrap(~ anc, nrow = 1) +
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

ggsave(
  plot = p1,
  filename = here("analysis/figures", "apoe_effects.png"),
  width = 12,
  height = 4,
  dpi = 600,
  device = png,
  type = "cairo"
)