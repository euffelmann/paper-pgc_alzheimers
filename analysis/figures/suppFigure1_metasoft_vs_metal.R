# =============================================================
#
#  suppFigure3_metasoft_vs_metal.R
#  Compare p-values from METAL (fixed-effects meta-analysis)
#  versus METASOFT RE2 (random-effects meta-analysis with
#  heterogeneity) for genome-wide significant loci.
#
# =============================================================

#### 00. Set-up ####

library(here)
library(dplyr)
library(data.table)
library(ggplot2)

#### 01. Load data ####
loci <- fread(here("analysis/risk_loci/main/main_combined_all_neff0.6_nsumstats1_loci.txt"))

main <- 
  fread(
    here("data/sumstats/meta/main/combined/main_combined_all_neff0.6_nsumstats1.txt.gz"), 
    select = c("variant_id", "p_value")
    ) %>%
  filter(variant_id %in% loci$index_variant_id)

soft <- 
  fread(
    here("analysis/metasoft/output/combined_out.txt"), 
    select = c("RSID", "PVALUE_RE2")
    ) %>%
  filter(RSID %in% loci$index_variant_id)


#### 02. Prepare data and perform linear regression ####
# Join METAL and METASOFT results
m <- main %>%
  inner_join(soft, by = c("variant_id" = "RSID")) %>%
  mutate(p_metal = as.numeric(p_value), p_metasoft = as.numeric(PVALUE_RE2)) %>%
  mutate(
    # Handle zero p-values (set to machine minimum)
    p_metal    = ifelse(p_metal == 0,    1e-320, p_metal),
    p_metasoft = ifelse(p_metasoft == 0, 1e-320, p_metasoft)
  ) %>%
  mutate(
    # Convert to -log10 scale
    logp_metal    = -log10(p_metal),
    logp_metasoft = -log10(p_metasoft)
  )

# Fit linear regression
fit <- lm(logp_metasoft ~ logp_metal, data = m)
b0  <- coef(fit)[1]
b1  <- coef(fit)[2]
r2  <- summary(fit)$r.squared

eq_label <- sprintf("y = %.3f + %.3f * x  ;  R² = %.3f", b0, b1, r2)

#### 03. Create scatter plot ####

p <- ggplot(m, aes(x = logp_metal, y = logp_metasoft)) +
  geom_point(size = 2) +
  geom_abline(slope = 1,  intercept = 0,  linetype = "dashed", color = "red") +  # y=x reference line
  geom_abline(slope = b1, intercept = b0, linetype = "solid",  color = "blue") +  # Fitted regression line
  annotate("text", x = -Inf, y = Inf, label = eq_label,
           hjust = -0.5, vjust = 5, size = 3.5, color = "blue") +
  labs(
    x = expression(-log[10](P[METAL])),
    y = expression(-log[10](P[METASOFT])),
    title = "METAL vs. METASOFT p-values"
  ) +
  theme_classic()

#### 04. Save figure ####

ggsave(
  plot = p,
  filename = here("analysis/metasoft/metal_vs_metasoft.png"),
  device = png,
  width = 5,
  height = 5,
  dpi = 320
  )
