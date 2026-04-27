# =============================================================
#
#  suppFigure12_h2l_curve.R
#  Plot SNP heritability on the liability scale (h2l) as a
#  function of disease prevalence. Shows how h2l estimates
#  from LDSC and SBayesRC vary with assumed prevalence values.
#
#  The plot demonstrates sensitivity of liability-scale
#  heritability to prevalence assumptions, with a reference
#  line at 5% (typical AD prevalence in older adults).
#
# =============================================================

#### 0. Set-up ####

library(dplyr)
library(ggplot2)
library(here)
library(data.table)
library(openxlsx)
source(here("R/h2l_h2o.R"))  # Functions for converting h2 between observed and liability scales


#### 1. Load and convert heritability estimates ####

## LDSC heritability estimates
# Load LDSC results and extract h2 for European ancestry
ldsc_h2 <- fread("analysis/ldsc/ldsc_h2.txt")
ldsc_h2l <- ldsc_h2 %>%
  filter(
    proxy_casec == "combined",
    analysis == "main",
    anc == "eur") %>%
  pull(h2)

# Convert from liability scale to observed scale at K=0.05 prevalence
ldsc_h2o <- h2l_to_h2o(K = 0.05, P = 0.5, h2l = ldsc_h2l); ldsc_h2o

# Verify conversion at different prevalence values
h2o_to_h2l(K = 0.05, P = 0.5, h2o = ldsc_h2o)  # K=5%
h2o_to_h2l(K = 0.2, P = 0.5, h2o = ldsc_h2o)   # K=20%

## SBayesRC heritability estimates
# Load SBayesRC results and extract h2 for European ancestry
sbayesrc_h2l <- read.xlsx(here("analysis/sbayesrc/h2_estimates_8Apr2026.xlsx"), sheet = "8Apr2026") %>%
  filter(
    Dataset == "main_case_control_eur_neff0.6_nsumstats1",
    h2.parameter == "hsq") %>%
  pull(Estimate)

# Convert from liability scale to observed scale at K=0.05 prevalence
sbayesrc_h2o <- h2l_to_h2o(K = 0.05, P = 0.5, h2l = sbayesrc_h2l); sbayesrc_h2o

# Verify conversion at different prevalence values
h2o_to_h2l(K = 0.05, P = 0.5, h2o = sbayesrc_h2o)  # K=5%
h2o_to_h2l(K = 0.2, P = 0.5, h2o = sbayesrc_h2o)   # K=20%

#### 2. Create curve plot showing h2l as a function of prevalence ####
png(paste0(here("analysis/figures/h2l_curve.png")), width = 16, height = 12, units = "cm", res = 500)
par(mar = c(5, 6, 4, 2) + 0.1)

# Define function to convert LDSC h2o to h2l at varying prevalence K
f1 <- function(K) {
  h2o_to_h2l(K = K, P = 0.5, h2o = ldsc_h2o)
}

# Plot LDSC curve: h2l as a function of disease prevalence (K)
curve(
  expr = f1,
  from = 0.01,      # Minimum prevalence: 1%
  to = 0.35,        # Maximum prevalence: 35%
  col = "#4575b4",  # Blue for LDSC
  xlab = "Disease Prevalence",
  ylab = expression(italic(h)[italic("liability")] ^ 2),
  main = " ",
  ylim = c(0, 0.4),
  lwd = 2
)

# Define function to convert SBayesRC h2o to h2l at varying prevalence K
f2 <- function(K) {
  h2o_to_h2l(K = K, P = 0.5, h2o = sbayesrc_h2o)
}

# Add SBayesRC curve to the same plot
curve(
  expr = f2,
  from = 0.01,
  to = 0.35,
  col = "#d73027",  # Red for SBayesRC
  lwd = 2,
  add = TRUE         # Add to existing plot
)

# Add vertical line at K=0.05 (typical AD prevalence assumption)
abline(v = 0.05, col = "grey", lty = 2)

# Add legend
legend(0.01, 0.39, legend = c("SBayesRC", "LDSC"),
       col = c("#d73027", "#4575b4"), lty = 1, cex = 0.8)

dev.off()
