# =============================================================
#
#  h2l_curve.R
#  Plot SNP heritability on the liability scale (h2l) as a
#  function of assumed disease prevalence, for LDSC and
#  SBayesRC EUR estimates. h2l is back-calculated from the
#  observed-scale h2o at each prevalence value using the
#  Lee et al. (2011) correction.
#
# =============================================================

#### 00. Set-up ####

library(dplyr)
library(ggplot2)
library(here)
library(data.table)
library(openxlsx)
source(here("R", "h2l_h2o.R"))

# Disease prevalence assumed when computing h2l from the GWAS
# (case fraction P = 0.5 is standard for case-control designs)
K_GWAS <- 0.05
P_GWAS <- 0.5


#### 01. Load h2l estimates and convert to observed scale ####

## LDSC — combined EUR main analysis
ldsc_h2l <- fread(here("analysis", "ldsc", "ldsc_h2.txt")) %>%
  filter(
    proxy_casec == "combined",
    analysis    == "main",
    anc         == "eur"
  ) %>%
  pull(h2)

ldsc_h2o <- h2l_to_h2o(K = K_GWAS, P = P_GWAS, h2l = ldsc_h2l)

## SBayesRC — EUR case-control main analysis
sbayesrc_h2l <- read.xlsx(
    here("analysis", "sbayesrc", "h2_estimates_14Apr2025.xlsx"),
    sheet = "14Apr2025"
  ) %>%
  filter(Dataset == "main_case_control_eur_neff0.6_nsumstats1") %>%
  pull(h2.estimate)

sbayesrc_h2o <- h2l_to_h2o(K = K_GWAS, P = P_GWAS, h2l = sbayesrc_h2l)


#### 02. Plot h2l curve across prevalence values ####

png(
  here("analysis", "figures", "h2l_curve.png"),
  width = 16, height = 12, units = "cm", res = 500
)

par(mar = c(5, 6, 4, 2) + 0.1)

curve(
  expr = function(K) h2o_to_h2l(K = K, P = P_GWAS, h2o = ldsc_h2o),
  from = 0.01, to = 0.35,
  col  = "#4575b4", lwd = 2,
  xlab = "Disease Prevalence",
  ylab = expression(italic(h)[italic("liability")]^2),
  main = " ",
  ylim = c(0, 0.4)
)

curve(
  expr = function(K) h2o_to_h2l(K = K, P = P_GWAS, h2o = sbayesrc_h2o),
  from = 0.01, to = 0.35,
  col  = "#d73027", lwd = 2,
  add  = TRUE
)

abline(v = K_GWAS, col = "grey", lty = 2)

legend(
  x      = 0.01, y = 0.39,
  legend = c("SBayesRC", "LDSC"),
  col    = c("#d73027", "#4575b4"),
  lty    = 1,
  cex    = 0.8
)

dev.off()
