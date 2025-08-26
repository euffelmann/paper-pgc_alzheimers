#### 0. set-up ####
library(dplyr)
library(ggplot2)
library(here)
library(data.table)
source(here("R/h2l_h2o.R"))


#### 1. load data ####
## ldsc
ldsc_h2 <- fread("analysis/ldsc/ldsc_h2.txt")
ldsc_h2l <- ldsc_h2 %>%
  filter(
    proxy_casec == "combined",
    analysis == "main",
    anc == "eur") %>%
  pull(h2)
ldsc_h2o <- h2l_to_h2o(K = 0.05, P = 0.5, h2l = ldsc_h2l); ldsc_h2o
h2o_to_h2l(K = 0.05, P = 0.5, h2o = ldsc_h2o)
h2o_to_h2l(K = 0.4, P = 0.5, h2o = ldsc_h2o)

## sbayesrc
sbayesrc_h2l <- read.xlsx(here("analysis/sbayesrc/h2_estimates_14Apr2025.xlsx"), sheet = "14Apr2025") %>%
  filter(Dataset == "main_case_control_eur_neff0.6_nsumstats1") %>%
  pull(h2.estimate)
sbayesrc_h2o <- h2l_to_h2o(K = 0.05, P = 0.5, h2l = sbayesrc_h2l); sbayesrc_h2o
h2o_to_h2l(K = 0.05, P = 0.5, h2o = sbayesrc_h2o)
h2o_to_h2l(K = 0.4, P = 0.5, h2o = sbayesrc_h2o)

#### 2. plot ####
png(paste0(here("analysis/figures/h2l_curve.png")), width = 16, height = 12, units = "cm", res = 500)
par(mar = c(5, 6, 4, 2) + 0.1)
f1 <- function(K) {
  h2o_to_h2l(K = K, P = 0.5, h2o = ldsc_h2o)
}
curve(
  expr = f1,
  from = 0.01,
  to = 0.35,
  col = "#4575b4",
  xlab = "Disease Prevalence",
  ylab = expression(italic(h)[italic("liability")] ^ 2),
  main = " ",
  ylim = c(0, 0.4),
  lwd = 2
)
f2 <- function(K) {
  h2o_to_h2l(K = K, P = 0.5, h2o = sbayesrc_h2o)
}
curve(
  expr = f2,
  from = 0.01,
  to = 0.35,
  col = "#d73027",
  lwd = 2,
  add = T
)
abline(v=0.05, col="grey", lty=2)
legend(0.01, 0.39, legend = c("SBayesRC", "LDSC"), col = c("#d73027", "#4575b4"), lty = c(1, 1, 5), cex = 0.8)
dev.off()
