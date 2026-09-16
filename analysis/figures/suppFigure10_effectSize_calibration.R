# =============================================================
#
#  suppFigure16_effectSize_calibration.R
#  Create scatter plots comparing effect sizes (beta) between
#  case-control and proxy phenotype analyses. Shows calibration
#  of effect estimates across phenotype definitions.
#
#  Two versions are created:
#    - All loci
#    - Excluding APOE locus
#
# =============================================================

#### 0. Set-up ####

library(data.table)
library(ggplot2)
library(svglite)

# Build "y=a+b*x\nR2=r" annotation label from a fitted lm
lm_label <- function(m) {
  intercept <- coef(m)[1]
  slope <- coef(m)[2]
  r2 <- summary(m)$r.squared
  paste0("y=", format(signif(intercept, 2), scientific = FALSE),
         ifelse(slope < 0, "-", "+"), format(round(abs(slope), 3), nsmall = 3), "*x",
         "\nR2=", format(round(r2, 4), nsmall = 4))
}

#### 1. Load and prepare effect size data ####

# Read supplementary data table containing lead SNP effect sizes
# from both case-control and proxy phenotype analyses
a <- fread("SuppData5.txt")

# Calculate 95% confidence intervals for both phenotypes
a$xmin <- a$ccbeta - (1.96 * a$ccstandard_error)
a$xmax <- a$ccbeta + (1.96 * a$ccstandard_error)
a$ymin <- a$proxybeta - (1.96 * a$proxystandard_error)
a$ymax <- a$proxybeta + (1.96 * a$proxystandard_error)

#### 2. Linear regression (all loci) ####

m <- lm(a$proxybeta ~ a$ccbeta)
summary(m)

#### 3. Create scatter plot with error bars (all loci) ####

p <- ggplot(data = a, aes(x = ccbeta, y = proxybeta)) +
  geom_point() +
  geom_errorbar(aes(ymin = ymin, ymax = ymax)) +  # Vertical error bars (proxy SE)
  geom_errorbarh(aes(xmin = xmin, xmax = xmax)) +  # Horizontal error bars (case-control SE)
  xlab("Case-control effect size") +
  ylab("Proxy effect size") +
  theme(axis.text.y = element_text(vjust = 1, size = 9, hjust = 1)) +
  guides(fill = guide_legend(title = "")) +
  theme_bw() +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid = element_line(linetype = "blank")) +
  theme(axis.text.x = element_text(vjust = 1, size = 9, hjust = 1, color = "black", angle = 45)) +
  theme(axis.text.y = element_text(vjust = 1, size = 9, hjust = 1, color = "black")) +
  geom_smooth(method = "lm") +  # Add linear regression line
  annotate("text", x = -0.75, y = 0.4, label = lm_label(m))
ggsave("casecontrolvsproxybeta.png", p)

#### 4. Repeat analysis excluding APOE locus ####

# Remove APOE-e4 variant (19:45411941:C:T)
a <- a[!a$LeadSNP == "19:45411941:C:T", ]
m <- lm(a$proxybeta ~ a$ccbeta)
summary(m)

#### 5. Create scatter plot excluding APOE ####

p <- ggplot(data = a, aes(x = ccbeta, y = proxybeta)) +
  geom_point() +
  geom_errorbar(aes(ymin = ymin, ymax = ymax)) +
  geom_errorbarh(aes(xmin = xmin, xmax = xmax)) +
  xlab("Case-control effect size") +
  ylab("Proxy effect size") +
  theme(axis.text.y = element_text(vjust = 1, size = 9, hjust = 1)) +
  guides(fill = guide_legend(title = "")) +
  theme_bw() +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid = element_line(linetype = "blank")) +
  theme(axis.text.x = element_text(vjust = 1, size = 9, hjust = 1, color = "black", angle = 45)) +
  theme(axis.text.y = element_text(vjust = 1, size = 9, hjust = 1, color = "black")) +
  geom_smooth(method = "lm") +
  annotate("text", x = 0, y = 0.4, label = lm_label(m))
ggsave("casecontrolvsproxybetanoAPOE.png", p)
