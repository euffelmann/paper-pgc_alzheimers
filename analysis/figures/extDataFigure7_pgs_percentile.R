# =============================================================
#
#  suppFigure15_pgs_percentile.R
#  Create bar plot showing case prevalence across polygenic score
#  (PGS) percentiles, comparing full PGS models with those excluding
#  the APOE region.
#
# =============================================================

#### 0. Set-up ####

library(data.table)
library(ggplot2)
library(svglite)

#### 1. Load and prepare PGS percentile data ####

# Read supplementary data table containing PGS percentile and case prevalence
a <- fread("SuppData19.txt")
colnames(a)[1] <- gsub(" ", "_", colnames(a)[1])

# Clean up labels
a$PGS_Percentile[a$PGS_Percentile == "Over65"] <- "Whole Sample"
a$PGS[a$PGS_Percentile == "Whole Sample"] <- NA

# Set factor levels for ordered display
a$PGS_Percentile <- factor(a$PGS_Percentile,
                            levels = c("Whole Sample", "0-1", "0-10", "11-20", "21-30",
                                       "31-40", "41-50", "51-60", "61-70", "71-80",
                                       "81-90", "91-100", "99-100"))
a$PGS <- factor(a$PGS, levels = c("Full", "noAPOE"))

#### 2. Create bar plot ####

p <- ggplot(a, aes(y = proportion, x = PGS_Percentile, fill = PGS)) + 
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  ylab("Case prevalence") +
  xlab("PGS percentile") +
  theme(axis.text.y = element_text(vjust = 1, size = 9, hjust = 1)) +
  guides(colour = guide_legend(title = "PGS model")) +
  scale_fill_manual(breaks = c("Full", "noAPOE"),
                    values = c("#83afcf", "#ece7f2")) +
  # Add horizontal line showing overall sample prevalence
  geom_hline(yintercept = a$proportion[a$PGS_Percentile == "Whole Sample"],
             linetype = "longdash", color = "black", size = 0.5) +
  theme_bw() +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid = element_line(linetype = "blank")) +
  theme(axis.text.x = element_text(vjust = 1, size = 9, hjust = 1, color = "black", angle = 45)) +
  theme(axis.text.y = element_text(vjust = 1, size = 9, hjust = 1, color = "black"))
p

#### 3. Save figure ####

svglite("SuppFigure9.svg")
p
dev.off()