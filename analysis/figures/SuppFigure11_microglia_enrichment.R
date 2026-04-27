# =============================================================
#
#  SuppFigure11_microglia_enrichment.R
#  Create bar plots showing enrichment of AD GWAS signals in
#  microglia cell types and microglia functional states across
#  different brain regions and datasets.
#
#  Two plots are generated:
#    - Microglia enrichment by dataset and brain region
#    - Microglia functional state enrichment
#
# =============================================================

#### 0. Set-up ####

library(data.table)
library(ggplot2)
library(svglite)

#### 1. Microglia enrichment by dataset and brain region ####

# Read supplementary data table containing microglia enrichment results
a <- fread("SuppData14.txt")
# Clean up dataset names
a$Dataset <- gsub("_Human_2022_level2", "", a$Dataset)
a$Dataset <- gsub("_level2_Adult", "", a$Dataset)
a$Dataset <- gsub("_level2", "", a$Dataset)
a$Dataset <- gsub("^.*_Sil", "Sil", a$Dataset)

# Assign brain region labels
a$Region <- "Hippocampus"  # Default region
a$Region[a$Dataset == "GSE168408_Human_Prefrontal_Cortex"] <- "Prefrontal cortex"
a$Region[a$Dataset == "PsychENCODE_Adult"] <- "Prefrontal cortex, Temporal cortex, cerebellum"
a$Region[a$Dataset == "Allen_Human_MTG"] <- "Middle temporal gyrus"

# Create bar plot for microglia enrichment
p <- ggplot(a, aes(x = -log10(p_value), y = reorder(Dataset, -p_value), fill = Region)) + 
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  ylab("-log10(P)") +
  xlab("Microglia") +
  theme(axis.text.y = element_text(vjust = 1, size = 9, hjust = 1)) +
  guides(fill = guide_legend(title = "")) +
  # Add Bonferroni-corrected significance threshold (0.05/336 tests)
  geom_vline(xintercept = -log10(0.05/336), linetype = "longdash",
             color = "black", size = 1) +
  theme_bw() +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid = element_line(linetype = "blank")) +
  theme(axis.text.x = element_text(vjust = 1, size = 9, hjust = 1, color = "black", angle = 45)) +
  theme(axis.text.y = element_text(vjust = 1, size = 9, hjust = 1, color = "black")) +
  scale_fill_manual(breaks = c("Hippocampus", "Prefrontal cortex",
                                "Prefrontal cortex, Temporal cortex, cerebellum",
                                "Middle temporal gyrus"),
                    values = c("#084081", "#4eb3d3", "#a8ddb5", "#f7fcf0"))

p
svglite("SuppFigure6a.svg")
p
dev.off()


#### 2. Microglia functional state enrichment ####

# Read supplementary data table containing microglia functional state results
a <- fread("SuppData13.txt")
a$label <- gsub("_", " ", a$label)

# Create bar plot for microglia functional states
p <- ggplot(a, aes(x = -log10(p_value), y = reorder(label, -p_value), fill = label)) + 
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  xlab("-log10(P)") +
  ylab("") +
  theme(axis.text.x = element_text(vjust = 1, size = 9, hjust = 1)) +
  guides(fill = "none") +
  # Add Bonferroni-corrected significance threshold (0.05/12 tests)
  geom_vline(xintercept = -log10(0.05/12), linetype = "longdash",
             color = "black", size = 1) +
  theme_bw() +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid = element_line(linetype = "blank")) +
  theme(axis.text.x = element_text(vjust = 1, size = 9, hjust = 1, color = "black")) +
  theme(axis.text.y = element_text(vjust = 1, size = 9, hjust = 1, color = "black")) +
  scale_fill_manual(breaks = c("Inflammatory I", "Inflammatory II", "Inflammatory III",
                                "Homeostatic", "Glycolytic", "Stress signature",
                                "Lipid processing", "Antiviral", "Neuronal surveillance",
                                "Phagocytic", "Cycling", "Ribosome biogenesis"),
                    values = c("#213564", "#213564", "#213564", "#0B519C", "#4392C6",
                               "#2272B6", "#6BAED6", "#9FCAE2", "#C6DBEF", "#DEEBF7",
                               "#F8FCFF", "#FBFDFF"))
p

svglite("SuppFigure6b.svg")
p
dev.off()