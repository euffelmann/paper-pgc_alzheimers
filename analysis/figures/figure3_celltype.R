# =============================================================
#
#  figure3_celltype.R
#  Create bar plot showing cell type enrichment analysis results
#  from MAGMA gene-set analysis. Shows directionality (up/down)
#  of enrichment for different brain cell types.
#
# =============================================================

#### 0. Set-up ####

library(data.table)
library(ggplot2)
library(svglite)

#### 1. Load and prepare cell type enrichment data ####

# Read supplementary data table containing cell type enrichment results
a <- fread("SuppData16.txt")

# Convert p-values to -log10 scale for visualization
a$ltp <- -log10(a$p_value)

# Add directionality: negative values indicate downregulation
a$ltp[a$Direction == "Down"] <- -1 * a$ltp[a$Direction == "Down"]

#### 2. Categorize and clean cell type labels ####

# Assign broad cell type categories
a$celltype <- "Neuron"  # Default category
a$celltype[a$FULL_NAME == "Microglia"] <- "Glial"
a$celltype[a$FULL_NAME == "Oligodendrocyte_Precursor_Cell"] <- "Glial"
a$celltype[a$FULL_NAME == "Endothelial_Cell"] <- "HEndothelial"
a$celltype[a$FULL_NAME == "Astrocytes"] <- "Glial"
a$celltype[a$FULL_NAME == "Oligodendrocytes"] <- "Glial"

# Shorten specific cell type names for plotting
a$FULL_NAME[a$FULL_NAME == "Oligodendrocyte_Precursor_Cell"] <- "OPC"

#### 3. Filter and order data ####

# Sort by p-value and keep most significant direction per cell type
a <- a[order(a$p_value), ]
a <- a[!duplicated(a$FULL_NAME), ]  # Keep first (most significant) entry per cell type
a <- a[order(a$celltype, a$p_value), ]

# Clean up cell type names for display
a$FULL_NAME <- gsub("_", " ", a$FULL_NAME)
a$FULL_NAME <- gsub("Neuron", "", a$FULL_NAME)
a$FULL_NAME <- gsub("Cell", "", a$FULL_NAME)
a$FULL_NAME <- factor(a$FULL_NAME, levels = a$FULL_NAME)

# Finalize cell type categories and ordering
a$celltype[a$celltype == "HEndothelial"] <- "Endothelial"
a$celltype <- factor(a$celltype, levels = c("Glial", "Endothelial", "Neuron"))

#### 4. Create bar plot ####

p <- ggplot(a, aes(y = ltp, x = FULL_NAME, fill = celltype)) +
  geom_bar(stat = "identity", color = "black", position = position_dodge()) +
  ylab("-log10(P) * Direction") +
  xlab("Cell type") +
  theme(axis.text.y = element_text(vjust = 1, size = 9, hjust = 1)) +
  guides(fill = guide_legend(title = "")) +
  # Add Bonferroni-corrected significance threshold lines (0.05/42 tests)
  geom_hline(yintercept = -log10(0.05/42), linetype = "longdash",
             color = "black", size = 1) +
  geom_hline(yintercept = -log10(0.05/42) * -1, linetype = "longdash",
             color = "black", size = 1) +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  theme_bw() +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid = element_line(linetype = "blank")) +
  theme(axis.text.x = element_text(vjust = 1, size = 9, hjust = 1, color = "black", angle = 45)) +
  theme(axis.text.y = element_text(vjust = 1, size = 9, hjust = 1, color = "black")) +
  scale_fill_manual(breaks = c("Glial", "Neuron", "Endothelial"),
                    values = c("#83AFCF", "#ECE7F2", "#AAD5B3"))
p

#### 5. Save figure ####

svglite("Figure3.svg")
p
dev.off()