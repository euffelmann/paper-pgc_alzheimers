# =============================================================
#
#  suppFigure6_7_heatmap.R
#  Create heatmaps showing p-values for genome-wide significant
#  loci across different ancestry groups and analysis strata.
#
#  Two heatmaps are generated:
#    - Novel loci (not previously reported)
#    - Known loci (replications from prior studies)
#
# =============================================================

#### 0. Set-up ####

library(data.table)
library(ggplot2)
library(patchwork)

#### 1. Load and prepare locus data ####

# Read supplementary data table containing p-values for all loci
# across multiple ancestry × stratum combinations
a <- fread("SuppData5c.txt")

# Simplify gene labels: keep first two genes only
a$Gene <- gsub("/.*", "/..", a$Gene)

#### 2. Prepare data for novel loci heatmap ####

# Filter to novel loci only
n <- a[a$Novel == 1, ]
n <- n[, -c("chromosome", "start", "end", "size_bp", "Novel")]

# Convert p-value columns to numeric
n[, 3:22] <- lapply(n[, 3:22], as.numeric)

# Reshape data to long format for ggplot
nl <- melt(n, id.vars = c("Gene", "locus"), measure.vars = colnames(n)[-c(1, 2)])
nl$p <- -log10(nl$value)

# Order by locus number and create display labels
nl <- nl[order(nl$locus), ]
nl$id <- paste("Locus", nl$locus, nl$Gene, sep = " ")
nl$id <- factor(nl$id, levels = unique(nl$id))

# Round p-values for display
nl$value <- signif(nl$value, digits = 3)

# Assign significance categories
nl$sig <- "none"  # P > 1e-5
nl$sig[nl$value < 1e-5] <- "sug"   # Suggestive: P < 1e-5
nl$sig[nl$value < 5e-8] <- "sig"   # Genome-wide significant: P < 5e-8

# Format p-values as text for heatmap display
nl$value <- as.character(nl$value)
nl$value[is.na(nl$value)] <- "NA"

#### 3. Create heatmap for novel loci ####

ggheatmap <- ggplot(nl, aes(variable, id, fill = sig)) +
  scale_y_discrete(limits = rev) +  # Reverse y-axis (loci top to bottom)
  scale_x_discrete(position = "top") +  # Put ancestry labels on top
  geom_tile(color = "black") +  # Create tiles with black borders
  scale_fill_manual(
    values = c("white", "#FF9900", "#FF3333"),
    breaks = c("none", "sug", "sig"),
    labels = c(">1e-5", "Suggestive (<1e-5)", "Significant (<5e-8)")
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(vjust = 0, size = 7, hjust = 0.1, angle = 45),
        axis.title.x = element_blank()) +
  theme(axis.text.y = element_text(vjust = 0.4, size = 7, hjust = 1.05),
        axis.title.y = element_blank()) +
  geom_text(aes(variable, id, label = value), color = "black", size = 1.8) +  # Add p-values as text
  guides(fill = "none") +
  theme(plot.background = element_rect(fill = 'white', colour = 'white'),
        axis.line = element_line(colour = "black"))

ggheatmap
ggsave("novelloci.png", ggheatmap, scale = 1.3, height = 200, width = 160,
       dpi = 300, units = "mm", device = 'png')

#### 4. Prepare data for known loci heatmap ####

# Filter to known/replicated loci only
n <- a[a$Novel == 0, ]
n <- n[, -c("chromosome", "start", "end", "size_bp", "Novel")]

# Convert p-value columns to numeric
n[, 3:22] <- lapply(n[, 3:22], as.numeric)

# Reshape data to long format
nl <- melt(n, id.vars = c("Gene", "locus"), measure.vars = colnames(n)[-c(1, 2)])
nl$p <- -log10(nl$value)

# Order and label
nl <- nl[order(nl$locus), ]
nl$id <- paste("Locus", nl$locus, nl$Gene, sep = " ")
nl$id <- factor(nl$id, levels = unique(nl$id))

# Round and categorize p-values
nl$value <- signif(nl$value, digits = 3)
nl$sig <- "none"
nl$sig[nl$value < 1e-5] <- "sug"
nl$sig[nl$value < 5e-8] <- "sig"

#### 5. Create heatmap for known loci ####

ggheatmap <- ggplot(nl, aes(variable, id, fill = sig)) +
  scale_y_discrete(limits = rev) +
  scale_x_discrete(position = "top") +
  geom_tile(color = "black") +
  scale_fill_manual(
    values = c("white", "#FF9900", "#FF3333"),
    breaks = c("none", "sug", "sig"),
    labels = c(">1e-5", "Suggestive (<1e-5)", "Significant (<5e-8)")
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(vjust = 0, size = 7, hjust = 0.1, angle = 45),
        axis.title.x = element_blank()) +
  theme(axis.text.y = element_text(vjust = 0.4, size = 7, hjust = 1.05),
        axis.title.y = element_blank()) +
  geom_text(aes(variable, id, label = value), color = "black", size = 1.8) +
  guides(fill = "none") +
  theme(plot.background = element_rect(fill = 'white', colour = 'white'),
        axis.line = element_line(colour = "black"))

ggheatmap
ggsave("knownloci.png", ggheatmap, scale = 1.3, height = 200, width = 160,
       dpi = 300, units = "mm", device = 'png')
