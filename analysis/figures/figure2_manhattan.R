# =============================================================
#
#  figure2_manhattan.R
#  Generate circular Manhattan plots for Figure 2 of the
#  PGC-ALZ3 manuscript showing genome-wide association results
#  across all chromosomes in a circular layout.
#
#  Two versions are created:
#    - Full scale (ylim 0-120)
#    - Zoomed view (ylim 0-75)
#
# =============================================================

#### 0. Set-up ####

library(here)
library(data.table)
library(dplyr)
library(CMplot)

#### 1. Load and filter summary statistics ####

# Load main combined meta-analysis results
temp <- fread(here("data/sumstats/meta/main/combined/main_combined_all_neff0.6_nsumstats1.txt.gz"))

# Prepare data for circular Manhattan plot
# - Convert p-values to numeric and handle zero values
# - Filter to suggestive significance (p < 5e-5) for visualization
# - Rename columns to match CMplot requirements
temp2 <- temp %>%
  mutate(
    p_value = as.numeric(p_value),
    p_value = ifelse(p_value == 0, 1e-320, p_value)) %>%  # Handle p=0 (set to machine minimum)
  filter(p_value < 5e-5) %>%  # Only plot suggestive hits for clarity
  select(variant_id, chromosome, base_pair_location, p_value) %>%
  rename(
    SNP = variant_id,
    pos = base_pair_location
  )

#### 2. Generate circular Manhattan plots ####

# Set output directory
setwd(here("analysis/figures"))

# Version 1: Full scale circular Manhattan plot (y-axis 0-120)
CMplot(
  as.data.frame(temp2),
  plot.type = "c",              # Circular plot
  r = 0.9,                      # Radius
  outward = TRUE,               # Plot points outward from center
  cir.chr.h = .1,               # Chromosome label height
  chr.den.col = "orange",       # Chromosome density color
  file = "pdf",                 # Output format
  col = c("#3690c0", "#bdbdbd"),# Alternating chromosome colors
  ylim = c(0, 120),             # Y-axis limit (full scale)
  cex = 0.3,                    # Point size
  dpi = 300,                    # Resolution
  chr.labels = seq(1, 22)       # Show autosomes only
)

# Version 2: Zoomed circular Manhattan plot (y-axis 0-75)
# This provides better resolution for lower -log10(p) values
temp3 <- temp2 %>%
  rename(p_value_75 = p_value)  # Rename for different y-axis scaling

CMplot(
  as.data.frame(temp3),
  plot.type = "c",
  r = 0.9,
  outward = TRUE,
  cir.chr.h = .1,
  chr.den.col = "orange",
  file = "pdf",
  col = c("#3690c0", "#fc8d59"),# Different color scheme for distinction
  ylim = c(0, 75),              # Y-axis limit (zoomed)
  cex = 0.3,
  dpi = 300,
  chr.labels = seq(1, 22)
)
