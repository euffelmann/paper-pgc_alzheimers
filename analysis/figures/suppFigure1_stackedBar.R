# =============================================================
#
#  suppFigure1_stackedBar.R
#  Create stacked bar charts showing the effective sample size
#  (neff) distribution across ancestry groups and phenotype types
#  (case-control vs proxy phenotypes).
#
#  Two versions are created:
#    - All ancestries
#    - Excluding European ancestry (for better visualization
#      of smaller ancestry groups)
#
# =============================================================

#### 0. Set-up ####

library(dplyr)
library(ggplot2)
library(svglite)
library(openxlsx)
library(here)
library(scales)

#### 1. Load and summarize cohort data ####

# Load cohort summary statistics
main <- readWorkbook(here("analysis/cohort_summary/cohort_summary.xlsx"), sheet = "main")
df <- main

# Compute effective sample size by ancestry and phenotype type
summary_df <- df %>%
  group_by(ancestry, proxy_casecontrol) %>%
  summarise(neff = sum(neff, na.rm = TRUE), .groups = "drop")

# Compute total neff per ancestry for bar labels
total_neff_df <- summary_df %>%
  group_by(ancestry) %>%
  summarise(total_neff = sum(neff), .groups = "drop")

#### 2. Order ancestries by sample size ####

# Reorder ancestry by total neff (largest to smallest)
ancestry_order <- total_neff_df %>%
  arrange(desc(total_neff)) %>%
  pull(ancestry)

# Apply factor levels to both data frames
summary_df$ancestry <- factor(toupper(as.character(summary_df$ancestry)), levels = toupper(ancestry_order))
total_neff_df$ancestry <- factor(toupper(as.character(total_neff_df$ancestry)), levels = toupper(ancestry_order))

#### 3. Create stacked bar plot (all ancestries) ####

p <- summary_df %>%
  mutate(proxy_casecontrol = case_when(
    proxy_casecontrol == "case_control" ~ "Case Control",
    proxy_casecontrol == "proxy" ~ "Proxy"
  )) %>%
ggplot(aes(x = ancestry, y = neff, fill = proxy_casecontrol, color = proxy_casecontrol)) +
  geom_bar(stat = "identity") +
  # Add total sample size labels on top of bars
  geom_text(
    data = total_neff_df,
    aes(x = ancestry, y = total_neff, label = comma(round(total_neff))),
    inherit.aes = FALSE,
    vjust = -0.5,
    fontface = "bold",
    size = 3.5
  ) +
  scale_fill_manual(
    values = c(
      "Case Control" = "#2b8cbe",
      "Proxy" = "#a6bddb"
    )
  ) +
  scale_color_manual(
    values = c(
      "Case Control" = "#2b8cbe",
      "Proxy" = "#a6bddb"
    )
  ) +
  scale_y_continuous(labels = label_comma()) +
  labs(x = "Ancestry", y = "Effective Sample Size", fill = "Type") +
  theme_classic(base_size = 14) +
  theme(legend.title = element_blank()) +
  guides(color = "none")

# Save as SVG
ggsave(here("analysis/figures/stacked_bar.svg"), plot = p, width = 6, height = 6)


#### 4. Create stacked bar plot (excluding European ancestry) ####

# This version excludes EUR to better visualize smaller ancestry groups
p <- summary_df %>%
  mutate(proxy_casecontrol = case_when(
    proxy_casecontrol == "case_control" ~ "Case Control",
    proxy_casecontrol == "proxy" ~ "Proxy"
  )) %>%
  filter(ancestry != "EUR") %>%  # Remove European ancestry
  ggplot(aes(x = ancestry, y = neff, fill = proxy_casecontrol, color = proxy_casecontrol)) +
  geom_bar(stat = "identity") +
  geom_text(
    data = total_neff_df[total_neff_df$ancestry != "EUR",],
    aes(x = ancestry, y = total_neff, label = comma(round(total_neff))),
    inherit.aes = FALSE,
    vjust = -0.5,
    fontface = "bold",
    size = 3.5
  ) +
  scale_fill_manual(
    values = c(
      "Case Control" = "#2b8cbe",
      "Proxy" = "#a6bddb"
    )
  ) +
  scale_color_manual(
    values = c(
      "Case Control" = "#2b8cbe",
      "Proxy" = "#a6bddb"
    )
  ) +
  scale_y_continuous(labels = label_comma()) +
  labs(x = "Ancestry", y = "Effective Sample Size", fill = "Type") +
  theme_classic(base_size = 14) +
  theme(legend.title = element_blank()) +
  guides(color = "none", fill = "none") +  # Remove legend for cleaner look
  ylab("")  # Remove y-axis label

ggsave(here("analysis/figures/stacked_bar_noEUR.svg"), plot = p, width = 6, height = 6)
