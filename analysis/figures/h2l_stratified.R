#### 0. set-up ####
library(dplyr)
library(ggplot2)
library(here)
library(data.table)
library(openxlsx)

#### 1. load data ####

## sbayesrc
sbay_h2 <- read.xlsx(here("analysis/sbayesrc/h2_estimates_14Apr2025.xlsx"), sheet = "14Apr2025") %>%
  filter(Dataset %in% c("fema_case_control_eur_neff0.6_nsumstats1",
                        "male_case_control_eur_neff0.6_nsumstats1",
                        "main_case_control_eur_neff0.6_nsumstats1",
                        "main_combined_eur_neff0.6_nsumstats1",
                        "main_proxy_eur_neff0.6_nsumstats1")) %>%
  mutate(
    method = "sbayesrc",
    sst = case_when(
      grepl("fema", Dataset) ~ "Female",
      grepl("male", Dataset) ~ "Male",
      grepl("main_case_control", Dataset) ~ "Case control",
      grepl("main_combined", Dataset) ~ "Case control + proxy",
      grepl("main_proxy", Dataset) ~ "Proxy"
    ),
    analysis = case_when(
      sst %in% c("Female", "Male") ~ "Sex-stratified (Case control)",
      sst %in% c("Case control", "Case control + proxy", "Proxy") ~ "Phenotype-stratified"
    )) %>%
  rename(h2 = h2.estimate, h2_se = Posterior.SE) %>%
  select(sst, analysis, h2, h2_se, method)

#### 2. make figure ####
dodge <- 0.5
p1 <- sbay_h2 %>%
  ggplot(aes(y = h2, x = sst)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = dodge),
    width = 0.5,
    colour = "black",,
    fill = "#2b8cbe",
    alpha = 0.7
  ) +
  geom_errorbar(
    aes(ymin = h2 - h2_se * 1.96, ymax = h2 + h2_se * 1.96),
    position = position_dodge(width = dodge),
    width = 0.2
  ) +
  facet_grid(~ analysis, scales = "free", space = "free") +
  geom_point(position = position_dodge(width = dodge)) +
  xlab("") +
  ylab(expression(italic(h)[liability]^2)) +
  theme_classic()

ggsave(
  plot = p1,
  filename = paste0("analysis/figures/h2l_stratified.png"),
  device = "png",
  width = 8,
  height = 4
)
