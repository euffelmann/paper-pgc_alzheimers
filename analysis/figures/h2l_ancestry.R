#### 0. set-up ####
library(dplyr)
library(ggplot2)
library(here)
library(data.table)
library(openxlsx)

#### 1. load data ####
## ldsc
ldsc_h2 <- fread(here("analysis/ldsc/ldsc_h2.txt")) %>%
  filter(
    proxy_casec == "case_control",
    analysis == "main",
    anc %in% c("eur", "afr", "eas")) %>%
  mutate(method = "ldsc") %>%
  select(anc, h2, h2_se, method)

## sbayesrc
sbay_h2 <- read.xlsx(here("analysis/sbayesrc/h2_estimates_14Apr2025.xlsx"), sheet = "14Apr2025") %>%
  filter(Dataset %in% c("main_case_control_eur_neff0.6_nsumstats1",
                        "main_case_control_eas_neff0.6_nsumstats1",
                        "main_case_control_afr_neff0.6_nsumstats1")) %>%
  mutate(
    anc = tolower(Ancestry),
    method = "sbayesrc") %>%
  rename(h2 = h2.estimate, h2_se = Posterior.SE) %>%
  select(anc, h2, h2_se, method)


h2 <- rbind(ldsc_h2, sbay_h2)


#### 2. make figure ####
dodge <- 0.5
p1 <- h2 %>%
  mutate(anc = toupper(anc)) %>%
  ggplot(aes(y = h2, x = reorder(anc, -h2), fill = method)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = dodge),
    width = 0.5,
    colour = "black",
    alpha = 0.7
  ) +
  geom_errorbar(
    aes(ymin = h2 - h2_se * 1.96, ymax = h2 + h2_se * 1.96),
    position = position_dodge(width = dodge),
    width = 0.2
  ) +
  geom_point(position = position_dodge(width = dodge)) +
  scale_fill_manual(
    values = c("#d73027", "#2b8cbe"),
    labels = c(
      "ldsc" = "LDSC",
      "sbayesrc" = "SBayesRC"
    ),
    guide = guide_legend(title = "Method")
  ) +
  xlab("Ancestry") +
  ylab(expression(italic(h)[liability]^2)) +
  theme_classic()

ggsave(
  plot = p1,
  filename = paste0("analysis/figures/h2l_ancestry.png"),
  device = "png",
  width = 8,
  height = 4
)
