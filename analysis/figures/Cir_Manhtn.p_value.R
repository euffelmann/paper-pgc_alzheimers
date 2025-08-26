#### 0. set-up ####
library(here)
library(data.table)
library(dplyr)
library(CMplot)

temp <- fread(here("data/sumstats/meta/main/combined/main_combined_all_neff0.6_nsumstats1.txt.gz"))
temp2 <- temp %>%
  mutate(
    p_value = as.numeric(p_value),
    p_value = ifelse(p_value == 0, 1e-320, p_value)) %>%
  filter(p_value < 5e-5) %>%
  select(variant_id, chromosome, base_pair_location, p_value) %>%
  rename(
    SNP = variant_id,
    pos = base_pair_location
  )


setwd(here("analysis/figures"))

CMplot(
  as.data.frame(temp2),
  plot.type = "c",
  r = 0.9,
  outward = TRUE,
  cir.chr.h = .1 ,
  chr.den.col = "orange",
  file = "pdf",
  col = c("#3690c0", "#bdbdbd"),
  ylim = c(0, 120),
  cex = 0.3,
  dpi = 300,
  chr.labels = seq(1, 22)
)


temp3 <- temp2 %>%
  rename(p_value_75 = p_value)

CMplot(
  as.data.frame(temp3),
  plot.type = "c",
  r = 0.9,
  outward = TRUE,
  cir.chr.h = .1 ,
  chr.den.col = "orange",
  file = "pdf",
  col = c("#3690c0", "#fc8d59"),
  ylim = c(0, 75),
  cex = 0.3,
  dpi = 300,
  chr.labels = seq(1, 22)
)
