#### 0. setup ####
library(data.table)
library(dplyr)
library(here)
library(ggplot2)
library(pROC)
library(forcats)
source(here("R/r2l_r2o.R"))


#### 1. plot incremental R2 ####
pgs_models <- c(
  "PC[1-10]",
  "APOE",
  "PGS",
  "APOE + PGS",
  "APOE + PGS + sex"
)

pgs_eval <- readRDS("analysis/pgs/pgs_eval.rds")[["sbayesrc"]] %>%
  filter(model != "APOE + PGS + sex") %>%
  group_by(cohort) %>%
  mutate(
    r2l_base = r2l[model == "PC[1-10]"], 
    r2l_delta = r2l - r2l_base    
  ) %>%
  ungroup() 

dodge <- 0.7
p1 <- pgs_eval %>%
  filter(model != "PC[1-10]") %>%
  mutate(
    model = factor(model, levels = pgs_models),
    cohort_label = case_when(
      cohort == "xstsa" ~ "STSA",
      cohort == "demge" ~ "DemGene",
      cohort == "gothe" ~ "Gothenburg",
      cohort == "twing" ~ "TwinGene",
      grepl("xukbb", cohort) ~ "UK Biobank",
    )) %>%
  ggplot(aes(y = r2l_delta, x = reorder(cohort_label, -r2l_delta), fill = model)) +
  stat_summary(
    fun = "mean", 
    geom = "bar", 
    position = position_dodge(width = dodge),
    width = 0.7, 
    colour = "black",
    alpha = 0.7
  ) +
  ylim(0, 0.19) +
  scale_fill_manual(
    values = c("#ece7f2", "#a6bddb", "#2b8cbe"),
    labels = c(
      "APOE" = "APOE e2+e4",
      "PGS" = expression(PGS[full]),
      "APOE + PGS" = bquote("APOE e2+e4 + " * PGS["excl. APOE"])
    ),
    guide = guide_legend(title = "Predictor")
  ) +
  guides(fill=guide_legend(title="Predictor")) +
  xlab("Cohort") +
  ylab(expression("Incremental" ~ italic(R)[liability]^2)) + 
  theme_classic()

ggsave(
  plot = p1,
  filename = paste0("analysis/figures/sbayesrc_r2l_increment.png"),
  device = "png",
  width = 8,
  height = 4
)


#### 2. plot total R2 ####

pgs_eval %>%
  mutate(model = factor(model, levels = pgs_models)) %>%
  ggplot(aes(y = r2l, x = model, fill = cohort)) +
  stat_summary(
    fun = "mean", 
    geom = "bar", 
    position = position_dodge(width = dodge),
    width = 0.7, 
    colour = "black",
    alpha = 0.7
  ) +
  theme_light() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))


#### 3. age [TEMP] ####
pgs_eval <- data.frame(
  cohort = as.character(),
  method = as.character(),
  model = as.character(),
  r2o = as.numeric(),
  r2l = as.numeric(),
  auc = as.numeric(),
  auc_analytic = as.numeric()
)

pgs_models <- c(
  "PC[1-10]",
  "age",
  "APOE",
  "PGS",
  "APOE + PGS",
  "APOE + PGS + sex",
  "APOE + PGS + sex + age"
)

pgs_model <- NULL
pgs_model[[pgs_models[1]]] <- "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
pgs_model[[pgs_models[2]]] <- "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + age"
pgs_model[[pgs_models[3]]] <- "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + APOE2 + APOE4"
pgs_model[[pgs_models[4]]] <- "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + pgs_yeaapoe"
pgs_model[[pgs_models[5]]] <- "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + APOE2 + APOE4 + pgs_nayapoe"
pgs_model[[pgs_models[6]]] <- "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + APOE2 + APOE4 + pgs_nayapoe + sex"
pgs_model[[pgs_models[7]]] <- "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + APOE2 + APOE4 + pgs_nayapoe + sex + age"


cohort <- "demge"

## pgs files
profile_yeaapoe <- fread(here("analysis/pgs/sbayesrc/plink_scores", paste0(cohort, "_eur_yeaapoe.profile"))) %>%
  dplyr::select(FID, IID, SCORESUM) %>%
  rename(pgs_yeaapoe = SCORESUM)
profile_nayapoe <- fread(here("analysis/pgs/sbayesrc/plink_scores", paste0(cohort, "_eur_nayapoe.profile"))) %>%
  dplyr::select(FID, IID, SCORESUM) %>%
  rename(pgs_nayapoe = SCORESUM)

## pheno file
pheno <- fread(here("analysis/pgs/testing_samples", paste0(cohort, ".cov"))) %>%
  filter(Phenotype %in% c(0, 1))

apoe2 <- fread(here("analysis/pgs/testing_samples", paste0(cohort, "_apoe2.raw")))[,c(1,2,7)]
colnames(apoe2)[3] <- "APOE2"

pheno <- pheno %>%
  inner_join(apoe2, by = c("FID", "IID")) %>%
  filter(!is.na(APOE2))

if (cohort == "demge") {pheno <- pheno %>% filter(age >= 40)}

## merge files
m <- pheno %>%
  inner_join(profile_yeaapoe, by = c("FID", "IID")) %>%
  inner_join(profile_nayapoe, by = c("FID", "IID"))

for (i in 1:length(pgs_models)) {
  
  formula_r2 <- as.formula(paste("Phenotype ~", pgs_model[[pgs_models[i]]]))
  r2o <- summary(lm(formula_r2, data = m))$adj.r.sq
  r2l <- prs_r2obs_to_r2liab(
    K = 0.05,
    P = mean(m$Phenotype),
    prs_r2obs = r2o
  )
  mod1 <- glm(formula_r2, data = m, family = "binomial")
  auc <- as.numeric(roc(response = m$Phenotype, predictor = predict(mod1), quiet = TRUE)$auc)
  
  
  ## r20 to auc
  n1    <- sum(m$Phenotype == 1); n2 <- (sum(m$Phenotype == 0))
  alpha <- (n1 + n2)^2 / (n1 * n2)
  d     <- (sqrt(alpha) * sqrt(r2o)) / sqrt(1 - r2o)
  auc_analytic   <- pnorm(d/sqrt(2), 0, 1)
  
  
  pgs_eval[nrow(pgs_eval) + 1, "cohort"] <- cohort
  pgs_eval[nrow(pgs_eval), "method"] <- "sbayesrc"
  pgs_eval[nrow(pgs_eval), "model"] <- pgs_models[i]
  pgs_eval[nrow(pgs_eval), "r2o"] <- r2o
  pgs_eval[nrow(pgs_eval), "r2l"] <- r2l
  pgs_eval[nrow(pgs_eval), "auc"] <- auc 
  pgs_eval[nrow(pgs_eval), "auc_analytic"] <- auc_analytic 
  
}

pgs_eval <- pgs_eval %>%
  group_by(cohort) %>%
  mutate(
    r2l_base = r2l[model == "PC[1-10]"], 
    r2l_delta = r2l - r2l_base    
  ) %>%
  ungroup() 

pgs_eval %>%
  filter(model != "PC[1-10]") %>%
  mutate(model = factor(model, levels = pgs_models)) %>%
  ggplot(aes(y = r2l_delta, x = model)) +
  stat_summary(
    fun = "mean", 
    geom = "bar", 
    position = position_dodge(width = dodge),
    width = 0.7, 
    colour = "black",
    alpha = 0.7
  ) +
  theme_light() +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

## total r2l
pgs_eval %>%
  mutate(model = factor(model, levels = pgs_models)) %>%
  ggplot(aes(y = r2l, x = model)) +
  stat_summary(
    fun = "mean", 
    geom = "bar", 
    position = position_dodge(width = dodge),
    width = 0.7, 
    colour = "black",
    fill = "black",
    alpha = 0.7
  ) +
  theme_light() +
  scale_x_discrete(
    guide = guide_axis(n.dodge = 2), 
    expand = expansion(add = c(.6, 1.2))
  )
                   

#### 4. SBayesRC-mult ####
pgs_models <- c(
  "PC[1-10]",
  "PGS_eur",
  "PGS_eas",
  "PGS_afr",
  "PGS_all"
)

pgs_eval <- readRDS("analysis/pgs/pgs_eval.rds")[["sbayesrc_mult"]] %>%
  group_by(cohort) %>%
  mutate(
    r2l_base = r2l[model == "PC[1-10]"], 
    r2l_delta = r2l - r2l_base    
  ) %>%
  ungroup() %>%
  mutate(ancestry = case_when(
    cohort %in% c("xstsa", "demge", "twing", "xukbb_eur_casec", "gothe") ~ "EUR",
    cohort == "xukbb_afr_casec" ~ "AFR",
    cohort == "BBJ" ~ "EAS",
    cohort == "xukbb_sas_casec" ~ "SAS"
  ))

dodge <- 0.7
p2 <- pgs_eval %>%
  filter(model %in% c("PGS_eur", "PGS_all")) %>%
  mutate(
    model = case_when(
      model == "PGS_eur" ~ "EUR",
      model == "PGS_all" ~ "MULTI"
    ),
    model = factor(model, levels = c("EUR", "MULTI")),
    ancestry = factor(ancestry, levels = c("EUR", "AFR","EAS", "SAS")),
    cohort_label = case_when(
      cohort == "xstsa" ~ "STSA",
      cohort == "demge" ~ "DemGene",
      cohort == "gothe" ~ "Gothenburg",
      cohort == "twing" ~ "TwinGene",
      grepl("xukbb", cohort) ~ "UK Biobank",
      cohort == "BBJ" ~ "BioBank Japan"
    )) %>%
  ggplot(aes(y = r2l_delta, x = fct_reorder(cohort_label, desc(r2l_delta), .fun='min'), fill = model)) +
  stat_summary(
    fun = "mean", 
    geom = "bar", 
    position = position_dodge(width = dodge),
    width = 0.7, 
    colour = "black",
    alpha = 0.7
  ) +
  facet_grid(~ancestry, scales = "free", space = "free") +
  ylim(0, 0.18) +
  scale_fill_manual(values = c("#ece7f2", "#2b8cbe")) +
  guides(fill=guide_legend(title="PGS Model")) +
  xlab("Cohort") +
  ylab(expression("Incremental" ~ italic(R)[liability]^2)) + 
  theme_classic()

ggsave(
  plot = p2,
  filename = paste0("analysis/figures/sbsayesrc_mult_r2l_increment.png"),
  device = "png",
  width = 8,
  height = 4
)
