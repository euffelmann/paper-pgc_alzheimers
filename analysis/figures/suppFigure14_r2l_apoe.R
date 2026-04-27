#### 0. set-up ####
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
library(pROC)

prs_r2obs_to_r2liab <- function(K,P,prs_r2obs){
  ## Lee et al. 2012 Genet Epidemiology
  t = -qnorm(K,mean=0,sd=1) # disease threshold
  z<-dnorm(t)               # height of the normal distribution at T
  i1<-z/K                   # mean liability of A1 (eg Falconer and Mackay) 
  k1=i1*(i1-t)              # reduction in variance in A1
  i0<--z/(1-K)              # mean liability of A0
  k0=i0*(i0-t)              # reduction in variance in A0
  
  theta = i1*(P-K)/(1-K)*(i1*(P-K)/(1-K)-t) # theta in equation (15) Lee et al. 2012 Genet Epidemiology
  cv = K*(1-K)/z^2*K*(1-K)/(P*(1-P))        # C in equation (15)
  R2 = prs_r2obs*cv/(1+prs_r2obs*theta*cv)
  
  return(R2)
}

prs_r2liab_to_r2obs <- function(K,P,prs_r2liab){
  ## Lee et al. 2012 Genet Epidemiology
  t = -qnorm(K,mean=0,sd=1) # disease threshold
  z<-dnorm(t)               # height of the normal distribution at T
  i1<-z/K                   # mean liability of A1 (eg Falconer and Mackay) 
  k1=i1*(i1-t)              # reduction in variance in A1
  i0<--z/(1-K)              # mean liability of A0
  k0=i0*(i0-t)              # reduction in variance in A0
  
  theta = i1*(P-K)/(1-K)*(i1*(P-K)/(1-K)-t) # theta in equation (15) Lee et al. 2012 Genet Epidemiology
  cv = K*(1-K)/z^2*K*(1-K)/(P*(1-P))        # C in equation (15)
  
  prs_r2obs = prs_r2liab / {cv-prs_r2liab*theta*cv}
  
  return(prs_r2obs)
}


## download inputs
system("dx download -r Dougs_data:/PGCALZ3/PGCALZ3_PRS/2026_04_10_pgsAnalysis_emil/testing_samples -o .")
system("dx download -r Dougs_data:/PGCALZ3/PGCALZ3_PRS/2026_04_10_pgsAnalysis_emil/plink_scores -o .")
system("dx download -r Dougs_data:/PGCALZ3/PGCALZ3_PRS/2026_04_10_pgsAnalysis_emil/pgs_eval.rds -o .")


#### 01. run regression ####
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
  "PGS_exclAPOE"
)

pgs_model <- NULL
pgs_model[[pgs_models[1]]] <- "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + pgs_nayapoe"


for (cohort in c("xstsa", "gothe", "twing", "demge", "xukbb_eur_casec")) { 
  
  ### load pgs scores 
  ## (yeaapoe = including APOE region; nayapoe = excluding APOE region)
  profile_yeaapoe <- fread(paste0("plink_scores/", cohort, "_eur_yeaapoe.profile")) %>%
    dplyr::select(FID, IID, SCORESUM) %>%
    rename(pgs_yeaapoe = SCORESUM)
  profile_nayapoe <- fread(paste0("plink_scores/", cohort, "_eur_nayapoe.profile")) %>%
    dplyr::select(FID, IID, SCORESUM) %>%
    rename(pgs_nayapoe = SCORESUM)
  
  ## pheno file
  pheno <- fread(paste0("testing_samples/", cohort, ".cov")) %>%
    rename(
      Phenotype = ifelse(grepl("xukbb", cohort), "AD", "Phenotype"),
      PC1 = ifelse(grepl("xukbb", cohort), "pop_pc1", "PC1"),
      PC2 = ifelse(grepl("xukbb", cohort), "pop_pc2", "PC2"),
      PC3 = ifelse(grepl("xukbb", cohort), "pop_pc3", "PC3"),
      PC4 = ifelse(grepl("xukbb", cohort), "pop_pc4", "PC4"),
      PC5 = ifelse(grepl("xukbb", cohort), "pop_pc5", "PC5"),
      PC6 = ifelse(grepl("xukbb", cohort), "pop_pc6", "PC6"),
      PC7 = ifelse(grepl("xukbb", cohort), "pop_pc7", "PC7"),
      PC8 = ifelse(grepl("xukbb", cohort), "pop_pc8", "PC8"),
      PC9 = ifelse(grepl("xukbb", cohort), "pop_pc9", "PC9"),
      PC10 = ifelse(grepl("xukbb", cohort), "pop_pc10", "PC10")) %>%
    filter(Phenotype %in% c(0, 1))
  
  apoe2 <- fread(paste0("testing_samples/", cohort, "_apoe2.raw"))[,c(1,2,7)]
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
    
    
    ## r2o to auc
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
}

pgs_eval_exclAPOE <- pgs_eval
temp <- readRDS("pgs_eval.rds")[["sbayesrc"]] %>%
  filter(method == "sbayesrc_new") 

#### plot
pgs_eval <- rbind(
  temp,
  pgs_eval_exclAPOE
)

pgs_models <- c(
  "PC[1-10]",
  "PGS_exclAPOE",
  "APOE",
  "PGS",
  "APOE + PGS",
  "APOE + PGS + sex"
)

pgs_eval <- pgs_eval %>%
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
    values = c("#f1eef6", "#bdc9e1", "#74a9cf", "#0570b0"),
    labels = c(
      "APOE" = "APOE e2+e4",
      "PGS" = expression(PGS[full]),
      "PGS_exclAPOE" = expression(PGS["excl. APOE"]),
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
  filename = paste0("supFig_r2l_apoe.png"),
  device = "png",
  width = 8,
  height = 4,
  type = "cairo"
)

system("dx upload supFig_r2l_apoe.png --path Dougs_data:/PGCALZ3/PGCALZ3_PRS/2026_04_10_pgsAnalysis_emil/")