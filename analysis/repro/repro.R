#### 0. setup ####
library(data.table)
library(dplyr)
library(here)
library(data.table)
library(ggplot2)
source(here("R/riskloc.R"))


#### 1. risk loci P (eadb vs. pgc-alz3)  ####

eadb <- fread(here("data/sumstats/combined/clean_xeadb_combined_noag_eur.txt.gz"))
pgc3 <- NULL
pgc3[["pgc3_noeadb"]] <- fread(here("data/sumstats/meta/eadb/combined/main_combined_all_neff0.6_nsumstats1.txt.gz"))
pgc3[["pgc3_full"]]   <- fread(here("data/sumstats/meta/main/combined/main_combined_all_neff0.6_nsumstats1.txt.gz"))

## define risk loci
riskloci_eadb <- riskloc(sumstats = eadb) %>%
  mutate(study = "eadb")

## initiate empty columns
riskloci_eadb$pgc3_min_p_value <- NA
riskloci_eadb$pgc3_ind_p_value <- NA

## iniate empty list to save dataframes
riskloci_ps <- list()

for (sst_set in c("pgc3_noeadb", "pgc3_full")) {
  
  for (i in 1:(nrow(riskloci_eadb))) {
    
    index <- riskloci_eadb$index_variant_id[i]
    start <- riskloci_eadb$start[i]
    end <- riskloci_eadb$end[i]
    chr <- riskloci_eadb$chromosome[i]
    
    ## smallest p_value in the same locus
    p_min <- pgc3[[sst_set]] %>%
      filter(base_pair_location >= start,
             base_pair_location <= end,
             chromosome == chr) %>%
      mutate(p_value = as.numeric(p_value)) %>%
      pull(p_value)
    
    ## p_value of the same variant
    if (index %in% pgc3[[sst_set]]$variant_id) {
      p_index <- pgc3[[sst_set]] %>%
        filter(variant_id == index) %>%
        mutate(p_value = as.numeric(p_value)) %>%
        pull(p_value)
    } else {
      p_index <- NA
    }
    
    riskloci_eadb$pgc3_min_p_value[i] <- min(p_min)
    riskloci_eadb$pgc3_ind_p_value[i] <- p_index
    
  }
  
  riskloci_ps[[sst_set]] <- riskloci_eadb
  
}

saveRDS(riskloci_ps, file = here("analysis/repro/eadb.rds"))


#### 2. risk loci P (pgc-alz2 vs. pgc-alz3)  ####

pgc2 <- fread(here("data_raw/published_sumstats/ad_wightman/PGCALZ2sumstatsExcluding23andMe.txt.gz")) %>%
  rename(
    chromosome = chr,
    base_pair_location = PosGRCh37,
    effect_allele = testedAllele,
    other_allele = otherAllele,
    p_value = p
  ) %>%
  mutate(
    variant_id = paste(
      chromosome,
      base_pair_location,
      pmin(effect_allele, other_allele),
      pmax(effect_allele, other_allele),
      sep = ":"
    ),
    effect_allele_frequency = 1 # placeholder
  ) 
pgc3[["pgc3_nopgc2"]]  <- fread(here("data/sumstats/meta/newc/combined/newc_combined_all_neff0.6_nsumstats1.txt.gz"))

## define risk loci
riskloci_pgc2 <- riskloc(sumstats = pgc2) %>%
  mutate(study = "pgc2")

## initiate empty columns
riskloci_pgc2$pgc3_min_p_value <- NA
riskloci_pgc2$pgc3_ind_p_value <- NA

## iniate empty list to save dataframes
riskloci_ps <- list()

for (sst_set in c("pgc3_nopgc2", "pgc3_full")) {
  
  for (i in 1:(nrow(riskloci_pgc2))) {
    
    index <- riskloci_pgc2$index_variant_id[i]
    start <- riskloci_pgc2$start[i]
    end <- riskloci_pgc2$end[i]
    chr <- riskloci_pgc2$chromosome[i]
    
    ## smallest p_value in the same locus
    p_min <- pgc3[[sst_set]] %>%
      filter(base_pair_location >= start,
             base_pair_location <= end,
             chromosome == chr) %>%
      mutate(p_value = as.numeric(p_value)) %>%
      pull(p_value)
    
    ## p_value of the same variant
    if (index %in% pgc3[[sst_set]]$variant_id) {
      p_index <- pgc3[[sst_set]] %>%
        filter(variant_id == index) %>%
        mutate(p_value = as.numeric(p_value)) %>%
        pull(p_value)
    } else {
      p_index <- NA
    }
    
    riskloci_pgc2$pgc3_min_p_value[i] <- min(p_min)
    riskloci_pgc2$pgc3_ind_p_value[i] <- p_index
    
  } 
  
  riskloci_ps[[sst_set]] <- riskloci_pgc2
  
}

saveRDS(riskloci_ps, file = here("analysis/repro/pgc2.rds"))


#### 3. risk loci P (all published loci vs. pgc-alz3) ####
published_loci <- fread(here("analysis/risk_loci/published_loci/ADconsensusLoci.txt")) %>%
  filter(
    Studies != "MegaMeta",
    Chromosome != 23
  ) %>%
  dplyr::rename(
    chromosome = Chromosome,
    start = Start_GRCh37,
    end = End_GRCh37,
    gene = Names
  ) %>%
  dplyr::select(
    chromosome, start, end, gene, Studies
  ) %>%
  mutate(source = "published") %>%
  arrange(chromosome, start)

pgc3 <- fread(here("data/sumstats/meta/main/combined/main_combined_all_neff0.6_nsumstats1.txt.gz"))

published_loci$pgc3_min_p_value <- NA
for (i in 1:nrow(published_loci)) {
  
  print(paste0("locus ", i))
  
  temp <- pgc3 %>%
    filter(
      base_pair_location >= published_loci$start[i],
      base_pair_location <= published_loci$end[i],
      chromosome == published_loci$chr[i]
    )
  
  if (nrow(temp) > 0) {
    p_min <- temp %>%
      mutate(p_value = as.numeric(p_value)) %>%
      pull(p_value)
    
    published_loci$pgc3_min_p_value[i] <- min(p_min, na.rm = T)
  }
  
}

## replicated loci from previous GWAS
for (gwas in c("Wightman 2021", "Bellenguez 2022")) {
  
  temp <- published_loci %>%
    filter(grepl(gwas, Studies))
  sum(temp$pgc3_min_p_value < 5e-8) / nrow(temp)
  temp %>% 
    filter(
      pgc3_min_p_value > 5e-8
    )
  temp %>% 
    filter(
      pgc3_min_p_value > 0.05 / nrow(temp)
    )
  
}


#### 4. risk loci BETA + SE (eadb vs. pgc-alz3)  ####

eadb <- fread(here("data/sumstats/combined/clean_xeadb_combined_noag_eur.txt.gz"))
pgc3 <- NULL
pgc3[["pgc3_noeadb"]] <- fread(here("data/sumstats/meta/eadb/combined/main_combined_all_neff0.6_nsumstats1.txt.gz"))
pgc3[["pgc3_full"]]   <- fread(here("data/sumstats/meta/main/combined/main_combined_all_neff0.6_nsumstats1.txt.gz"))

## define risk loci
riskloci_eadb <- riskloc(sumstats = eadb) %>%
  mutate(study = "eadb")

riskloci_eadb$beta_eadb <- NA
riskloci_eadb$standard_error_eadb <- NA
riskloci_eadb$effect_allele_eadb <- NA

riskloci_eadb$beta_pgc3_noeadb <- NA
riskloci_eadb$standard_error_pgc3_noeadb <- NA
riskloci_eadb$effect_allele_pgc3_noeadb <- NA

riskloci_eadb$beta_pgc3_full <- NA
riskloci_eadb$standard_error_pgc3_full <- NA
riskloci_eadb$effect_allele_pgc3_full <- NA

for (i in 1:nrow(riskloci_eadb)) {
  
  print(i)
  
  snp_id <- riskloci_eadb$index_variant_id[i]
  
  riskloci_eadb$beta_eadb[i] <- eadb$beta[eadb$variant_id == snp_id]
  riskloci_eadb$standard_error_eadb[i] <- eadb$standard_error[eadb$variant_id == snp_id]
  riskloci_eadb$effect_allele_eadb[i] <- eadb$effect_allele[eadb$variant_id == snp_id]
  
  ## extract betas and standard_errors
  if (any(pgc3[["pgc3_noeadb"]]$variant_id == snp_id)) {
    
    riskloci_eadb$beta_pgc3_noeadb[i] <- pgc3[["pgc3_noeadb"]]$beta[pgc3[["pgc3_noeadb"]]$variant_id ==snp_id]
    riskloci_eadb$standard_error_pgc3_noeadb[i] <- pgc3[["pgc3_noeadb"]]$standard_error[pgc3[["pgc3_noeadb"]]$variant_id == snp_id]
    riskloci_eadb$effect_allele_pgc3_noeadb[i] <- pgc3[["pgc3_noeadb"]]$effect_allele[pgc3[["pgc3_noeadb"]]$variant_id == snp_id] 
    
  }
  
  ## extract betas and standard_errors
  if (any(pgc3[["pgc3_full"]]$variant_id == snp_id)) {
    
    riskloci_eadb$beta_pgc3_full[i] <- pgc3[["pgc3_full"]]$beta[pgc3[["pgc3_full"]]$variant_id == snp_id]
    riskloci_eadb$standard_error_pgc3_full[i] <- pgc3[["pgc3_full"]]$standard_error[pgc3[["pgc3_full"]]$variant_id == snp_id]
    riskloci_eadb$effect_allele_pgc3_full[i] <- pgc3[["pgc3_full"]]$effect_allele[pgc3[["pgc3_full"]]$variant_id == snp_id] 
    
  }
}

riskloci_eadb <- riskloci_eadb %>%
  na.omit() %>%
  mutate(
    beta_pgc3_noeadb = ifelse(effect_allele_pgc3_noeadb != effect_allele_eadb, beta_pgc3_noeadb * -1, beta_pgc3_noeadb),
    beta_pgc3_full = ifelse(effect_allele_pgc3_full != effect_allele_eadb, beta_pgc3_full * -1, beta_pgc3_full),
    effect_allele_pgc3_noeadb = ifelse(effect_allele_pgc3_noeadb != effect_allele_eadb, effect_allele_eadb, effect_allele_pgc3_noeadb),
    effect_allele_pgc3_full = ifelse(effect_allele_pgc3_full != effect_allele_eadb, effect_allele_eadb, effect_allele_pgc3_full)
    )

any(riskloci_eadb$effect_allele_eadb != riskloci_eadb$effect_allele_pgc3_noeadb)
any(riskloci_eadb$effect_allele_eadb != riskloci_eadb$effect_allele_pgc3_full)

saveRDS(riskloci_eadb, file = here("analysis/repro/eadb_beta_se.rds"))

riskloci_eadb <- readRDS(here("analysis/repro/eadb_beta_se.rds"))

r_eadb_pgc3noeadb <- cor(
  riskloci_eadb$beta_eadb,
  riskloci_eadb$beta_pgc3_noeadb,
  method = "pearson",
  use = "complete.obs"
)

p1 <- ggplot(riskloci_eadb, aes(x = beta_eadb, y = beta_pgc3_noeadb)) +
  geom_point() +
  geom_errorbar(
    aes(
      ymin = beta_pgc3_noeadb - 1.96 * standard_error_pgc3_noeadb,
      ymax = beta_pgc3_noeadb + 1.96 * standard_error_pgc3_noeadb
    ),
    width = 0
  ) +
  geom_errorbar(
    aes(
      xmin = beta_eadb - 1.96 * standard_error_eadb,
      xmax = beta_eadb + 1.96 * standard_error_eadb
    ),
    width = 0
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "lightgrey") +
  annotate("text", x = 0, y = 0.6, 
           label = bquote(italic(r) == .(round(r_eadb_pgc3noeadb, 2))), 
           size = 4, colour = "#084594") +
  labs(x = "Bellenguez et al. (2022)", y = "PGC-ALZ3 (excl. overlap)") +
  ggtitle(expression("Log-Odds" %+-% " Standard Error")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(-0.2, 0.85) +
  ylim(-0.2, 0.85)

ggsave(plot = p1, filename = here("analysis/repro/pgc2noeadb_vs_eadb_beta_se.png"), width = 6, height = 5)

r_eadb_pgc3full <- cor(
  riskloci_eadb$beta_eadb,
  riskloci_eadb$beta_pgc3_full,
  method = "pearson",
  use = "complete.obs"
)

p2 <- ggplot(riskloci_eadb, aes(x = beta_eadb, y = beta_pgc3_full)) +
  geom_point() +
  geom_errorbar(
    aes(
      ymin = beta_pgc3_full - 1.96 * standard_error_pgc3_full,
      ymax = beta_pgc3_full + 1.96 * standard_error_pgc3_full
    ),
    width = 0
  ) +
  geom_errorbar(
    aes(
      xmin = beta_eadb - 1.96 * standard_error_eadb,
      xmax = beta_eadb + 1.96 * standard_error_eadb
    ),
    width = 0
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "lightgrey") +
  annotate("text", x = 0, y = 0.6, 
           label = bquote(italic(r) == .(round(r_eadb_pgc3full, 2))), 
           size = 4, colour = "#084594") +
  labs(x = "Bellenguez et al. (2022)", y = "PGC-ALZ3") +
  ggtitle(expression("Log-Odds" %+-% " Standard Error")) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlim(-0.2, 0.85) +
  ylim(-0.2, 0.85)

ggsave(plot = p2, filename = here("analysis/repro/pgc2full_vs_eadb_beta_se.png"), width = 6, height = 5)
