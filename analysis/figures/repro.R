# =============================================================
#
#  repro.R
#  Reproducibility analysis: compare EADB (Bellenguez 2022)
#  and PGC-ALZ2 (Wightman 2021) risk loci against PGC-ALZ3.
#
#  For each prior study, two quantities are extracted per locus
#  in PGC-ALZ3: the minimum p-value within the locus window,
#  and the p-value / effect size of the exact index variant.
#
#  Table of Contents:
#    00. Set-up
#    01. EADB loci — p-values in PGC-ALZ3
#    02. PGC-ALZ2 loci — p-values in PGC-ALZ3
#    03. EADB loci — effect sizes in PGC-ALZ3
#    04. Plot beta correlations (EADB vs. PGC-ALZ3)
#
# =============================================================

#### 00. Set-up ####

library(data.table)
library(dplyr)
library(here)
library(ggplot2)
source(here("R", "riskloc.R"))

FIG_W <- 6
FIG_H <- 5

# -------------------------------------------------------------
# Helper: for each locus in `riskloci`, extract the minimum
# p-value and the index-variant p-value from `sst`
# -------------------------------------------------------------
extract_locus_pvalues <- function(riskloci, sst) {
  riskloci$pgc3_min_p_value <- NA_real_
  riskloci$pgc3_ind_p_value <- NA_real_

  for (i in seq_len(nrow(riskloci))) {
    locus_snps <- sst %>%
      filter(
        chromosome         == riskloci$chromosome[i],
        base_pair_location >= riskloci$start[i],
        base_pair_location <= riskloci$end[i]
      ) %>%
      mutate(p_value = as.numeric(p_value))

    if (nrow(locus_snps) > 0) {
      riskloci$pgc3_min_p_value[i] <- min(locus_snps$p_value, na.rm = TRUE)
    }

    index_snp <- sst %>%
      filter(variant_id == riskloci$index_variant_id[i]) %>%
      mutate(p_value = as.numeric(p_value))

    if (nrow(index_snp) > 0) {
      riskloci$pgc3_ind_p_value[i] <- index_snp$p_value
    }
  }

  riskloci
}

# -------------------------------------------------------------
# Helper: beta vs. beta scatter with 95% CI error bars and
# Pearson r annotation
# -------------------------------------------------------------
beta_scatter <- function(data, x, y, x_se, y_se, x_label, y_label) {
  r <- cor(data[[x]], data[[y]], method = "pearson", use = "complete.obs")

  ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
    geom_point() +
    geom_errorbar(
      aes(
        ymin = .data[[y]] - 1.96 * .data[[y_se]],
        ymax = .data[[y]] + 1.96 * .data[[y_se]]
      ),
      width = 0
    ) +
    geom_errorbarh(
      aes(
        xmin = .data[[x]] - 1.96 * .data[[x_se]],
        xmax = .data[[x]] + 1.96 * .data[[x_se]]
      ),
      height = 0
    ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "lightgrey") +
    annotate(
      "text", x = 0, y = 0.6,
      label  = bquote(italic(r) == .(round(r, 2))),
      size   = 4, colour = "#084594"
    ) +
    labs(x = x_label, y = y_label) +
    ggtitle(expression("Log-Odds" %+-% " Standard Error")) +
    xlim(-0.2, 0.85) +
    ylim(-0.2, 0.85) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
}


#### 01. EADB loci — p-values in PGC-ALZ3 ####

eadb <- fread(here("data", "sumstats", "combined", "clean_xeadb_combined_noag_eur.txt.gz"))

pgc3_eadb <- list(
  pgc3_noeadb = fread(here("data", "sumstats", "meta", "eadb",
                            "combined", "main_combined_all_neff0.6_nsumstats1.txt.gz")),
  pgc3_full   = fread(here("data", "sumstats", "meta", "main",
                            "combined", "main_combined_all_neff0.6_nsumstats1.txt.gz"))
)

riskloci_eadb <- riskloc(sumstats = eadb) %>%
  mutate(study = "eadb")

riskloci_ps_eadb <- lapply(pgc3_eadb, function(sst) {
  extract_locus_pvalues(riskloci_eadb, sst)
})

saveRDS(riskloci_ps_eadb, file = here("analysis", "repro", "eadb.rds"))


#### 02. PGC-ALZ2 loci — p-values in PGC-ALZ3 ####

pgc2 <- fread(here(
    "data_raw", "published_sumstats", "ad_wightman",
    "PGCALZ2sumstatsExcluding23andMe.txt.gz"
  )) %>%
  rename(
    chromosome          = chr,
    base_pair_location  = PosGRCh37,
    effect_allele       = testedAllele,
    other_allele        = otherAllele,
    p_value             = p
  ) %>%
  mutate(
    variant_id = paste(
      chromosome, base_pair_location,
      pmin(effect_allele, other_allele),
      pmax(effect_allele, other_allele),
      sep = ":"
    ),
    effect_allele_frequency = 1   # placeholder required by riskloc()
  )

pgc3_pgc2 <- list(
  pgc3_nopgc2 = fread(here("data", "sumstats", "meta", "newc",
                            "combined", "newc_combined_all_neff0.6_nsumstats1.txt.gz")),
  pgc3_full   = pgc3_eadb[["pgc3_full"]]
)

riskloci_pgc2 <- riskloc(sumstats = pgc2) %>%
  mutate(study = "pgc2")

riskloci_ps_pgc2 <- lapply(pgc3_pgc2, function(sst) {
  extract_locus_pvalues(riskloci_pgc2, sst)
})

saveRDS(riskloci_ps_pgc2, file = here("analysis", "repro", "pgc2.rds"))


#### 03. EADB loci — effect sizes in PGC-ALZ3 ####

# Extract beta, SE, and effect allele for the index variant of each EADB locus
# from EADB itself and from both PGC-ALZ3 versions, then align to EADB allele coding.

extract_index_stats <- function(sst, riskloci, suffix) {
  cols <- c("beta", "standard_error", "effect_allele")
  out  <- setNames(
    as.data.frame(matrix(NA_real_, nrow = nrow(riskloci), ncol = 3)),
    paste0(cols, "_", suffix)
  )
  out[[paste0("effect_allele_", suffix)]] <- NA_character_

  for (i in seq_len(nrow(riskloci))) {
    row <- sst[sst$variant_id == riskloci$index_variant_id[i], ]
    if (nrow(row) > 0) {
      out[[paste0("beta_",             suffix)]][i] <- as.numeric(row$beta[1])
      out[[paste0("standard_error_",   suffix)]][i] <- as.numeric(row$standard_error[1])
      out[[paste0("effect_allele_",    suffix)]][i] <- row$effect_allele[1]
    }
  }
  out
}

riskloci_eadb_beta <- riskloci_eadb %>%
  bind_cols(extract_index_stats(eadb,                   riskloci_eadb, "eadb")) %>%
  bind_cols(extract_index_stats(pgc3_eadb$pgc3_noeadb,  riskloci_eadb, "pgc3_noeadb")) %>%
  bind_cols(extract_index_stats(pgc3_eadb$pgc3_full,    riskloci_eadb, "pgc3_full")) %>%
  na.omit() %>%
  # Align pgc3 betas to EADB effect allele coding
  mutate(
    beta_pgc3_noeadb = ifelse(
      effect_allele_pgc3_noeadb != effect_allele_eadb, -beta_pgc3_noeadb, beta_pgc3_noeadb
    ),
    beta_pgc3_full = ifelse(
      effect_allele_pgc3_full != effect_allele_eadb, -beta_pgc3_full, beta_pgc3_full
    ),
    effect_allele_pgc3_noeadb = effect_allele_eadb,
    effect_allele_pgc3_full   = effect_allele_eadb
  )

saveRDS(riskloci_eadb_beta, file = here("analysis", "repro", "eadb_beta_se.rds"))


#### 04. Plot beta correlations (EADB vs. PGC-ALZ3) ####

riskloci_eadb_beta <- readRDS(here("analysis", "repro", "eadb_beta_se.rds"))

p1 <- beta_scatter(
  data    = riskloci_eadb_beta,
  x       = "beta_eadb",          x_se = "standard_error_eadb",
  y       = "beta_pgc3_noeadb",   y_se = "standard_error_pgc3_noeadb",
  x_label = "Bellenguez et al. (2022)",
  y_label = "PGC-ALZ3 (excl. overlap)"
)

ggsave(
  here("analysis", "repro", "pgc3noeadb_vs_eadb_beta_se.png"),
  plot = p1, width = FIG_W, height = FIG_H
)

p2 <- beta_scatter(
  data    = riskloci_eadb_beta,
  x       = "beta_eadb",       x_se = "standard_error_eadb",
  y       = "beta_pgc3_full",  y_se = "standard_error_pgc3_full",
  x_label = "Bellenguez et al. (2022)",
  y_label = "PGC-ALZ3"
)

ggsave(
  here("analysis", "repro", "pgc3full_vs_eadb_beta_se.png"),
  plot = p2, width = FIG_W, height = FIG_H
)
