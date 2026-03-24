riskloc <- function(sumstats, clump_kb = 250, p_th = 5e-8) {

  # clusters genome-wide significant variants into loci by merging
  # windows of +/- clump_kb around each significant variant;
  # returns one row per locus with the lead variant (lowest p-value)
  #
  # args:
  #   sumstats   sumstats df (requires: p_value, base_pair_location,
  #              chromosome, variant_id, effect_allele_frequency)
  #   clump_kb   window size in kb around each significant variant
  #   p_th       significance threshold

  library(here)
  library(dplyr)
  library(data.table)
  library(openxlsx)

  loci <- sumstats %>%
    filter(as.numeric(p_value) < p_th) %>%
    mutate(
      start   = ifelse(base_pair_location <= clump_kb * 1e3, 0, base_pair_location - clump_kb * 1e3),
      end     = base_pair_location + clump_kb * 1e3,
      p_value = as.numeric(p_value),
      maf     = pmin(effect_allele_frequency, 1 - effect_allele_frequency)
    ) %>%
    rename(index_variant_id = variant_id) %>%
    select(chromosome, index_variant_id, start, base_pair_location, maf, end, p_value) %>%
    arrange(chromosome, base_pair_location)

  if (nrow(loci) <= 1) {
    return(loci)
  }

  # assign locus ids by merging overlapping windows;
  # track running max end to correctly handle nested intervals
  loci$locus[1] <- 1
  current_end <- loci$end[1]
  for (i in 2:nrow(loci)) {
    if (loci$start[i] <= current_end & loci$chromosome[i] == loci$chromosome[i - 1]) {
      loci$locus[i] <- loci$locus[i - 1]
      current_end <- max(current_end, loci$end[i])
    } else {
      loci$locus[i] <- loci$locus[i - 1] + 1
      current_end <- loci$end[i]
    }
  }

  # collapse each locus to its lead variant
  loci_lst <- lapply(1:max(loci$locus), function(i) {
    temp     <- loci[loci$locus == i, ]
    lead_idx <- which.min(temp$p_value)

    df <- data.frame(
      chromosome         = temp$chromosome[1],
      index_variant_id   = temp$index_variant_id[lead_idx],
      start              = min(temp$start),
      base_pair_location = temp$base_pair_location[lead_idx],
      maf                = temp$maf[lead_idx],
      end                = max(temp$end),
      p_value            = temp$p_value[lead_idx],
      locus              = temp$locus[1]
    )

    # pin APOE lead to APOE e4 variant (rs429358) when p_value == 0
    if (df$chromosome == 19 & df$start <= 45411941 & df$end >= 45411941 & df$p_value == 0) {
      df$index_variant_id   <- "19:45411941:C:T"
      df$base_pair_location <- 45411941
      df$maf                <- temp$maf[temp$index_variant_id == "19:45411941:C:T"]
    }

    df
  })

  do.call(rbind, loci_lst)
}
