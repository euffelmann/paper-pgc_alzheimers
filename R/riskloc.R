riskloc <- function(sumstats, clump_kb = 250, p_th = 5e-8)  {
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # #             function to define risk loci                              # #
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 
  ## - sumstats:    sumstats df for which to define risk loci
  ## - clump_kb:    size of windows around index SNPs, i.e., index SNP +- clump_kb

  #### set-up ####
  library(here)
  library(dplyr)
  library(data.table)
  library(openxlsx)
    
  loci <- sumstats %>%
    filter(as.numeric(p_value) < p_th) %>%
    mutate(
      start = ifelse(base_pair_location <= clump_kb * 1e3, 0, base_pair_location - clump_kb * 1e3),
      end = base_pair_location + clump_kb * 1e3,
      p_value = as.numeric(p_value),
      maf = pmin(effect_allele_frequency, 1 - effect_allele_frequency)
    ) %>%
    rename(index_variant_id = variant_id) %>%
    select(chromosome, index_variant_id, start, base_pair_location, maf, end, p_value) %>%
    arrange(chromosome, base_pair_location)
  
  if (nrow(loci) <= 1) {
    return(loci)
  } else {
    ## define loci to be merged
    loci$locus[1] <- 1
    for (i in 2:nrow(loci)) {
      if (loci$start[i] <= loci$end[i - 1] & loci$chromosome[i] == loci$chromosome[i - 1]) {
        loci$locus[i] <- loci$locus[i - 1]
      } else {
        loci$locus[i] <- loci$locus[i - 1] + 1
      }
    }
    
    ## merge loci
    loci_lst <- NULL
    for (i in 1:max(loci$locus)) {
      
      temp <- loci[loci$locus == i, ]
      
      df <- data.frame(
        chromosome = temp$chromosome[1],
        index_variant_id = temp$index_variant_id[temp$p_value == min(temp$p_value)][1],
        start = min(temp$start),
        base_pair_location = temp$base_pair_location[temp$p_value == min(temp$p_value)][1],
        maf = temp$maf[temp$p_value == min(temp$p_value)][1],
        end = max(temp$end),
        p_value = temp$p_value[temp$p_value == min(temp$p_value)][1],
        locus = temp$locus[1]
      )
      
      ## redefine APOE locus if p_values = 0
      if (any(df$chromosome == 19) & min(df$start) <= 45411941 & max(df$end) >= 45411941 & any(df$p_value == 0)) {
        df <- data.frame(
          chromosome = 19,
          index_variant_id = "19:45411941:C:T",
          start = min(temp$start),
          base_pair_location = 45411941,
          maf = temp$maf[temp$index_variant_id == "19:45411941:C:T"],
          end = max(temp$end),
          p_value = 0,
          locus = temp$locus[1]
        )
      }
      
      loci_lst[[i]] <- df
      
    }
    
    loci_df <- do.call(rbind, loci_lst)
    
    return(loci_df)
  }
} 

