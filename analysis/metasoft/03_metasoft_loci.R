# ------------------------------------------------------------
# identify genome-wide significant loci from metasoft output
# compares loci from fixed-effects (fe) and re2 tests
# ------------------------------------------------------------

library(here)
library(data.table)
library(dplyr)


# -- load data -----------------------------------------------

sst <- fread(
  here("analysis/metasoft/output/combined_out.txt"),
  select = c("RSID", "PVALUE_FE", "PVALUE_RE", "PVALUE_RE2", "I_SQUARE", "Q", "PVALUE_Q")
)


# -- function: define risk loci ------------------------------

# clusters genome-wide significant variants into loci by merging
# windows of +/- clump_kb around each significant variant;
# returns one row per locus with the lead variant (lowest p-value)
riskloc <- function(sumstats, clump_kb = 250, p_th = 5e-8, p_col = "PVALUE_RE2") {

  loci <- sumstats %>%
    mutate(
      chromosome         = as.integer(sub(":.*", "", RSID)),
      base_pair_location = as.integer(sub("^[^:]+:([^:]+):.*", "\\1", RSID)),
      p_value            = as.numeric(.data[[p_col]]),
      start              = ifelse(base_pair_location <= clump_kb * 1e3, 0, base_pair_location - clump_kb * 1e3),
      end                = base_pair_location + clump_kb * 1e3
    ) %>%
    filter(p_value < p_th) %>%
    rename(index_variant_id = RSID) %>%
    select(chromosome, index_variant_id, start, base_pair_location, end, p_value) %>%
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
      end                = max(temp$end),
      p_value            = temp$p_value[lead_idx],
      locus              = temp$locus[1]
    )

    # pin APOE lead to APOE e4 variant (rs429358) when p_value == 0
    if (df$chromosome == 19 & df$start <= 45411941 & df$end >= 45411941 & df$p_value == 0) {
      df$index_variant_id   <- "19:45411941:C:T"
      df$base_pair_location <- 45411941
    }

    df
  })

  do.call(rbind, loci_lst)
}


# -- run for fe and re2 tests --------------------------------

fe_loci  <- riskloc(sst, p_col = "PVALUE_FE") %>% mutate(test = "fe")
re2_loci <- riskloc(sst)                        %>% mutate(test = "re2")


# -- merge overlapping loci across tests ---------------------

# combine and sort so the interval-merging loop below works correctly
temp <- rbindlist(list(fe = fe_loci, re2 = re2_loci)) %>%
  arrange(chromosome, base_pair_location)

# assign shared locus ids to regions that overlap between fe and re2
temp$locus_unique[1] <- 1
current_end <- temp$end[1]
for (i in 2:nrow(temp)) {
  if (temp$start[i] <= current_end & temp$chromosome[i] == temp$chromosome[i - 1]) {
    temp$locus_unique[i] <- temp$locus_unique[i - 1]
    current_end <- max(current_end, temp$end[i])
  } else {
    temp$locus_unique[i] <- temp$locus_unique[i - 1] + 1
    current_end <- temp$end[i]
  }
}

# one row per locus; lead variant is whichever has the lowest p-value
# across both tests; fe/re2 p-values are then looked up from sst for
# that variant, so both are always populated regardless of significance
loci_df <- temp %>%
  group_by(locus_unique) %>%
  summarize(
    chromosome       = first(chromosome),
    index_variant_id = index_variant_id[which.min(p_value)],
    start            = min(start),
    end              = max(end),
    test             = paste(sort(unique(test)), collapse = ";")
  ) %>%
  rename(locus = locus_unique) %>%
  mutate(size_bp = end - start) %>%
  left_join(
    sst %>% select(RSID, p_value_fe = PVALUE_FE, p_value_re2 = PVALUE_RE2),
    by = c("index_variant_id" = "RSID")
  )
