flip <- function(sumstats, chain) {
  ### liftOver can cause strand issues, because it only updates chromosome
  ### positions without taking into account that strand information can change 
  ### between genome builds. flip() flips alleles in sumstats based on strand 
  ### information in a relevant chain file.
  
  ## - sumstats: dataframe of summary statistics with updated chromosome and 
  ##             base pair location after using liftOver. assumes that the 
  ##             file has a column labeled variant_id_hgxx, containing the
  ##             chromosome and base pair location pre-liftOver (i.e. in this
  ##             case liftOver was performed from hg38 to hg19/hg37)
  ## - chain:    directory of file containing the mappings of chromosome and 
  ##             base pair location, as well as strand between 2 genome builds.
  ##             should be created with CrossMap software to make the ucsc and
  ##             emblem chain files human readable.
  ##             e.g. CrossMap viewchain hg38ToHg19.over.chain.gz > hg38ToHg19_ucsc.chain
  ##             see https://crossmap.readthedocs.io/en/latest/#how-crossmap-works
  
  ## library to join sumstats and chain by region
  if (!require(fuzzyjoin)) {
    
    stop("fuzzyjoin not installed")
    
  } else {
    
    ## load and adjust chain file and extract relevant columns
    chain_df <- fread(chain) %>%
      # filter for regions where the strand has changed
      filter(V4 != V8) %>%
      # there are chromosomes with suffixes don't seem to be relevant,
      # because removing them still corrects all strand flips. Hence, I filter them
      # out. More info: https://genome.ucsc.edu/FAQ/FAQdownloads.html#downloadFix
      filter(V5 %in% paste0("chr", c(seq(1,22), "X")),
             V1 %in% paste0("chr", c(seq(1,22), "X"))) %>%
      mutate(chr_hg19 = ifelse(V5 == "chrX", "chr23", V5),
             chr_hg38 = ifelse(V1 == "chrX", "chr23", V1)) %>%
      mutate(chr_hg19 = as.integer(str_extract(chr_hg19, "(?<=chr)\\d+")),
             chr_hg38 = as.integer(str_extract(chr_hg38, "(?<=chr)\\d+"))) %>%
      rename(
        base_pair_location_start_hg38 = V2,
        base_pair_location_end_hg38 = V3,
        strand_hg38 = V4,
        base_pair_location_start_hg19 = V6,
        base_pair_location_end_hg19 = V7,
        strand_hg19 = V8
      ) %>%
      select(
        chr_hg38,
        base_pair_location_start_hg38,
        base_pair_location_end_hg38,
        strand_hg38,
        chr_hg19,
        base_pair_location_start_hg19,
        base_pair_location_end_hg19,
        strand_hg19
      )
    
    ## extract previous chromosome and base pair location
    sumstats <- sumstats %>%
      mutate(
        chromosome_hg38 = as.integer(str_extract(variant_id_hg38, "^\\d+")),
        base_pair_location_hg38 = as.integer(str_match(variant_id_hg38, "^\\d+:(\\d+):")[, 2]) # [,2] to select the capture group
      )
    
    ## perform a range-based join for new hg19 positions
    sumstats_with_chain_hg19 <- genome_join(
      sumstats[!is.na(sumstats$base_pair_location),],
      chain_df,
      by = c(
        "chromosome" = "chr_hg19",
        "base_pair_location" = "base_pair_location_start_hg19",
        "base_pair_location" = "base_pair_location_end_hg19"
      ),
      mode = "inner"
    )
    
    ## perform a range-based join for old hg38 positions
    sumstats_with_chain_hg38 <- genome_join(
      sumstats[!is.na(sumstats$base_pair_location_hg38),],
      chain_df,
      by = c(
        "chromosome_hg38" = "chr_hg38",
        "base_pair_location_hg38" = "base_pair_location_start_hg38",
        "base_pair_location_hg38" = "base_pair_location_end_hg38"
      ),
      mode = "inner"
    )
    
    ## join to find SNPs that have matched both hg19 and hg38 regions
    # this is necessary because regions in hg19 in the chain file can be
    # overlapping, such that a SNP may fall within a region whose strand is flipped
    # but would actually be mapped to a region that hasn't. Hence, I need to find
    # the rows in the chain file where the SNP's hg38 base pair location falls inside
    # the chain's hg38 region, and the SNP's hg19 base pair locations that falls inside
    # the chains hg19 region. 
    sumstats_with_chain <- sumstats_with_chain_hg19 %>%
      inner_join(
        sumstats_with_chain_hg38,
        by = c(
          "variant_id",
          "chr_hg38",
          "base_pair_location_start_hg38",
          "base_pair_location_end_hg38",
          "chr_hg19",
          "base_pair_location_start_hg19",
          "base_pair_location_end_hg19"
        )
      )
    
    ## flip alleles
    sumstats_flip <- sumstats %>%
      mutate(flip = ifelse(variant_id %in% sumstats_with_chain$variant_id, 1, 0)) %>%
      mutate(
        effect_allele = case_when(
          flip == 1 & effect_allele == "C" ~ "G",
          flip == 1 & effect_allele == "G" ~ "C",
          flip == 1 & effect_allele == "A" ~ "T",
          flip == 1 & effect_allele == "T" ~ "A",
          flip == 0 ~ effect_allele),
        other_allele = case_when(
          flip == 1 & other_allele == "C" ~ "G",
          flip == 1 & other_allele == "G" ~ "C",
          flip == 1 & other_allele == "A" ~ "T",
          flip == 1 & other_allele == "T" ~ "A",
          flip == 0 ~ other_allele))
    
    return(sumstats_flip)
    
  }
}

