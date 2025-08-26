process_sumstats <- function(sumstats_file,
                             outfile,
                             ancestry, 
                             proxy_casecontrol, 
                             outdir,
                             dosage_or_hardcalls,
                             maf_filter = 0.01) {
    
  ## function to process/clean sumstats
  ## input:
  ##  sumstats_file:        sumstats file in GWAS catalog format
  ##  outfile:              name of file to save
  ##  ancestry:             ancestry of GWAS population
  ##  proxy_casecontrol:    sumstats of proxy or case-control phenotype
  ##  outdir:               output directory
  ##  dosage_or_hardcalls:  GWAS based on dosage or hard calls

  ## only run if cleaned sumstats do not exist
  if (!file.exists(here(outdir, paste0("clean_", outfile, ".txt.gz")))) {
    
    print(paste0("processing: ", outfile))
    
    ## load qc_summary
    cohort <- stringr::str_extract(outfile, ".*(?=_case_control|_proxy|_combined)")
    qc_summary <- fread(here("data", "sumstats", "qc_summary", paste0("qc_summary_", cohort, ".txt")), data.table = F)
    
    ## add file name to qc_summary if it's not already there
    if (!outfile %in% qc_summary$file) {
      qc_summary[nrow(qc_summary) + 1, "file"] <- outfile 
    }
    index <- which(qc_summary$file == outfile)
    
    ## load sumstats
    sumstats <- fread(sumstats_file, na.strings = c("", "NA"))
    
    ## basic qc
    if (sumstats$chromosome[1] == "chr1") {
      sumstats$chromosome <- as.integer(sub("chr", "", sumstats$chromosome))
    }
    if ("X" %in% sumstats$chromosome) {
      sumstats$chromosome[sumstats$chromosome == "X"] <- 23
      sumstats$chromosome <- as.integer(sumstats$chromosome)
    }
    sumstats$p_value <- as.numeric(sumstats$p_value)
    
    ## flip alleles for regions with strand changes between hg38 and hg19 if liftOver was used
    # summary statistics are always reported on the + strand (imputation forces this)
    # after using liftOver, some variants will be on the - strand, hence the flipping
    # NOTE: there can be warning messages "Detected an unexpected many-to-many relationship between `x` and `y`."
    # this is due to duplicated in the sumstats that will be removed later in this script
    if ("variant_id_hg38" %in% colnames(sumstats)) {
      print("flipping alleles due to strand conflicts ...")
      sumstats <- flip(sumstats, here("data", "reference_data", "chains", "hg38ToHg19_ucsc.chain"))
      print(paste0("flipped ", sum(sumstats$flip), " alleles"))
    }
    
    ## load reference file
    if (ancestry == "eur") {
      ref_file <- fread(here( "data", "reference_data", "hrc_eur", "hrc_eur.frq")) %>%
        mutate(variant_id = paste0(CHR, ":", POS, ":", pmin(A1, A2), ":", pmax(A1, A2)))
    }
    if (ancestry == "eas") {
      ref_file <- fread(here( "data", "reference_data", "jpb_eas", "jpb_eas.frq")) %>%
        mutate(variant_id = paste0(CHR, ":", POS, ":", pmin(A1, A2), ":", pmax(A1, A2)))
    }
    if (ancestry == "afr") {
      ref_file <- fread(here( "data", "reference_data", "agr_afr", "agr_afr.frq")) %>%
        mutate(variant_id = paste0(CHR, ":", POS, ":", pmin(A1, A2), ":", pmax(A1, A2)))
    }
    if (ancestry == "amr") {
      ref_file <- fread(here( "data", "reference_data", "g1000_amr", "g1000_amr.frq")) %>%
        mutate(variant_id = paste0(CHR, ":", POS, ":", pmin(A1, A2), ":", pmax(A1, A2)))
    }
    if (ancestry == "sas") {
      ref_file <- fread(here( "data", "reference_data", "g1000_sas", "g1000_sas.frq")) %>%
        mutate(variant_id = paste0(CHR, ":", POS, ":", pmin(A1, A2), ":", pmax(A1, A2)))
    }
    
    ## clean sumstats
    temp <- sumstats %>%
      {nrow(.)->> qc_summary[index, "initial_number"];.} %>%
      
      ## change variant_id to alphabetical order of alleles
      # --> important for metal that the variant_ids are the same between
      #     input sumstats. The other_allele and effect_allele may differ
      #     between sumstats, leading to different IDs although they're the 
      #     same SNPs. Sorting it alphabetically fixes this. Also there may
      #     be tri-allelic SNPs, where chromosome and base_pair_location may
      #     be the same in two datasets, but the minor allele is different.
      mutate(variant_id = paste0(chromosome, ":", base_pair_location, ":",
                                 pmin(effect_allele, other_allele), ":",
                                 pmax(effect_allele, other_allele))) %>%
      
      ## remove missing values
      {sum(
        is.na(.$chromosome) |
        is.na(.$base_pair_location) |
        is.na(.$effect_allele) |
        is.na(.$other_allele) |
        is.na(.$beta) |
        is.na(.$standard_error) |
        is.na(.$p_value) |
        is.na(.$effect_allele_frequency)
      ) ->> qc_summary[index, "miss_filter"];.} %>%
      filter(
        !is.na(chromosome), 
        !is.na(base_pair_location),
        !is.na(effect_allele),
        !is.na(other_allele),
        !is.na(beta),
        !is.na(standard_error),
        !is.na(p_value),
        !is.na(effect_allele_frequency)
      ) %>%
      
      ## remove indels 
      {sum(!(
        (.$effect_allele %in% c("A", "C", "T", "G")) & 
        (.$other_allele %in% c("A", "C", "T", "G")))
        ) ->> qc_summary[index, "indel_filter"];.} %>%
      filter(
        (effect_allele %in% c("A", "C", "T", "G")) & 
        (other_allele %in% c("A", "C", "T", "G"))
      ) %>%
      
      ## remove duplicates (keeping only one occurrence)
      {sum(duplicated(.$variant_id)) ->> qc_summary[index, "dupli_filter"];.} %>%
      {print(paste0("removed ", sum(duplicated(.$variant_id)), " duplicated variants"));.} %>%
      filter(duplicated(variant_id) == FALSE) %>%
      
      ## remove multiallelic SNPs
      # in GWAS data such SNPs can flag genotyping errors
      mutate(cptid = paste(chromosome, base_pair_location, sep = ":")) %>%
      {sum(.$cptid %in% .$cptid[duplicated(.$cptid)]) ->> qc_summary[index, "multi_filter"];.} %>%
      {print(paste0("detected ", sum(.$cptid %in% .$cptid[duplicated(.$cptid)]), " multi-allelic variants"));.} %>%
      filter(!(cptid %in% cptid[duplicated(.$cptid)])) %>%
      
      ## remove low maf
      rowwise() %>%
      mutate(maf = min(effect_allele_frequency, 1 - effect_allele_frequency)) %>%
      ungroup() %>%
      {sum(.$maf <= maf_filter) ->> qc_summary[index, "maf_filter"];.} %>%
      {print(paste0("removed ", sum(.$maf <= maf_filter), " variants with maf <= ", maf_filter));.} %>%
      filter(maf >= maf_filter) %>%
      
      ## remove monomorphic SNPs
      {sum(.$effect_allele == .$other_allele) ->> qc_summary[index, "mono_filter"];.} %>%
      {print(paste0("removed ", sum(.$effect_allele == .$other_allele), " monomorphic snps"));.} %>%
      filter(effect_allele != other_allele) %>%
      
      ## remove nonsense values
      {sum(
        .$p_value > 1 | 
        .$p_value < 0 | 
        .$standard_error <= 0 | 
        .$standard_error == Inf | 
        .$beta == Inf
        ) ->> qc_summary[index, "nonsense_filter"];.} %>%
      {print(paste0(
        "removed ", 
        sum(.$p_value > 1 | 
            .$p_value < 0 | 
            .$standard_error <= 0 | 
            .$standard_error == Inf | 
            .$beta == Inf),
        " snps with nonsense values"));.} %>%
      filter(
        p_value < 1,
        p_value >= 0,
        standard_error > 0,
        standard_error != Inf,
        beta != Inf
      ) %>%
      
      ## removed failed liftOver SNPs
      {sum(is.na(.$base_pair_location)) ->> qc_summary[index, "lift_filter"];.} %>%
      {print(paste0("removed ", sum(is.na(.$base_pair_location)), " snps that failed liftOver"));.} %>%
      filter(!is.na(base_pair_location)) %>%
      
      ## remove ambiguous SNPs with maf >= 0.4
      {sum(
        ((.$effect_allele == "A" & .$other_allele == "T") & .$maf >= 0.4 |
         (.$effect_allele == "T" & .$other_allele == "A") & .$maf >= 0.4 |
         (.$effect_allele == "C" & .$other_allele == "G") & .$maf >= 0.4 |
         (.$effect_allele == "G" & .$other_allele == "C") & .$maf >= 0.4)
      ) ->> qc_summary[index, "ambig_filter"];.} %>%
      filter(
        !((effect_allele == "A" & other_allele == "T") & maf >= 0.4 |
          (effect_allele == "T" & other_allele == "A") & maf >= 0.4 |
          (effect_allele == "C" & other_allele == "G") & maf >= 0.4 |
          (effect_allele == "G" & other_allele == "C") & maf >= 0.4
        )
      ) %>%
      
      ## record SNPs that failed merging with the reference file
      left_join(ref_file, c("variant_id" = "variant_id")) %>%
      {sum(is.na(.$A1)) ->> qc_summary[index, "failed_merge_ref"];.} %>%
      {print(paste0(sum(is.na(.$A1)), " snps that failed merging with the ref_file (not removed)"));.} %>%
      # 0 = not present in ref_file, 1 = present in ref_file
      mutate(in_ref = ifelse(is.na(A1), 0, 1)) %>%
      
      ## remove SNPs where alleles frequencies differ wildly
      mutate(effect_allele_frequency_ref = ifelse(effect_allele == A1, MAF, 1 - MAF)) %>%
      {sum(!(abs(.$effect_allele_frequency - .$effect_allele_frequency_ref) <= 0.2), na.rm = T) ->> qc_summary[index, "frq_filter"];.} %>%
      {print(paste0("removed ", sum(!(abs(.$effect_allele_frequency - .$effect_allele_frequency_ref) <= 0.2), na.rm = T), " snps whose alleles freqs don't match the ref_file"));.} %>%
      filter(is.na(.$A1) | (abs(effect_allele_frequency - effect_allele_frequency_ref) <= 0.2)) 
    
    ## essential columns
    columns <-
      c(
        "variant_id",
        "chromosome",
        "base_pair_location",
        "other_allele",
        "effect_allele",
        "beta",
        "standard_error",
        "p_value",
        "effect_allele_frequency",
        "maf",
        "neff",
        "n_case",
        "n_control",
        "in_ref",
        "rsid"
      )
    
    if ("flip" %in% colnames(sumstats)) {
      columns <- c(columns, "flip")
    }
    
    ## for combined (case-control + proxy) samples
    if ("n_proxy_case" %in% colnames(sumstats)) {
      columns <- c(columns, "n_proxy_case", "n_proxy_control")
    }
    
    ## check if dosage or hard calls were used, and set info filter accordingly
    if (grepl("dosage", tolower(dosage_or_hardcalls), fixed = TRUE)) {
      info_thresh <- 0.6
    } else {
      info_thresh <- 0.8
    }
    
    ## info filter
    if ("info" %in% colnames(temp)) {
      temp <- temp %>%
        {sum(!(.$info >= info_thresh), na.rm = T) ->> qc_summary[index, "info_filter"];.} %>%
        {print(paste0("removed ", sum(!(.$info >= info_thresh), na.rm = T), " snps where info < ", info_thresh));.} %>%
        filter(.$info >= info_thresh)
      
      columns <- c(columns, "info")
    } else {
      temp %>%
        {"NA" ->> qc_summary[index, "info_filter"];.}
    }

    ## add rsids
    if (!"rsid" %in% colnames(sumstats)) {
      ## load g1000_afr (largest set of variants) to get rsids
      rsid_file <- fread(here( "data", "reference_data", "g1000_afr", "g1000_afr.frq")) %>%
        select(CHR, POS, SNP) %>%
        rename(rsid = SNP)
      
      temp <- temp %>%
        left_join(rsid_file, by = c("chromosome" = "CHR", "base_pair_location" = "POS"))
    }
    
    ## double betas and SEs for proxy cohorts
    # this ensures the effect sizes are on the same scale between true
    # case-control and proxy GWAS (see Marioni et al. 2018 and Liu et al. 2017)
    if (proxy_casecontrol == "proxy") {
      temp <- temp %>%
        mutate(beta = beta * 2,
               standard_error = standard_error * 2)
    }
    
    ## change precision of p_value if some are = 0
    if (any(temp$p_value == 0)) {
      
      # package to change precision
      library(Rmpfr)
      
      temp_subset <- temp %>%
        filter(p_value == 0)
      
      # loop through p_values and change precision
      for (i in 1:nrow(temp_subset)) {
        p_low <- 2*Rmpfr::pnorm(mpfr((abs(temp_subset$beta[i]) / temp_subset$standard_error[i]), precBits=100), lower.tail=FALSE, log.p = FALSE)
        p_low <- capture.output(str(p_low, give.head=FALSE, digits.d=5))
        p_low <- gsub(" ", "", p_low)
        
        temp_subset$p_low[i] <- p_low
      }
      
      temp_subset <- temp_subset %>%
        select(variant_id, p_low)
      
      temp <- temp %>%
        left_join(temp_subset, by = c("variant_id")) %>%
        mutate(p_value = ifelse(!is.na(p_low), p_low, p_value))
    }
    
    ## select columns
    clean_sumstats <- temp %>%
      select(all_of(columns))
    
    qc_summary[index, "cleaned_number"] <- nrow(clean_sumstats)
        
    ## save sumstats file
    print("saving sumstats")
    fwrite(
      x = clean_sumstats,
      file = here(
        outdir,
        paste0("clean_", outfile, ".txt.gz")
      ),
      sep = "\t",
      na = "NA",
      quote = FALSE)
    
    ## save QC summary
    print(qc_summary)
    
    qc_summary %>%
      arrange(file) %>%
      fwrite(
        file = here("data", "sumstats", "qc_summary", paste0("qc_summary_", cohort, ".txt")),
        sep = "\t",
        na = "NA",
        quote = FALSE
        )
  } else {
    print(paste0(
      here(outdir, paste0("clean_", outfile, ".txt.gz")),
      ": already exists!"
    ))
  }
}
