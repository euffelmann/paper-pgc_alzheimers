run_metal <- function(sumstats,
                      metal_dir = here("src/METAL-2020-05-05/build/bin/metal"),
                      rsid_file = here("data/reference_data/g1000_afr/g1000_afr.frq"),
                      out_dir = ".",
                      outfile,
                      scheme = "STDERR",
                      average_frq = "ON",
                      opts,
                      cols = "beta",
                      filter_sumstats = F,
                      neff_threshold = 0,
                      min_sumstats = 2) {
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  ### run genome-wide association meta analysis with metal
  ### ref: Willer et al. 2010 (Bioinformatics)
  
  ## - sumstats:        vector of full sumstats file names
  ## - metal_dir:       path to metal .exe
  ## - out_dir:         output sumstats
  ## - outfile:         file name of output
  ## - scheme:          standard error- or N-weighted meta analysis
  ##                    default is standard error weighted
  ##                    {SAMPLESIZE | STDERR}
  ## - average_frq:     calculate average allele frequency across input sumstats
  ##                    {ON | OFF}
  ## - opts:            additional input to metal, such as to filter input or 
  ##                    create an additional variable in the output
  ## - cols:            use betas or z-scores as input for N-weighted meta-analysis
  ##                    {beta | z_score}
  ## - neff_threshold:  percentage of total neff. SNPs with neff below this
  ##                    will be removed
  ## - min_cohorts:     SNPs that are observed in fewer cohorts than min_cohorts
  ##                    will be removed

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  ## set-up (load libraries and install when necessary)
  if (!require(here)) { install.packages("here"); library(here) }
  if (!require(dplyr)) { install.packages("dplyr"); library(dplyr) }
  if (!require(data.table)) { install.packages("data.table"); library(data.table) }
  if (!require(stringr)) { install.packages("stringr"); library(stringr) }
  
  if (!file.exists(paste0(out_dir, "/", outfile, ".txt.gz"))) {
    
    #### metal settings ####
    system(
      paste0(
        "
        echo -e '
        SCHEME ", scheme," 
        ", opts, "
        AVERAGEFREQ ", average_frq, "
        ' > ", out_dir, "/temp_metal.sh"
      )
    ) 
    
    if (scheme == "STDERR") {
      
      cols <- list(
        "MARKER" = "variant_id", 
        "WEIGHT" = "neff", 
        "ALLELE" = "effect_allele other_allele", 
        "FREQ" = "effect_allele_frequency", 
        "EFFECT" = "beta",
        "STDERR" = "standard_error",
        "PVAL" = "p_value"
      )
      
    } else if (scheme == "SAMPLESIZE" & cols == "z_score") {
      
      cols <- list(
        "MARKER" = "variant_id", 
        "WEIGHT" = "neff", 
        "ALLELE" = "effect_allele other_allele", 
        "FREQ" = "effect_allele_frequency", 
        "EFFECT" = "z_score",
        "PVAL" = "p_value"
      )
      
    } else if (scheme == "SAMPLESIZE" & cols == "beta") {
      
      cols <- list(
        "MARKER" = "variant_id", 
        "WEIGHT" = "neff", 
        "ALLELE" = "effect_allele other_allele", 
        "FREQ" = "effect_allele_frequency", 
        "EFFECT" = "beta",
        "PVAL" = "p_value"
      )
      
    }
    
    #### process sumstats ####
    for (file in sumstats) {
      
      if (scheme == "STDERR") {
        
        system(
          paste0(
            "
          echo -e '
          MARKER  ", cols[["MARKER"]] , "
          WEIGHT  ", cols[["WEIGHT"]] , "
          ALLELE  ", cols[["ALLELE"]] , "
          FREQ    ", cols[["FREQ"]] , "
          EFFECT  ", cols[["EFFECT"]] , "
          STDERR  ", cols[["STDERR"]] , "
          PVAL    ", cols[["PVAL"]] , "
          
          PROCESS ", file, " ' ", " >> ", out_dir, "/temp_metal.sh"
          )
        )
        
      } else {
        
        system(
          paste0(
            "
          echo -e '
          MARKER  ", cols[["MARKER"]] , "
          WEIGHT  ", cols[["WEIGHT"]] , "
          ALLELE  ", cols[["ALLELE"]] , "
          FREQ    ", cols[["FREQ"]] , "
          EFFECT  ", cols[["EFFECT"]] , "
          PVAL    ", cols[["PVAL"]] , "
          
          PROCESS ", file, " ' ", " >> ", out_dir, "/temp_metal.sh"
          )
        )
        
      }
      
    }
    
    #### run metal ####
    system(
      paste0(
        "
        echo -e '
        OUTFILE ", out_dir, "/", outfile, " .txt 
        ANALYZE
        QUIT
        
        ' >> ", out_dir, "/temp_metal.sh",
        "
      ", metal_dir, " ",  out_dir, "/temp_metal.sh"
      )
    )
    
    ## remove metal script
    system(paste0("rm ", out_dir, "/temp_metal.sh"))
    
    ## remove "1" from outfile
    system(paste0("mv ", out_dir, "/", outfile, "1.txt ", out_dir, "/", outfile, ".txt"))
    system(paste0("mv ", out_dir, "/", outfile, "1.txt.info ", out_dir, "/", outfile, ".txt.info"))
    
    
    #### sanity checks ####
    
    ## check correct input files
    info_files <- system(paste0("grep '# --> Input File' ", out_dir, "/", outfile, ".txt.info"),
                         intern = T)
    n_files <- length(info_files)
    for (i in 1:n_files) {
      info_files[i] <- gsub(pattern = paste0("# --> Input File ", i, " : "),
                            replacement = "",
                            info_files[i])
      error_messages <- NULL
      if (!info_files[i] %in% sumstats) {
        error_messages[i] <- paste0("This file should not have been meta-analyzed: ", info_files[i])
      }
      if (!is.null(error_messages)) {
        stop(print(error_messages))
      }
    }
    
    
    #### re-format output  ####
    print("reformatting sumstats ...")
    
    if (scheme == "STDERR") {
      reformated_sumstats <- fread(paste0(out_dir, "/", outfile, ".txt")) %>%
        mutate(
          chromosome = as.integer(sub(":.*", "", MarkerName)),
          base_pair_location = as.integer(sub("^[^:]*:([^:]*):.*", "\\1", MarkerName)),
          Allele2 = toupper(Allele2),
          Allele1 = toupper(Allele1)
        ) %>%
        rename(
          other_allele = Allele2,
          effect_allele = Allele1,
          effect_allele_frequency = Freq1,
          beta = Effect,
          standard_error = StdErr,
          p_value = `P-value`,
          direction = Direction
        ) %>%
        mutate(
          variant_id = paste0(chromosome, ":", base_pair_location, ":",
                              pmin(effect_allele, other_allele), ":",
                              pmax(effect_allele, other_allele))
        ) %>%
        select(
          -MarkerName,
          -FreqSE
        ) %>%
        relocate(
          variant_id,
          chromosome,
          base_pair_location,
          other_allele
        ) %>%
        arrange(chromosome, base_pair_location)
    } else if (scheme == "SAMPLESIZE") {
      reformated_sumstats <- fread(paste0(out_dir, "/", outfile, ".txt")) %>%
        mutate(
          chromosome = as.integer(sub(":.*", "", MarkerName)),
          base_pair_location = as.integer(sub("^[^:]*:([^:]*):.*", "\\1", MarkerName)),
          Allele2 = toupper(Allele2),
          Allele1 = toupper(Allele1)
        ) %>%
        rename(
          other_allele = Allele2,
          effect_allele = Allele1,
          effect_allele_frequency = Freq1,
          z_score = Zscore,
          p_value = `P-value`,
          direction = Direction
        ) %>%
        mutate(
          variant_id = paste0(chromosome, ":", base_pair_location, ":",
                              pmin(effect_allele, other_allele), ":",
                              pmax(effect_allele, other_allele))
        ) %>%
        select(
          -MarkerName,
          -FreqSE,
          -Weight
        ) %>%
        relocate(
          variant_id,
          chromosome,
          base_pair_location,
          other_allele
        ) %>%
        arrange(chromosome, base_pair_location)
    }
    
    
    ## add rsids
    # load g1000_afr (largest set of variants) to get rsids
    print("adding rsids ...")
    rsid_file <- fread(rsid_file) %>%
      select(CHR, POS, SNP) %>%
      rename(rsid = SNP)
    
    reformated_sumstats <- reformated_sumstats %>%
      left_join(rsid_file, by = c("chromosome" = "CHR", "base_pair_location" = "POS"))
    
    #### filter meta-analyzed sumstats  ####
    if (length(sumstats) > 2 & filter_sumstats == T) {
      print("filtering sumstats ...")
      
      ## filter SNPs
      pre_filter_snps <- NULL
      filter_snps <- NULL
      post_filter_snps <- NULL
      reformated_sumstats <- reformated_sumstats %>%
        {nrow(.) ->> pre_filter_snps;.} %>%
        mutate(n_sumstats = nchar(direction) - str_count(direction, "\\?")) %>%
        {sum(!((.$neff >= round(neff_threshold * max(.$neff))) & .$n_sumstats >= min_sumstats)) ->> filter_snps;.} %>%
        filter((neff >= round(neff_threshold * max(.$neff))) & (n_sumstats >= min_sumstats)) %>%
        {nrow(.) ->> post_filter_snps;.}
      
      ## record number of SNPs
      system(
      paste0("

echo -e 'Number of SNPs before filtering: ", pre_filter_snps, " SNPs
Number of SNPs with neff < ", neff_threshold * 100, "% of max(neff) and that are not present in at least ", min_sumstats, " input sumstats: ", filter_snps, " SNPs
Number of SNPs after filtering: ", post_filter_snps, " SNPs
' >> ", out_dir, "/", outfile, ".txt.info"))
      
    }
    
    
    #### save sumstats ####
    print(paste0("saving sumstats: ", out_dir, "/", outfile, ".txt.gz"))
    reformated_sumstats %>%
      fwrite(
        file = paste0(out_dir, "/", outfile, ".txt.gz"),
        sep = "\t",
        na = "NA",
        quote = FALSE
      )
    
    ## remove non-formatted sumstats files
    system(paste0("rm ", out_dir, "/", outfile, ".txt"))

  } else {
    
    print(paste0("file already exists: ", out_dir, "/", outfile, ".txt.gz"))
  
  }
  
}


