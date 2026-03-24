load_and_qc_sumstats <- function(sumstats_file) {
  ## function to clean sumstats to prevent errors when plotting raw sumstats
  
  if (!exists("sumstats")) {
    sumstats <- fread(sumstats_file)
    
    if (sumstats$chromosome[1] == "chr1") {
      sumstats$chromosome <- as.integer(sub("chr", "", sumstats$chromosome))
    }
    
    if ("X" %in% sumstats$chromosome) {
      sumstats$chromosome[sumstats$chromosome == "X"] <- 23
      sumstats$chromosome <- as.integer(sumstats$chromosome)
    }
    
    sumstats <- sumstats[!is.na(sumstats$base_pair_location), ]
    
    sumstats$p_value <- as.numeric(sumstats$p_value)
    
    return(sumstats)
  }
}

qc_plots <- function(sumstats_file, outfile, qc_dir, ancestry) {
  ## function to several qc plots for visual inspection
  ## input:
  ##    sumstats_file:  raw sumstats data frame in GWAS catalog format
  ##    outfile:        name of output file
  ##    qc_dir:         directory to save output files
  ##    ancestry:       reference file with allele frequencies
  
  ## set variables
  width <- 30
  height <- 15
  res <- 300
  
  print(paste0("Evaluating: ", outfile))
  
  #### qq plot ####
  if (!file.exists(here(qc_dir, paste0(outfile, "_qq_plot.png")))) {
    print(paste0("plotting qq for: ", outfile))
    if (!exists("sumstats")) {
      sumstats <- load_and_qc_sumstats(sumstats_file)
    }
    
    png(
      filename = here(
        qc_dir,
        paste0(outfile, "_qq_plot.png")
      ),
      width = height,
      height = height,
      units = "cm",
      res = res
    )
    qq(sumstats$p_value[sumstats$p_value > 0])
    dev.off()
  }
  
  #### manhattan plot ####
  if (!file.exists(here(qc_dir, paste0(outfile, "_manhattan_plot.png")))) { 
    print(paste0("plotting manhattan for: ", outfile))
    
    if (!exists("sumstats")) {
      sumstats <- load_and_qc_sumstats(sumstats_file)
    }
    
    png(
      filename = here(
        qc_dir,
        paste0(outfile, "_manhattan_plot.png")
      ),
      width = width,
      height = height,
      units = "cm",
      res = res
    )
    manhattan(
      sumstats[sumstats$p_value > 0 & sumstats$p_value < 1, ],
      chr = "chromosome",
      bp = "base_pair_location",
      snp = "variant_id",
      p = "p_value"
    )
    dev.off()
    }
  
  #### histogram of BETAs ###
  if (!file.exists(here(qc_dir, paste0(outfile, "_beta_histogram_plot.png")))) { 
    print(paste0("plotting betas for: ", outfile))
    
    if (!exists("sumstats")) {
      sumstats <- load_and_qc_sumstats(sumstats_file)
    }
    
    png(
      filename = here(
        qc_dir,
        paste0(outfile, "_beta_histogram_plot.png")
      ),
      width = width,
      height = height,
      units = "cm",
      res = res
    )
    hist(sumstats$beta)
    dev.off()
  }

  #### boxplots of BETAs per chromosome ####
  if (!file.exists(here(qc_dir, paste0(outfile, "_beta_box_per_chr_plot.png")))) { 
    print(paste0("plotting betas per chr for: ", outfile))
    if (!exists("sumstats")) {
      sumstats <- load_and_qc_sumstats(sumstats_file)
    }
    
    png(
      filename = here(
        qc_dir,
        paste0(outfile, "_beta_box_per_chr_plot.png")
      ),
      width = width,
      height = height,
      units = "cm",
      res = res
    )
    boxplot(sumstats$beta ~ sumstats$chromosome)
    dev.off()
  }
  
  #### boxplot of SEs per chromosome ####
  if (!file.exists(here(qc_dir, paste0(outfile, "_se_box_per_chr_plot.png")))) { 
    print(paste0("plotting SEs for: ", outfile))
    
    if (!exists("sumstats")) {
      sumstats <- load_and_qc_sumstats(sumstats_file)
    }
    
    png(
      filename = here(
        qc_dir,
        paste0(outfile, "_se_box_per_chr_plot.png")
      ),
      width = width,
      height = height,
      units = "cm",
      res = res
    )
    boxplot(sumstats$standard_error ~ sumstats$chromosome)
    dev.off()
  }

  #### histogram of SEs ####
  if (!file.exists(here(qc_dir, paste0(outfile, "_se_histogram_plot.png")))) { 
    print(paste0("plotting SEs per chr for: ", outfile))
    if (!exists("sumstats")) {
      sumstats <- load_and_qc_sumstats(sumstats_file)
    }
    
    png(
      filename = here(
        qc_dir,
        paste0(outfile, "_se_histogram_plot.png")
      ),
      width = width,
      height = height,
      units = "cm",
      res = res
    )
    hist(sumstats$standard_error)
    dev.off()
  }
  
  #### histogram of P values ####
  if (!file.exists(here(qc_dir, paste0(outfile, "_p_histogram_plot.png")))) { 
    print(paste0("plotting p_values for: ", outfile))
    
    if (!exists("sumstats")) {
      sumstats <- load_and_qc_sumstats(sumstats_file)
    }
    
    png(
      filename = here(
        qc_dir,
        paste0(outfile, "_p_histogram_plot.png")
      ),
      width = width,
      height = height,
      units = "cm",
      res = res
    )
    hist(sumstats$p_value)
    dev.off()
  }
  
  #### histogram of maf ####
  if (!file.exists(here(qc_dir, paste0(outfile, "_maf_histogram_plot.png")))) { 
    print(paste0("plotting mafs for: ", outfile))
    
    if (!exists("sumstats")) {
      sumstats <- load_and_qc_sumstats(sumstats_file)
    }
    
    # get maf
    temp <- sumstats %>%
      rowwise() %>%
      mutate(maf = min(effect_allele_frequency, 1-effect_allele_frequency, na.rm = TRUE)) 
    png(
      filename = here(
        qc_dir,
        paste0(outfile, "_maf_histogram_plot.png")
      ),
      width = width,
      height = height,
      units = "cm",
      res = res
    )
    hist(temp$maf)
    dev.off()
  }
  
  #### histogram of info ####
  if (!file.exists(here(qc_dir, paste0(outfile, "_info_histogram_plot.png")))) { 
    
    ## check if info is available
    header <- system(
      paste0(
        "gunzip -c ", sumstats_file, " | head -1"
      ), 
      intern = T
    )
    
    if (grepl("info", header, fixed = TRUE)) {
      print(paste0("plotting info scores for: ", outfile))
      
      if (!exists("sumstats")) {
        sumstats <- load_and_qc_sumstats(sumstats_file)
      }
      
      png(
        filename = here(
          qc_dir,
          paste0(outfile, "_info_histogram_plot.png")
        ),
        width = width,
        height = height,
        units = "cm",
        res = res
      )
      hist(sumstats$info)
      dev.off()
      
    }
  }
  
  #### histogram of neff ####
  if (!file.exists(here(qc_dir, paste0(outfile, "_neff_histogram_plot.png")))) { 
    print(paste0("plotting neffs for: ", outfile))
    
    if (!exists("sumstats")) {
      sumstats <- load_and_qc_sumstats(sumstats_file)
    }
    
    png(
      filename = here(
        qc_dir,
        paste0(outfile, "_neff_histogram_plot.png")
      ),
      width = width,
      height = height,
      units = "cm",
      res = res
    )
    hist(sumstats$neff)
    dev.off()
  }
  
  ####  P-Z plot ####
  if (!file.exists(here(qc_dir, paste0(outfile, "_pz_plot.png")))) { 
    print(paste0("plotting pz_plot for: ", outfile))
    
    if (!exists("sumstats")) {
      sumstats <- load_and_qc_sumstats(sumstats_file)
    }
    
    temp <- sumstats %>%
      mutate(
        p_value_z = 2 * pnorm(q = (abs(sumstats$beta) / sumstats$standard_error), lower.tail = F)
        )
    png(
      filename = here(
        qc_dir,
        paste0(outfile, "_pz_plot.png")
      ),
      width = height,
      height = height,
      units = "cm",
      res = res
    )
    plot(-log10(temp$p_value_z), -log10(temp$p_value))
    abline(a = 0, b = 1)
    dev.off()
  }
  
  #### effect_allele_frequency plot (eaf_plot) ####
  
  ## determine which reference panels to use
  if (ancestry == "eur" | ancestry == "all") { reference_panels <- c("g1000_eur", "hrc_eur") }
  if (ancestry == "eas") { reference_panels <- c("g1000_eas", "jpb_eas") } 
  if (ancestry == "afr") { reference_panels <- c("g1000_afr", "agr_afr") } 
  if (ancestry == "amr") { reference_panels <- c("g1000_amr", "g1000_eur") } 
  if (ancestry == "sas") { reference_panels <- c("g1000_sas") } 
  
  if (exists("reference_panels")) {
    
    for (reference_panel in c(reference_panels)) {
      
      if (!file.exists(here(qc_dir, paste0(outfile, "_", reference_panel, "_eaf_plot.png")))) { 
        
        print(paste0("plotting eaf_plot for: ", outfile, "_", reference_panel, "_eaf_plot.png"))
        
        ref_file <-
          fread(here(
            "data",
            "reference_data",
            reference_panel,
            paste0(reference_panel, ".frq")
          )) %>%
          mutate(variant_id = paste0(CHR, ":", POS, ":", pmin(A1, A2), ":", pmax(A1, A2)))
        
        if (!exists("sumstats")) {
          sumstats <- load_and_qc_sumstats(sumstats_file)
        }
        
        temp <- sumstats %>%
          left_join(ref_file, by = c("variant_id" = "variant_id")) %>%
          mutate(effect_allele_frequency_ref = ifelse(effect_allele == A1, MAF, 1 - MAF))
        
        png(
          filename = here(
            qc_dir,
            paste0(outfile, "_", reference_panel, "_eaf_plot.png")),
          width = height,
          height = height,
          units = "cm",
          res = res
        )
        
        plot(temp$effect_allele_frequency_ref, temp$effect_allele_frequency)
        abline(a = 0, b = 1, col = "red")
        abline(a = 0.2, b = 1, col = "orange")
        abline(a = -0.2, b = 1, col = "orange")
        dev.off()
      }
    }
  }
}


