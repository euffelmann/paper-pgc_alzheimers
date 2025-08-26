#### 00. setup ####
library(data.table)
library(dplyr)
library(here)

extract_snps <- function(locus_nr, loci, analysis, ancestry, proxy_cc_combined, params) {
    
  variant_id <- loci$index_variant_id[loci$locus == locus_nr]
  plot_sst_dir <- here("analysis/figures/forest_plots", analysis, proxy_cc_combined, ancestry, "sumstats")
  plot_sst_file <- paste0(gsub(pattern = ":", replacement = "_", x = variant_id), ".txt")
  
  if (!file.exists(paste0(plot_sst_dir, "/", plot_sst_file))) {
    
    snps <- NULL
    for (j in 1:nrow(params)) {
      
      file <- params$sumstats_file[j]
      name <- params$name[j]; print(name)
      casec_proxy <- params$casec_proxy[j]
      ancestry <- params$ancestry[j]
      
      temp <- fread(cmd = paste0("gunzip -c ", file, " | grep '", variant_id, "'"))
      if (nrow(temp) == 0) {
        
        temp <- data.frame(
          variant_id = variant_id,
          chromosome = loci$chromosome[loci$locus == locus_nr],
          base_pair_location = loci$base_pair_location[loci$locus == locus_nr],
          other_allele = NA,
          effect_allele = NA,
          beta = NA,
          standard_error = NA,
          p_value = NA,
          effect_allele_frequency = NA,
          maf = NA,
          neff = NA,
          n_case = NA,
          n_control = NA,
          cohort = paste(name, casec_proxy, ancestry, sep = "|")
          
        )
        
      } else {
        
        sumstats_colnames <- colnames(fread(cmd = paste0("gunzip -c ", params$sumstats_file[j], " | head -1")))
        names(temp) <- sumstats_colnames
        temp <- temp %>%
          mutate(
            cohort = paste(name, casec_proxy, ancestry, sep = "|"),
            maf = pmin(effect_allele_frequency, 1 - effect_allele_frequency)
          ) %>%
          select(
            variant_id,
            chromosome,
            base_pair_location,
            other_allele,
            effect_allele,
            beta,
            standard_error,
            p_value,
            effect_allele_frequency,
            maf,
            neff,
            n_case,
            n_control,
            cohort
          )
      }
      
      snps <- rbind(snps, temp)
      print(snps)
      
    }
    
    return(snps)
    
  } else {
    
    print(paste0("SNPs have already been extracted and saved in: ",
                 plot_sst_dir, 
                 "/", 
                 plot_sst_file)
    )
    
  }
}
save_forest <- function(locus_nr, loci, analysis, ancestry, proxy_cc_combined, meta_sumstats) {
  
  variant_id <- loci$index_variant_id[loci$locus == locus_nr]
  index_variant_id <- gsub(pattern = ":", replacement = "_", variant_id)
  
  if (!file.exists(here("analysis/figures/forest_plots", analysis, proxy_cc_combined, ancestry, paste0(index_variant_id, ".pdf")))) {
    
    sumstats <- fread(here("analysis/figures/forest_plots", analysis, proxy_cc_combined, ancestry, "sumstats", paste0(index_variant_id, ".txt")))
    
    ## align effect sizes to combined meta-analysis sumstats
    temp <- meta_sumstats %>%
      filter(variant_id == !!variant_id)
    
    for (j in 1:nrow(sumstats)) {
      
      if ((sumstats$effect_allele[j] != temp$effect_allele) & !is.na(sumstats$effect_allele[j])) {
        
        ## flip allele
        sumstats$beta[j] <- sumstats$beta[j] * -1
        sumstats$effect_allele[j] <- temp$effect_allele
        sumstats$other_allele[j] <- temp$other_allele
        sumstats$effect_allele_frequency[j] <- 1 - sumstats$effect_allele_frequency[j]
        
      }
    }
    
    m.gen <- meta::metagen(
      TE = beta,
      seTE = standard_error,
      pval = p_value,
      studlab = cohort,
      data = sumstats,
      common = TRUE,
      random = F, 
      title = sumstats$variant_id[1]
    )
    
    square_colors <- ifelse(m.gen$TE < 0, "#ca0020", "#0571b0")
    
    options(na.action = "na.pass")
    
    pdf(
      here(
        "analysis/figures/forest_plots/", analysis, proxy_cc_combined, ancestry,
        paste0(index_variant_id, ".pdf")),
      width = 12,
      height = 11,
      family = "Courier"
    )
    
    meta::forest(
      m.gen,
      prediction = FALSE,
      print.tau2 = FALSE,
      leftcols = c("studlab", "TE", "seTE", "pval"),
      leftlabs = c("Cohort", "log(OR)", "Standard Error", "P-value"),
      xlim = c(min(min(sumstats$beta, na.rm = T) - 0.1, -0.1), 
               max(max(sumstats$beta, na.rm = T) + 0.1, 0.1)),
      digits = 4,
      scientific.pval = T,
      lab.NA = "NA",
      ref = 0,
      col.square = square_colors,
      col.square.lines = square_colors,
      col.diamond = "black"
    )
    grid::grid.text(
      paste0("effect allele:", temp$effect_allele, " frequency:", temp$effect_allele_frequency, " neff:", format(temp$neff, big.mark = ",")), 
      x = 0.5, 
      y = 0.96, 
      gp = grid::gpar(fontsize = 12, fontface = "bold"))
    grid::grid.text(
      sumstats$variant_id[1], 
      x = 0.5, 
      y = 0.98, 
      gp = grid::gpar(fontsize = 16, fontface = "bold"))
    dev.off()
    
  }
}

#### 01. main ####
## (main_combined_all_neff0.6_nsumstats1.txt.gz

##### a. SNP stats from each GWAS for each meta-analysis locus #####

## individual sumstats
params_cc <- data.frame(
  sumstats_file = list.files(
    here("data", "sumstats", "case_control"),
    pattern = "\\.gz$",
    full.names = T
  ),
  ancestry = c(
    sub(".*_(eur|eas|afr|amr|sas).*", "\\1", list.files(here("data", "sumstats", "case_control"), pattern = "\\.gz$"))) 
  ) %>%
  filter(
    # filter for main analysis or those that didn't have age available
    grepl("main", sumstats_file) | 
      (
        grepl("noag", sumstats_file) & 
          (
            grepl("mvplo", sumstats_file) |
              grepl("biovu", sumstats_file) |
              grepl("grace", sumstats_file) | 
              grepl("xstsa", sumstats_file) | 
              grepl("twing", sumstats_file) | 
              grepl("gothe", sumstats_file)
          )
      )
  )

params_proxy <- data.frame(
  sumstats_file = list.files(
    here("data", "sumstats", "proxy_meta"),
    pattern = "\\.gz$",
    full.names = T
  ),
  ancestry = c(
    sub(".*_(eur|eas|afr|amr|sas).*", "\\1", list.files(here("data", "sumstats", "proxy_meta"), pattern = "\\.gz$"))) 
) %>%
  filter(
    grepl("main", sumstats_file) | grepl("mvplo", sumstats_file),
    !grepl("inref", sumstats_file)
  )

params <- rbind(params_cc, params_proxy) %>%
  mutate(
    name = sub(".*clean_([a-zA-Z0-9]{5}).*", "\\1", sumstats_file),
    casec_proxy = ifelse(grepl("case_control", sumstats_file), "casec", "proxy")
  ) %>%
  arrange(casec_proxy, factor(ancestry, levels = c("eur", "afr", "eas", "amr", "sas")), name)

analysis <- "main"
proxy_cc_combined <- "combined"

for (anc in c("all", "eur", "afr", "eas", "amr", "sas")) {
  
  locus_file <- here("analysis/risk_loci",analysis,paste(analysis, proxy_cc_combined,anc,"neff0.6_nsumstats1_loci.txt",sep = "_"))
  loci <- fread(locus_file)
  temp_params <- params %>%
    filter(ancestry == anc)
  
  ## extract index SNPs for each locus
  for (i in 1:max(loci$locus)) {
    
    snps <- extract_snps(locus_nr = i, loci = loci, analysis = analysis, ancestry = anc, proxy_cc_combined = proxy_cc_combined, params = temp_params)
    variant_id <- loci$index_variant_id[loci$locus == i]
    
    if (is.data.frame(snps))
    fwrite(
      x = snps,
      file = here(
        "analysis/figures/forest_plots/",
        analysis,
        proxy_cc_combined,
        anc,
        "sumstats",
        paste0(gsub(
          pattern = ":", replacement = "_", variant_id
        ), ".txt")
      ),
      quote = F,
      sep = "\t",
      na = NA
    )
  }
}



##### b. forest plot #####

analysis <- "main"
proxy_cc_combined <- "combined"

for (anc in c("all", "eur", "afr", "eas", "amr", "sas")) { 
  
  locus_file <- here(
    "analysis/risk_loci",
    analysis,
    paste(
      analysis,
      proxy_cc_combined,
      anc,
      "neff0.6_nsumstats1_loci.txt",
      sep = "_"
    )
  )
  loci <- fread(locus_file)
  meta_sumstats <- fread(here(
    "data/sumstats/meta/",
    analysis,
    proxy_cc_combined,
    paste0(
      analysis,
      "_",
      proxy_cc_combined,
      "_",
      anc,
      "_neff0.6_nsumstats1.txt.gz"
    )
  ))
  
  ## save plots
  for (i in 1:max(loci$locus)) { 
    
    save_forest(
      locus_nr = i,
      loci = loci,
      analysis = analysis,
      ancestry = anc,
      proxy_cc_combined = proxy_cc_combined,
      meta_sumstats = meta_sumstats
    )
  }
}


