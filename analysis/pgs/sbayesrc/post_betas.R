#### 0. set-up ####
library(data.table)
library(dplyr)
library(here)
library(stringr)
source(here("R/snp_modifyBuild_local.R"))


#### 1. adjust posterior mean beta sumstats (hg38) ####
for (cohort in c("demge", "gothe", "twing", "xstsa", "jpbyu", "main")) {
  
  if (cohort %in% c("demge", "gothe", "twing", "xstsa")) {
    ancs <- "eur"
  } else if (cohort == "jpbyu") {
    ancs <- "eas"
  } else if (cohort == "main") {
    ancs <- c("eur", "eas", "afr")
  }
  
  for (anc in ancs) {
        
      for (q in c("_q0.9", "_q0.95")) { 
        
        if (anc != "afr") {q <- ""}
        
        ## only run if one of the output files don't exist
        if (
          !file.exists(here("analysis/pgs/sbayesrc/post_betas",paste0(cohort,"_combined_",anc,"_neff0.6_nsumstats1_lia_rc_multichain", q, "_idhg38_nayapoe.snpRes"))) | 
          !file.exists(here("analysis/pgs/sbayesrc/post_betas",paste0(cohort,"_combined_",anc,"_neff0.6_nsumstats1_lia_rc_multichain", q, "_idhg38_yeaapoe.snpRes")))
        ) {
        
          print(paste0("adjust posterior means for ", cohort, "_", anc))
          
          infile <- here(
            "analysis/pgs/sbayesrc/jz",
            paste0(
              cohort,
              "_combined_",
              anc,
              "_neff0.6_nsumstats1_lia_rc_multichain", 
              q, 
              ".snpRes.gz"
            )
          )
          
          post_betas <- fread(here(infile)) %>%
            rename(
              chr = Chrom,
              pos = Position) %>%
            mutate(variant_id_hg19 = paste(chr, pos, pmin(A1, A2), pmax(A1, A2), sep = ":")) %>%
            snp_modifyBuild_local(
              liftOver = here("src", "liftOver"),
              from = "hg19",
              to = "hg38",
              check_reverse = T
            ) %>%
            mutate(
              variant_id = paste(
                chr,
                pos,
                pmin(A1, A2),
                pmax(A1, A2),
                sep = ":"
              )
            )
          
          post_betas %>%
            ## remove APOE locus
            filter(!(chr == 19 & pos > (44908684 - 5e6) & pos < (44908684 + 5e6))) %>%
            fwrite(
              file = here(
                "analysis/pgs/sbayesrc/post_betas", 
                paste0(cohort, "_combined_", anc, "_neff0.6_nsumstats1_lia_rc_multichain", q, "_idhg38_nayapoe.snpRes")),
              quote = F,
              sep = "\t",
              col.names = T,
              na = "-9"
            )
          
          post_betas %>%
            fwrite(
              file = here(
                "analysis/pgs/sbayesrc/post_betas", 
                paste0(cohort, "_combined_", anc, "_neff0.6_nsumstats1_lia_rc_multichain", q, "_idhg38_yeaapoe.snpRes")),
              quote = F,
              sep = "\t",
              col.names = T,
              na = "-9"
            ) 
      }
    }
  }
}


#### 2. adjust posterior mean beta sumstats (hg37) ####
for (cohort in c("xukbb", "main")) {
  
  if (cohort == "xukbb") {
    ancs <- c("eur", "afr")
  } else if (cohort == "main") {
    ancs <- c("eur", "eas", "afr")
  }
  
  for (anc in ancs) {
    
    for (q in c("_q0.9", "_q0.95")) { 
      
      if (anc != "afr") {q <- ""}
      
      ## only run if one of the output files don't exist
      if (
        !file.exists(here("analysis/pgs/sbayesrc/post_betas",paste0(cohort, "_combined_", anc, "_neff0.6_nsumstats1_lia_rc_multichain", q, "_idhg37_yeaapoe.snpRes"))) | 
        !file.exists(here("analysis/pgs/sbayesrc/post_betas",paste0(cohort, "_combined_", anc, "_neff0.6_nsumstats1_lia_rc_multichain", q, "_idhg37_nayapoe.snpRes")))
      ) {
        
        print(paste0("adjust posterior means for ", cohort, "_", anc))
        
        infile <- here(
          "analysis/pgs/sbayesrc/jz",
          paste0(
            cohort,
            "_combined_",
            anc,
            "_neff0.6_nsumstats1_lia_rc_multichain", 
            q, 
            ".snpRes.gz"
          )
        )
        
        post_betas <- fread(here(infile)) %>%
          mutate(
            variant_id = paste(
              Chrom,
              Position,
              pmin(A1, A2),
              pmax(A1, A2),
              sep = ":"
            )
          )
        
        post_betas %>%
          ## remove APOE locus
          filter(!(Chrom == 19 & Position > (45411941 - 5e6) & Position < (45411941 + 5e6))) %>%
          fwrite(
            file = here(
              "analysis/pgs/sbayesrc/post_betas", 
              paste0(cohort, "_combined_", anc, "_neff0.6_nsumstats1_lia_rc_multichain", q, "_idhg37_nayapoe.snpRes")),
            quote = F,
            sep = "\t",
            col.names = T,
            na = "-9"
          )
        
        post_betas %>%
          fwrite(
            file = here(
              "analysis/pgs/sbayesrc/post_betas", 
              paste0(cohort, "_combined_", anc, "_neff0.6_nsumstats1_lia_rc_multichain", q, "_idhg37_yeaapoe.snpRes")),
            quote = F,
            sep = "\t",
            col.names = T,
            na = "-9"
          ) 
      }
    }
  }
}

