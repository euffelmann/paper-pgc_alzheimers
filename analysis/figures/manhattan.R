#### 0. setup ####
library(here)
library(data.table)
library(ggplot2)
library(dplyr)
library(openxlsx)
source(here("R/manhattan_f.R"))

#### 1. individual manhattan plots ####
for (anc in c("all", "eur", "afr", "eas", "sas", "amr")) {
  
  for (casec in c("combined", "case_control", "proxy")) {
    
    if (casec == "combined") {casec_out <- "combi"}
    if (casec == "case_control") {casec_out <- "casec"}
    if (casec == "proxy") {casec_out <- "proxy"}
    
    if (anc == "all") {anc_loci <- "All"} else {anc_loci <- toupper(anc)}
    
    if (!file.exists(here("analysis/figures/manhattans", paste0("manhattan_main_", casec_out, "_", anc, "_flames.png"))) &
         file.exists(here("data/sumstats/meta/main", casec, paste0("main_", casec, "_", anc, "_neff0.6_nsumstats1.txt.gz")))) {
      
      print(paste0("loading sumstats: main_", casec, "_", anc, "_neff0.6_nsumstats1.txt.gz"))
      
      sst <- fread(here("data/sumstats/meta/main", casec, paste0("main_", casec, "_", anc, "_neff0.6_nsumstats1.txt.gz")))
      
      if (any(sst$p_value <= 5e-8)) {
        
        ### # # # #  gene names # # # #
        selected_genes <- read.xlsx(here("analysis/gene_table/manual_selection/2025_06_23_gene_table.xlsx")) %>%
          dplyr::select(locus, manual_genes, novel) %>%
          mutate(
            gene = sapply(strsplit(manual_genes, "/"), function(genes) {
              if(length(genes) > 2) {
                paste0(paste(genes[1:2], collapse = "/"), "/...")
              } else {
                paste(genes, collapse = "/")
              }
            }
            )
          ) %>%
          dplyr::select(-manual_genes)
        
        loci_info <- readRDS(here("analysis/risk_loci/risk_loci_table.rds"))[["unique_loci_main_combined"]] %>%
          dplyr::select(locus, index_variant_id, ancestry)
        
        genes <- selected_genes %>%
          left_join(loci_info, by = "locus") %>%
          filter(grepl(x = ancestry, anc)) %>%
          rowwise() %>%
          mutate(
            variant_id = {
              ancestries <- stringr::str_split(ancestry, ";")[[1]]
              variants <- stringr::str_split(index_variant_id, ";")[[1]]
              idx <- which(ancestries == anc)
              if (length(idx) > 0 && idx <= length(variants)) {
                variants[[idx]]
              } else {
                NA_character_
              }
            }
          ) %>%
          ungroup()
        ### # # # # # # # # # # # # # # #
       
        if (!file.exists(here("analysis/figures/manhattans", paste0("manhattan_main_", casec_out, "_", anc, "_with_genes.pdf"))) |
            !file.exists(here("analysis/figures/manhattans", paste0("manhattan_main_", casec_out, "_", anc, "_with_genes.png"))) |
            !file.exists(here("analysis/figures/manhattans", paste0("manhattan_main_", casec_out, "_", anc, "_no_genes.pdf"))) |
            !file.exists(here("analysis/figures/manhattans", paste0("manhattan_main_", casec_out, "_", anc, "_no_genes.png")))) {
          
          print(paste0("making plots for: manhattan_main_", casec_out, "_", anc))
          
          p <- manhattan_f(sst, annotate_loci = T, loci = genes)
          
          ggsave(
            plot = p,
            filename = here("analysis/figures/manhattans", paste0("manhattan_main_", casec_out, "_", anc, "_with_genes.pdf")),
            device = cairo_pdf,
            width = 8,
            height = 4
          )
          
          ggsave(
            plot = p,
            filename = here("analysis/figures/manhattans", paste0("manhattan_main_", casec_out, "_", anc, "_with_genes.png")),
            device = "png",
            width = 8,
            height = 4
          )
          
          p <- manhattan_f(sst, annotate_loci = F)
          
          ggsave(
            plot = p,
            filename = here("analysis/figures/manhattans", paste0("manhattan_main_", casec_out, "_", anc, "_no_genes.pdf")),
            device = cairo_pdf,
            width = 8,
            height = 4
          )
          
          ggsave(
            plot = p,
            filename = here("analysis/figures/manhattans", paste0("manhattan_main_", casec_out, "_", anc, "_no_genes.png")),
            device = "png",
            width = 8,
            height = 4
          ) 
          
          
        } else {
          
          print(paste0("plots already exist for: manhattan_main_", casec_out, "_", anc))
          
        }
        
      } else {
        
        if (!file.exists(here("analysis/figures/manhattans", paste0("manhattan_main_", casec_out, "_", anc, "_no_genes.pdf"))) |
            !file.exists(here("analysis/figures/manhattans", paste0("manhattan_main_", casec_out, "_", anc, "_no_genes.png")))) {
          
          print(paste0("making plots for: manhattan_main_", casec_out, "_", anc))
          
          p <- manhattan_f(sst, annotate_loci = F)
          
          ggsave(
            plot = p,
            filename = here("analysis/figures/manhattans", paste0("manhattan_main_", casec_out, "_", anc, "_no_genes.pdf")),
            device = cairo_pdf,
            width = 8,
            height = 4
          )
          
          ggsave(
            plot = p,
            filename = here("analysis/figures/manhattans", paste0("manhattan_main_", casec_out, "_", anc, "_no_genes.png")),
            device = "png",
            width = 8,
            height = 4
          ) 
          
        } else {
          
          print(paste0("plots already exist for: manhattan_main_", casec_out, "_", anc))
          
        }
      } 
    }
  }
}


