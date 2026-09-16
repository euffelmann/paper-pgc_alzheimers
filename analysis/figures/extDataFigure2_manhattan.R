# =============================================================
#
#  manhattan.R
#  Generate Manhattan plots for all ancestry × case/control
#  strata of the main AD GWAS meta-analysis.
#
#  For each stratum, two versions are saved (PDF + PNG):
#    - _with_genes:  genome-wide significant loci annotated
#                    with gene names
#    - _no_genes:    unannotated (always produced)
#
# =============================================================

#### 00. Set-up ####

library(here)
library(data.table)
library(ggplot2)
library(dplyr)
library(openxlsx)
source(here("R", "manhattan_f.R"))

# Plot dimensions (inches)
FIG_W <- 8
FIG_H <- 4

# Output directory
OUT_DIR <- here("analysis", "figures", "manhattans")

# -------------------------------------------------------------
# Helper: save a plot as both PDF and PNG
# -------------------------------------------------------------
save_plot <- function(p, path_stem) {
  ggsave(paste0(path_stem, ".pdf"), plot = p,
         device = cairo_pdf, width = FIG_W, height = FIG_H)
  ggsave(paste0(path_stem, ".png"), plot = p,
         device = "png",       width = FIG_W, height = FIG_H)
}

# -------------------------------------------------------------
# Helper: shorten a "/" -separated gene list to at most 2 names
# -------------------------------------------------------------
shorten_genes <- function(gene_string) {
  genes <- strsplit(gene_string, "/")[[1]]
  if (length(genes) > 2) {
    paste0(paste(genes[1:2], collapse = "/"), "/...")
  } else {
    paste(genes, collapse = "/")
  }
}


#### 01. Individual Manhattan plots ####

# Short codes used in output file names
casec_codes <- c(
  "combined"     = "combi",
  "case_control" = "casec",
  "proxy"        = "proxy"
)

for (anc in c("all", "eur", "afr", "eas", "sas", "amr")) {
  for (casec in names(casec_codes)) {

    casec_out <- casec_codes[[casec]]
    stem      <- file.path(OUT_DIR, paste0("manhattan_main_", casec_out, "_", anc))

    # Skip if sumstats file does not exist
    sst_file <- here(
      "data", "sumstats", "meta", "main", casec,
      paste0("main_", casec, "_", anc, "_neff0.6_nsumstats1.txt.gz")
    )
    if (!file.exists(sst_file)) next

    # Determine which output files are missing
    has_genes_missing <- !file.exists(paste0(stem, "_with_genes.pdf")) |
                         !file.exists(paste0(stem, "_with_genes.png"))
    no_genes_missing  <- !file.exists(paste0(stem, "_no_genes.pdf"))   |
                         !file.exists(paste0(stem, "_no_genes.png"))

    if (!has_genes_missing & !no_genes_missing) {
      message("Skipping (all outputs exist): ", basename(stem))
      next
    }

    message("Loading sumstats: ", basename(sst_file))
    sst <- fread(sst_file)

    # Always produce the unannotated plot if missing
    if (no_genes_missing) {
      message("Plotting (no genes): ", basename(stem))
      save_plot(manhattan_f(sst, annotate_loci = FALSE), paste0(stem, "_no_genes"))
    }

    # Produce the annotated plot only when genome-wide significant hits exist
    if (has_genes_missing) {
      if (!any(sst$p_value <= 5e-8)) {
        message("No GWS hits — skipping annotated plot for: ", basename(stem))
        next
      }

      # Load gene annotations
      selected_genes <- read.xlsx(
          here("analysis", "gene_table", "manual_selection", "2025_06_23_gene_table.xlsx")
        ) %>%
        select(locus, manual_genes, novel) %>%
        mutate(gene = sapply(manual_genes, shorten_genes)) %>%
        select(-manual_genes)

      loci_info <- readRDS(here("analysis", "risk_loci", "risk_loci_table.rds"))[["unique_loci_main_combined"]] %>%
        select(locus, index_variant_id, ancestry)

      # For multi-ancestry loci, extract the index variant for this specific ancestry
      genes <- selected_genes %>%
        left_join(loci_info, by = "locus") %>%
        filter(grepl(anc, ancestry)) %>%
        rowwise() %>%
        mutate(variant_id = {
          ancs     <- stringr::str_split(ancestry,        ";")[[1]]
          variants <- stringr::str_split(index_variant_id, ";")[[1]]
          idx      <- which(ancs == anc)
          if (length(idx) > 0 && idx <= length(variants)) variants[[idx]] else NA_character_
        }) %>%
        ungroup()

      message("Plotting (with genes): ", basename(stem))
      save_plot(manhattan_f(sst, annotate_loci = TRUE, loci = genes), paste0(stem, "_with_genes"))
    }
  }
}
