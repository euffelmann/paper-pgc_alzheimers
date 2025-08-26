#### 0. setup ####
library(here)
library(data.table)
library(ggplot2)
library(dplyr)
library(openxlsx)
source(here("R/manhattan_f.R"))

#### 1. manhattan plot ####
casec <- "combined"
anc <- "all"

if (casec == "combined") {casec_out <- "combi"}
if (casec == "case_control") {casec_out <- "casec"}
if (casec == "proxy") {casec_out <- "proxy"}

print(paste0("loading sumstats: main_", casec, "_", anc, "_neff0.6_nsumstats1.txt.gz"))

sst <- fread(here("data/sumstats/meta/main", casec, paste0("main_", casec, "_", anc, "_neff0.6_nsumstats1.txt.gz"))) %>%
  mutate(p_value = as.numeric(p_value)) %>%
  filter(p_value <= 5e-5)


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

p <- manhattan_f(sst, annotate_loci = T, loci = genes)

ggsave(
  plot = p,
  filename = here("analysis/figures/manhattans", paste0("manhattan_main_", casec_out, "_", anc, "_with_genes.svg")),
  device = "svg",
  width = 8,
  height = 4
)
