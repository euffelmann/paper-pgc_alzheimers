# =============================================================
#
#  githubFigure_forest_plots.R
#  Generate per-locus forest plots for all 127 genome-wide
#  significant loci (48 novel + 79 known) from the PGC-ALZ3
#  meta-analysis.
#
#  Outputs are saved to forest_plots/novel_loci/ and
#  forest_plots/known_loci/ for upload to the GitHub repository
#  (not included in the manuscript supplementary information).
#
# =============================================================

#### 0. Set-up ####

library(data.table)
library(dplyr)
library(here)
library(meta)
library(grid)

base_dir <- "/home/emilu/projects/pgc_alzheimers/analysis/figures/forest_plots"

#### 1. Functions ####

extract_snps <- function(locus_nr, loci, analysis, ancestry, proxy_cc_combined, params, output_dir) {

  if (nrow(params) == 0) {
    stop("params has 0 rows for ancestry '", ancestry, "' — check your ancestry filter.")
  }

  variant_id <- loci$index_variant_id[loci$locus == locus_nr]
  plot_sst_dir <- file.path(output_dir, "sumstats")
  plot_sst_file <- paste0(gsub(pattern = ":", replacement = "_", x = variant_id), ".txt")

  if (!file.exists(file.path(plot_sst_dir, plot_sst_file))) {

    snps <- NULL
    for (j in seq_len(nrow(params))) {

      file <- params$sumstats_file[j]
      name <- params$name[j]; print(name)
      casec_proxy <- params$casec_proxy[j]
      ancestry_j <- params$ancestry[j]

      temp <- fread(cmd = paste0("gunzip -c ", file, " | awk 'NR==1 || /", variant_id, "/'"))
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
          cohort = paste(name, casec_proxy, ancestry_j, sep = "|")
        )

      } else {

        sumstats_colnames <- colnames(fread(cmd = paste0("gunzip -c ", params$sumstats_file[j], " | head -1")))
        names(temp) <- sumstats_colnames
        temp <- temp %>%
          mutate(
            cohort = paste(name, casec_proxy, ancestry_j, sep = "|"),
            maf = pmin(effect_allele_frequency, 1 - effect_allele_frequency)
          ) %>%
          select(
            variant_id, chromosome, base_pair_location, other_allele, effect_allele,
            beta, standard_error, p_value, effect_allele_frequency, maf, neff,
            n_case, n_control, cohort
          )
      }

      snps <- rbind(snps, temp)
    }

    return(snps)

  } else {
    print(paste0("SNPs have already been extracted and saved in: ", plot_sst_dir, "/", plot_sst_file))
    return(NULL)
  }
}

save_forest <- function(locus_nr, loci, analysis, ancestry, proxy_cc_combined, meta_sumstats, output_dir, locus_label = NULL) {

  variant_id <- loci$index_variant_id[loci$locus == locus_nr]
  index_variant_id <- gsub(pattern = ":", replacement = "_", variant_id)

  ## zero-padded locus number (SD2 label if supplied, else local locus_nr) so files sort 001, 002, ... 127
  locus_num <- if (!is.null(locus_label)) locus_label else locus_nr
  locus_str <- formatC(locus_num, width = 3, flag = "0")

  out_file <- file.path(output_dir, paste0("locus_", locus_str, "_", index_variant_id, ".pdf"))

  if (!file.exists(out_file)) {

    sumstats <- fread(file.path(output_dir, "sumstats", paste0(index_variant_id, ".txt")))

    temp <- meta_sumstats %>%
      filter(variant_id == !!variant_id)

    for (j in seq_len(nrow(sumstats))) {
      if ((sumstats$effect_allele[j] != temp$effect_allele) & !is.na(sumstats$effect_allele[j])) {
        sumstats$beta[j] <- sumstats$beta[j] * -1
        sumstats$effect_allele[j] <- temp$effect_allele
        sumstats$other_allele[j] <- temp$other_allele
        sumstats$effect_allele_frequency[j] <- 1 - sumstats$effect_allele_frequency[j]
      }
    }

    sumstats$casec_proxy <- sapply(strsplit(sumstats$cohort, "\\|"), `[`, 2)
    sumstats$casec_proxy <- factor(
      dplyr::recode(sumstats$casec_proxy, casec = "case_control", proxy = "proxy"),
      levels = c("proxy", "case_control")
    )

    m.gen <- meta::metagen(
      TE = beta,
      seTE = standard_error,
      pval = p_value,
      studlab = cohort,
      subgroup = casec_proxy,
      data = sumstats,
      common = TRUE,
      random = FALSE,
      title = sumstats$variant_id[1]
    )

    square_colors <- ifelse(m.gen$TE < 0, "#ca0020", "#0571b0")
    options(na.action = "na.pass")

    n_cohorts <- nrow(sumstats)
    n_subgroups <- length(unique(sumstats$casec_proxy))
    n_rows_est <- n_cohorts + n_subgroups * 2 + 8
    row_height_in <- 0.22
    title_margin_in <- 0.6
    bottom_margin_in <- 0.35
    pdf_height <- max(11, n_rows_est * row_height_in) + title_margin_in + bottom_margin_in

    pdf(out_file, width = 12, height = pdf_height, family = "Courier")

    grid::pushViewport(grid::viewport(
      y = grid::unit(bottom_margin_in, "inches"),
      height = grid::unit(1, "npc") - grid::unit(title_margin_in + bottom_margin_in, "inches"),
      just = "bottom"
    ))
    meta::forest(
      m.gen,
      prediction = FALSE,
      print.tau2 = FALSE,
      fontsize = 9,
      leftcols = c("studlab", "TE", "seTE", "pval"),
      leftlabs = c("Cohort", "log(OR)", "Standard Error", "P-value"),
      xlim = c(min(min(sumstats$beta, na.rm = T) - 0.1, -0.1),
               max(max(sumstats$beta, na.rm = T) + 0.1, 0.1)),
      digits = 4,
      scientific.pval = TRUE,
      lab.NA = "NA",
      ref = 0,
      col.square = square_colors,
      col.square.lines = square_colors,
      col.diamond = "black",
      text.common = "Combined common effect model",
      text.common.w = c("Proxy common effect model", "Case-control common effect model"),
      print.subgroup.name = FALSE,
      new = FALSE
    )
    grid::popViewport()

    title_text <- if (!is.null(locus_label)) {
      paste0("Locus ", locus_label, "   |   ", sumstats$variant_id[1])
    } else {
      sumstats$variant_id[1]
    }

    grid::grid.text(
      title_text,
      x = 0.5,
      y = grid::unit(1, "npc") - grid::unit(0.2, "inches"),
      gp = grid::gpar(fontsize = 16, fontface = "bold")
    )
    grid::grid.text(
      paste0("effect allele:", temp$effect_allele, " frequency:", round(temp$effect_allele_frequency, 4),
             " neff:", format(temp$neff, big.mark = ",")),
      x = 0.5,
      y = grid::unit(1, "npc") - grid::unit(0.45, "inches"),
      gp = grid::gpar(fontsize = 12, fontface = "bold")
    )

    dev.off()
  }
}

#### 2. Locus lists ####

analysis <- "main"
proxy_cc_combined <- "combined"

novel_locus_nrs_sd2 <- c(1, 2, 3, 6, 15, 16, 17, 20, 21, 25, 28, 32, 37, 39, 40, 41, 42, 45, 50,
                         55, 57, 59, 60, 61, 62, 68, 72, 75, 76, 78, 79, 80, 85, 86, 90, 92, 94,
                         98, 101, 104, 113, 114, 115, 116, 117, 118, 124, 125)

novel_lead_snps <- c(
  "1:21134597:A:G", "1:31225678:A:T", "1:150149016:A:G", "1:201012937:C:T",
  "3:56273512:A:C", "3:136626697:G:T", "3:152177852:A:G", "3:185155925:A:T",
  "3:196231308:A:G", "4:159857539:A:G", "5:122705110:A:G", "5:176982662:C:T",
  "6:45407654:A:G", "6:71097590:C:T", "6:153371503:C:T", "7:1575972:A:G",
  "7:6069653:C:G", "7:17911038:C:T", "7:111580166:C:T", "8:103577024:C:G",
  "8:134064152:C:T", "9:21416492:C:T", "9:80718995:A:G", "9:93563884:C:T",
  "9:95626730:A:C", "10:113463315:A:G", "11:65599656:C:T", "12:901559:A:G",
  "12:56865338:A:G", "12:123349410:C:T", "12:124167021:C:T", "13:113672477:C:T",
  "15:38850330:C:T", "15:40589218:C:G", "15:69573077:G:T", "15:94802545:A:C",
  "16:23560398:A:G", "16:72207315:C:T", "16:87225431:A:G", "17:2302387:C:T",
  "19:5908969:A:G", "19:10837677:A:G", "19:18513594:C:T", "19:19575965:A:G",
  "19:33370952:A:G", "19:41734958:G:T", "20:4142809:A:G", "20:48578734:A:T"
)

known_locus_nrs_sd2 <- c(4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 18, 19, 22, 23, 24, 26, 27, 29, 30,
                         31, 33, 34, 35, 36, 38, 43, 44, 46, 47, 48, 49, 51, 52, 53, 54, 56, 58,
                         63, 64, 65, 66, 67, 69, 70, 71, 73, 74, 77, 81, 82, 83, 84, 87, 88, 89,
                         91, 93, 95, 96, 97, 99, 100, 102, 103, 105, 106, 107, 108, 109, 110,
                         111, 112, 119, 120, 121, 122, 123, 126, 127)

known_lead_snps <- c(
  "1:155135036:A:G", "1:161187665:C:T", "1:207684192:G:T", "2:37476688:C:G",
  "2:65608363:A:G", "2:106385448:A:G", "2:127891427:A:C", "2:135609247:A:G",
  "2:203719028:A:G", "2:233981912:C:G", "3:154746310:C:T", "3:183995854:A:G",
  "4:1006200:A:G", "4:11014822:A:G", "4:40198842:A:G", "5:14678996:A:G",
  "5:86223195:C:T", "5:139720400:C:T", "5:150432388:C:T", "5:156526331:A:G",
  "5:179624032:A:G", "6:28348158:G:T", "6:32589158:A:G", "6:41129207:C:T",
  "6:47483653:C:G", "7:7855379:C:T", "7:12264737:C:T", "7:28156887:C:T",
  "7:37698841:A:T", "7:54941328:C:T", "7:99971834:A:G", "7:143110762:A:G",
  "8:11702122:C:G", "8:27464929:A:G", "8:95969322:C:T", "8:126580618:C:T",
  "8:145158607:A:G", "9:107661129:A:C", "10:11720308:A:G", "10:61581307:A:C",
  "10:82273079:C:T", "10:98061083:C:T", "10:124179299:C:T", "11:47380340:G:T",
  "11:60021948:A:G", "11:85867875:A:G", "11:121435587:C:T", "12:113731926:G:T",
  "14:53414145:C:T", "14:92938855:A:G", "14:106194461:A:G", "14:107100641:A:G",
  "15:59000957:C:T", "15:63569902:C:T", "15:64351270:A:T", "15:79234957:A:G",
  "16:19740237:A:G", "16:30068354:C:T", "16:31134059:A:G", "16:70694000:A:C",
  "16:79608831:C:G", "16:81942028:C:G", "16:90160829:A:G", "17:1623675:A:G",
  "17:5133128:C:T", "17:17747514:C:T", "17:42430244:C:T", "17:44101563:C:T",
  "17:47407071:A:C", "17:56404349:A:G", "17:61545779:C:T", "19:1050874:A:G",
  "19:45411941:C:T", "19:49229323:A:G", "19:51727962:A:C", "19:54814234:C:T",
  "20:396156:A:G", "20:54995699:C:T", "21:27481343:A:C"
)

stopifnot(length(novel_locus_nrs_sd2) == 48, length(novel_lead_snps) == 48)
stopifnot(length(known_locus_nrs_sd2) == 79, length(known_lead_snps) == 79)

## config for the two categories — this drives everything downstream
categories <- list(
  novel = list(
    output_dir  = file.path(base_dir, "novel_loci"),
    locus_nrs   = novel_locus_nrs_sd2,
    lead_snps   = novel_lead_snps,
    label       = "novel"
  ),
  known = list(
    output_dir  = file.path(base_dir, "known_loci"),
    locus_nrs   = known_locus_nrs_sd2,
    lead_snps   = known_lead_snps,
    label       = "known"
  )
)

#### 3. Build params ####

params_cc <- data.frame(
  sumstats_file = list.files(here("data", "sumstats", "case_control"), pattern = "\\.gz$", full.names = T),
  ancestry = sub(".*_(eur|eas|afr|amr|sas).*", "\\1",
                 list.files(here("data", "sumstats", "case_control"), pattern = "\\.gz$"))
) %>%
  filter(
    grepl("main", sumstats_file) |
      (grepl("noag", sumstats_file) &
         (grepl("mvplo", sumstats_file) | grepl("biovu", sumstats_file) | grepl("grace", sumstats_file) |
            grepl("xstsa", sumstats_file) | grepl("twing", sumstats_file) | grepl("gothe", sumstats_file)))
  )

params_proxy <- data.frame(
  sumstats_file = list.files(here("data", "sumstats", "proxy_meta"), pattern = "\\.gz$", full.names = T),
  ancestry = sub(".*_(eur|eas|afr|amr|sas).*", "\\1",
                 list.files(here("data", "sumstats", "proxy_meta"), pattern = "\\.gz$"))
) %>%
  filter(grepl("main", sumstats_file) | grepl("mvplo", sumstats_file), !grepl("inref", sumstats_file))

params <- rbind(params_cc, params_proxy) %>%
  mutate(
    name = sub(".*clean_([a-zA-Z0-9]{5}).*", "\\1", sumstats_file),
    casec_proxy = ifelse(grepl("case_control", sumstats_file), "casec", "proxy")
  ) %>%
  arrange(casec_proxy, factor(ancestry, levels = c("eur", "afr", "eas", "amr", "sas")), name)

temp_params_all <- params
temp_params_eur <- params %>% filter(ancestry == "eur")

#### 4. Load loci and meta sumstats ####

loci_all <- fread(here("analysis/risk_loci", analysis,
                       paste(analysis, proxy_cc_combined, "all", "neff0.6_nsumstats1_loci.txt", sep = "_")))
loci_eur <- fread(here("analysis/risk_loci", analysis,
                       paste(analysis, proxy_cc_combined, "eur", "neff0.6_nsumstats1_loci.txt", sep = "_")))

meta_sumstats_all <- fread(here(
  "data/sumstats/meta/", analysis, proxy_cc_combined,
  paste0(analysis, "_", proxy_cc_combined, "_all_neff0.6_nsumstats1.txt.gz")
))
meta_sumstats_eur <- fread(here(
  "data/sumstats/meta/", analysis, proxy_cc_combined,
  paste0(analysis, "_", proxy_cc_combined, "_eur_neff0.6_nsumstats1.txt.gz")
))

#### 5. Extract sumstats and save forest plots ####

for (cat_name in names(categories)) {

  cat <- categories[[cat_name]]
  dir.create(file.path(cat$output_dir, "sumstats"), recursive = TRUE, showWarnings = FALSE)

  message("\n==== Processing ", cat_name, " loci (", length(cat$locus_nrs), " total) ====")

  lookup <- data.frame(
    sd2_locus_nr = cat$locus_nrs,
    lead_snp     = cat$lead_snps
  )
  lookup$local_locus_nr_all <- loci_all$locus[match(lookup$lead_snp, loci_all$index_variant_id)]
  lookup$local_locus_nr_eur <- loci_eur$locus[match(lookup$lead_snp, loci_eur$index_variant_id)]

  lookup$stratum <- ifelse(!is.na(lookup$local_locus_nr_all), "all",
                           ifelse(!is.na(lookup$local_locus_nr_eur), "eur", NA))
  lookup$local_locus_nr <- ifelse(lookup$stratum == "all",
                                  lookup$local_locus_nr_all,
                                  lookup$local_locus_nr_eur)

  unresolved <- lookup %>% filter(is.na(stratum))
  if (nrow(unresolved) > 0) {
    print(unresolved)
    stop(nrow(unresolved), " ", cat_name, " lead SNP(s) not found in either the 'all' or 'eur' loci tables — needs manual investigation.")
  }

  message("All ", nrow(lookup), " ", cat_name, " loci resolved (", sum(lookup$stratum == "all"), " from 'all', ",
          sum(lookup$stratum == "eur"), " from 'eur'):")
  print(lookup[, c("sd2_locus_nr", "lead_snp", "stratum", "local_locus_nr")])

  for (i in seq_len(nrow(lookup))) {

    stratum  <- lookup$stratum[i]
    locus_nr <- lookup$local_locus_nr[i]
    sd2_nr   <- lookup$sd2_locus_nr[i]

    loci_use          <- if (stratum == "all") loci_all else loci_eur
    params_use        <- if (stratum == "all") temp_params_all else temp_params_eur
    meta_sumstats_use <- if (stratum == "all") meta_sumstats_all else meta_sumstats_eur

    variant_id <- loci_use$index_variant_id[loci_use$locus == locus_nr]
    message("[", cat_name, "/", stratum, "] SD2 locus ", sd2_nr, " / local locus ", locus_nr, " (", variant_id, ") — extracting sumstats...")

    snps <- tryCatch(
      extract_snps(
        locus_nr = locus_nr, loci = loci_use, analysis = analysis,
        ancestry = stratum, proxy_cc_combined = proxy_cc_combined, params = params_use,
        output_dir = cat$output_dir
      ),
      error = function(e) { message("  ERROR extracting locus ", locus_nr, ": ", e$message); NULL }
    )

    if (!is.null(snps)) {
      fwrite(
        snps,
        file = file.path(cat$output_dir, "sumstats", paste0(gsub(":", "_", variant_id), ".txt")),
        quote = FALSE, sep = "\t", na = NA
      )
    }

    message("[", cat_name, "/", stratum, "] SD2 locus ", sd2_nr, " — plotting...")
    tryCatch(
      save_forest(
        locus_nr = locus_nr, loci = loci_use, analysis = analysis,
        ancestry = stratum, proxy_cc_combined = proxy_cc_combined, meta_sumstats = meta_sumstats_use,
        output_dir = cat$output_dir,
        locus_label = sd2_nr
      ),
      error = function(e) message("  ERROR plotting locus ", locus_nr, ": ", e$message)
    )
  }

  message("Done with ", cat_name, ". Check: ", cat$output_dir)
}

message("\nAll done. novel_loci -> ", categories$novel$output_dir, " | known_loci -> ", categories$known$output_dir)
