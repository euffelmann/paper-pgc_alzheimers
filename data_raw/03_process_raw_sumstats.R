# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                         #
#         script to process raw summary statistics        #
#                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# - - - - - Table of Contents - - - - - - - - - - - - - - #
# - - - - - 00. set-up - - - - - - - - - - - - - - - - -  #
# - - - - - 01. QC plots of raw data - - - - - - - - - -  #
# - - - - - 02. Perform QC of raw data - - - - - - - - -  #
# - - - - - 03. QC plots of cleaned data - - - - - - - -  #

#### 00. set-up ####

## load packages
library(data.table)
library(dplyr)
library(purrr)
library(here)
library(qqman)
library(stringr)
source(here("R", "qc_plots.R"))
source(here("R", "process_sumstats.R"))
source(here("R", "flip.R"))

cohort <- as.character(commandArgs(T)[1])

#
#
#
#

#### 01. QC plots of raw data ####

params <- data.frame(
  sumstats_file = c(
    list.files(here("data_raw", "sumstats", "case_control"), pattern = "\\.gz$", full.names = T),
    list.files(here("data_raw", "sumstats", "proxy"), pattern = "\\.gz$", full.names = T),
    list.files(here("data_raw", "sumstats", "combined"), pattern = "\\.gz$", full.names = T)
  ),
  outfile = c(
    sub("\\.txt\\.gz$", "", list.files(here("data_raw", "sumstats", "case_control"), pattern = "\\.gz$")),
    sub("\\.txt\\.gz$", "", list.files(here("data_raw", "sumstats", "proxy"), pattern = "\\.gz$")),
    sub("\\.txt\\.gz$", "", list.files(here("data_raw", "sumstats", "combined"), pattern = "\\.gz$"))
  ),
  qc_dir = here("data_raw", "sumstats", "qc_plots"),
  ancestry = c(
    sub(".*_(eur|eas|afr|amr|sas).*", "\\1", list.files(here("data_raw", "sumstats", "case_control"), pattern = "\\.gz$")),
    sub(".*_(eur|eas|afr|amr|sas).*", "\\1", list.files(here("data_raw", "sumstats", "proxy"), pattern = "\\.gz$")),
    sub(".*_(eur|eas|afr|amr|sas).*", "\\1", list.files(here("data_raw", "sumstats", "combined"), pattern = "\\.gz$"))
  )) %>%
  mutate(
    cohort = str_extract(outfile, ".*(?=_case_control|_proxy|_combined)")
  ) %>%
  filter(cohort == !!cohort)

params %>%
  pmap(~qc_plots(
    sumstats_file = ..1,
    outfile = ..2,
    qc_dir = ..3,
    ancestry = ..4
  ))

#
#
#
#

#### 02. Perform QC of raw data #### 

## make dataframe of parameters
params <- data.frame(
  sumstats_file = c(
    list.files(here("data_raw", "sumstats", "case_control"), pattern = "\\.gz$", full.names = T),
    list.files(here("data_raw", "sumstats", "proxy"), pattern = "\\.gz$", full.names = T),
    list.files(here("data_raw", "sumstats", "combined"), pattern = "\\.gz$", full.names = T)
  ),
  outfile = c(
    sub("\\.txt\\.gz$", "", list.files(here("data_raw", "sumstats", "case_control"), pattern = "\\.gz$")),
    sub("\\.txt\\.gz$", "", list.files(here("data_raw", "sumstats", "proxy"), pattern = "\\.gz$")),
    sub("\\.txt\\.gz$", "", list.files(here("data_raw", "sumstats", "combined"), pattern = "\\.gz$"))
  ),
  ancestry = c(
    sub(".*_(eur|eas|afr|amr|sas).*", "\\1", list.files(here("data_raw", "sumstats", "case_control"), pattern = "\\.gz$")),
    sub(".*_(eur|eas|afr|amr|sas).*", "\\1", list.files(here("data_raw", "sumstats", "proxy"), pattern = "\\.gz$")),
    sub(".*_(eur|eas|afr|amr|sas).*", "\\1", list.files(here("data_raw", "sumstats", "combined"), pattern = "\\.gz$"))
    )
  ) %>%
  mutate(
    proxy_casecontrol = ifelse(grepl("proxy", sumstats_file), "proxy", "case_control"),
    proxy_casecontrol = ifelse(grepl("combined", sumstats_file), "combined", proxy_casecontrol),
    outdir = here("data", "sumstats", proxy_casecontrol),
    cohort = str_extract(outfile, ".*(?=_case_control|_proxy|_combined)")
  ) %>%
  filter(cohort == !!cohort) 


# add whether dosage or hard calls were used
for (i in 1:nrow(params)) {
  params$dosage_or_hardcalls[i] <- system(paste0(
    "grep 'dosageOrHardcalls:' ",
    here(
      "data_raw",
      "sumstats",
      params$proxy_casecontrol[i],
      paste0(params$outfile[i], "_readme.txt")
    )
  ), intern = T) %>% str_remove("- dosageOrHardcalls: ")
}

## create qc_summary file
if (!file.exists(here("data", "sumstats", "qc_summary", paste0("qc_summary_", cohort, ".txt")))) {
  data.frame(file = as.character()) %>%
    fwrite(
      file = here("data", "sumstats", "qc_summary", paste0("qc_summary_", cohort, ".txt")),
      quote = F,
      sep = "\t",
      na = "NA",
      row.names = F
    )
}

## run processing function
params %>%
  pmap(~ process_sumstats(
    sumstats_file = ..1,
    outfile = ..2,
    ancestry = ..3, 
    proxy_casecontrol = ..4, 
    outdir = ..5,
    dosage_or_hardcalls = ..6,
    maf_filter = 0.001
  ))


#
#
#
#

#### 03. QC plots of cleaned data ####

params <- data.frame(
  sumstats_file = c(
    list.files(here("data", "sumstats", "case_control"), pattern = "\\.gz$", full.names = T),
    list.files(here("data", "sumstats", "proxy"), pattern = "\\.gz$", full.names = T),
    list.files(here("data", "sumstats", "combined"), pattern = "\\.gz$", full.names = T)
  ),
  outfile = c(
    sub("\\.txt\\.gz$", "", list.files(here("data", "sumstats", "case_control"), pattern = "\\.gz$")),
    sub("\\.txt\\.gz$", "", list.files(here("data", "sumstats", "proxy"), pattern = "\\.gz$")),
    sub("\\.txt\\.gz$", "", list.files(here("data", "sumstats", "combined"), pattern = "\\.gz$"))
  ),
  qc_dir = here("data", "sumstats", "qc_plots"),
  ancestry = c(
    sub(".*_(eur|eas|afr|amr|sas).*", "\\1", list.files(here("data", "sumstats", "case_control"), pattern = "\\.gz$")),
    sub(".*_(eur|eas|afr|amr|sas).*", "\\1", list.files(here("data", "sumstats", "proxy"), pattern = "\\.gz$")),
    sub(".*_(eur|eas|afr|amr|sas).*", "\\1", list.files(here("data", "sumstats", "combined"), pattern = "\\.gz$"))
  )) %>%
  mutate(
    cohort = str_extract(outfile, ".*(?=_case_control|_proxy|_combined)")
  ) %>%
  filter(cohort == paste0("clean_", !!cohort))

params %>%
  pmap(~qc_plots(
    sumstats_file = ..1,
    outfile = ..2,
    qc_dir = ..3,
    ancestry = ..4
  ))

#
#
#
#