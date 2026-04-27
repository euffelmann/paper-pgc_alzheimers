# =============================================================
#
#  cohort_summary.R
#  Generate comprehensive sample size summary tables for all
#  cohorts contributing to the PGC-ALZ3 meta-analysis.
#
#  This script:
#    1. Scans summary statistic files to extract sample sizes
#    2. Calculates effective sample size (neff) for each cohort
#    3. Organizes cohorts by analysis type and ancestry
#    4. Creates Excel workbook with separate sheets per analysis
#
#  Output: cohort_summary.xlsx with sheets for:
#    - main: Primary analysis (mixed age criteria)
#    - fema: Female-specific analysis
#    - male: Male-specific analysis
#    - noag: Analysis without age restrictions
#    - apoe: APOE-stratified analysis
#
# =============================================================

#### 1. Set-up ####
library(data.table)
library(here)
library(dplyr)
library(stringr)
library(openxlsx)
source(here("R", "neff_f.R"))  # Function to calculate effective sample size


#### 2. Process case-control cohorts ####

# Create data frame with metadata extracted from filenames
all_cc <- tibble(
  # List all case-control summary statistic files
  sumstats_file = list.files(
    here("data_raw", "sumstats", "case_control"),
    pattern = "\\.txt$", full.names = F),

  # Extract phenotype type from filename
  proxy_casecontrol = sub(
    ".*_(proxy|case_control).*", "\\1",
    list.files(here("data_raw", "sumstats", "case_control"), pattern = "\\.txt$")),

  # Extract ancestry from filename
  ancestry = sub(
    ".*_(eur|eas|afr|amr|sas).*", "\\1",
    list.files(here("data_raw", "sumstats", "case_control"), pattern = "\\.txt$")),

  # Extract analysis type from filename
  analysis = sub(
    ".*_(main|fema|male|noag|apoe).*", "\\1",
    list.files(here("data_raw", "sumstats", "case_control"), pattern = "\\.txt$")),

  # Initialize sample size columns
  max_n_case = NA,
  max_n_control = NA,
  neff = NA
) %>%
  mutate(
    # Extract cohort abbreviation from filename
    cohort = stringr::str_extract(sumstats_file, ".*(?=_case_control|_proxy)"),
    max_n_case = NA,
    max_n_control = NA,
    neff = NA
    )

## Extract case and control counts from file headers
# Loop through each file and extract sample sizes from header
for (i in 1:nrow(all_cc)) {
  # Extract case count from "- caseCount: N" line in file header
  all_cc$max_n_case[i] <- system(
    paste0("grep 'caseCount' ", here("data_raw/sumstats/case_control/"), all_cc$sumstats_file[i]),
    intern = T) %>%
    gsub("- caseCount: ", "", .) %>%
    as.integer()

  # Extract control count from "- controlCount: N" line in file header
  all_cc$max_n_control[i] <- system(
    command = paste0("grep 'controlCount' ", here("data_raw/sumstats/case_control/"), all_cc$sumstats_file[i]),
    intern = T) %>%
    gsub("- controlCount: ", "", .) %>%
    as.integer()
}

# Calculate effective sample size for case-control cohorts
# neff = 4 / (1/n_case + 1/n_control)
all_cc <- all_cc %>%
  mutate(
    neff = neff_f(n_case = max_n_case, n_control = max_n_control)
  )


#### 3. Process proxy phenotype cohorts ####

# Create data frame with metadata extracted from filenames
all_proxy <- data.frame(
  # List all proxy phenotype summary statistic files
  sumstats_file = list.files(
    here("data_raw", "sumstats", "proxy"),
    pattern = "\\.txt$", full.names = F),

  # Extract phenotype type from filename (all are "proxy")
  proxy_casecontrol = sub(
    ".*_(proxy|proxy).*", "\\1",
    list.files(here("data_raw", "sumstats", "proxy"), pattern = "\\.txt$")),

  # Extract ancestry from filename
  ancestry = sub(
    ".*_(eur|eas|afr|amr|sas).*", "\\1",
    list.files(here("data_raw", "sumstats", "proxy"), pattern = "\\.txt$")),

  # Extract analysis type from filename
  analysis = sub(
    ".*_(main|fema|male|noag|apoe).*", "\\1",
    list.files(here("data_raw", "sumstats", "proxy"), pattern = "\\.txt$")),

  # Initialize columns
  parent = NA,  # Will be "father" or "mother" for proxy phenotypes
  max_n_case = NA,
  max_n_control = NA,
  neff = NA
) %>%
  mutate(
    # Extract cohort abbreviation
    cohort = stringr::str_extract(sumstats_file, ".*(?=_proxy|_proxy)"),

    # Extract parent type (father/mother) for proxy phenotypes
    parent = ifelse(proxy_casecontrol == "proxy",
                    sub(".*_(father|mother).*", "\\1",
                        list.files(here("data_raw", "sumstats", "proxy"), pattern = "\\.txt$")),
                    NA),
    max_n_case = NA,
    max_n_control = NA,
    neff = NA
  )

## Extract case and control counts from file headers
# Loop through each file and extract sample sizes
for (i in 1:nrow(all_proxy)) {
  # Extract case count (individuals with affected parent)
  all_proxy$max_n_case[i] <- system(
    paste0("grep 'caseCount' ", here("data_raw/sumstats/proxy/"), all_proxy$sumstats_file[i]),
    intern = T) %>%
    gsub("- caseCount: ", "", .) %>%
    as.integer()

  # Extract control count (individuals with unaffected parent)
  all_proxy$max_n_control[i] <- system(
    command = paste0("grep 'controlCount' ", here("data_raw/sumstats/proxy/"), all_proxy$sumstats_file[i]),
    intern = T) %>%
    gsub("- controlCount: ", "", .) %>%
    as.integer()
}

# Calculate effective sample size for proxy phenotypes
# Proxy phenotypes have lower neff due to weaker genetic signal
all_proxy <- all_proxy %>%
  mutate(neff = neff_f(
    n_case = max_n_case,
    n_control = max_n_control,
    proxy_casecontrol = "proxy"  # Applies proxy-specific neff calculation
  ))

## Combine mother and father proxy cohorts
# For cohorts with separate mother/father analyses, combine them

# Vectorize neff_f for use with mutate
neff_f_vectorized <- Vectorize(neff_f)

# Combine mother and father data by summing counts
combined_proxy <- all_proxy %>%
  # Remove parent identifier from filename to group mother + father
  mutate(sumstats_file = gsub("_father_|_mother_", "_", sumstats_file)) %>%
  group_by(sumstats_file) %>%
  summarise(
    proxy_casecontrol = first(proxy_casecontrol),
    ancestry = first(ancestry),
    analysis = first(analysis),
    # Sum cases and controls across mother and father cohorts
    max_n_case = sum(max_n_case, na.rm = TRUE),
    max_n_control = sum(max_n_control, na.rm = TRUE),
    neff = sum(neff, na.rm = TRUE),
    cohort = first(cohort)
  ) %>%
  ungroup() %>%

  ## Special handling for MVP-Logue cohort
  # MVP-Logue used the same control set for both mother and father analyses,
  # so we divide by 2 to avoid double-counting controls
  mutate(max_n_control = case_when(
    cohort == "mvplo" ~ max_n_control / 2,
    .default = max_n_control)) %>%

  # Recalculate neff after combining cohorts
  mutate(
    neff = neff_f_vectorized(
      n_case = max_n_case,
      n_control = max_n_control,
      proxy_casecontrol = proxy_casecontrol
    )
  )

#### 4. Create summary tables per analysis type ####
# Load cohort name abbreviations lookup table
abbr <- as_tibble(fread(here("analysis/cohort_summary/abbr.txt")))

# Initialize list to store tables for each analysis
df_list <- list()

## Main analysis table
# Mix of "main" and "noag" (no age restriction) to maximize sample size
# Some cohorts didn't have age data, so we use their "noag" version
df_list[["main"]] <- rbind(all_cc, combined_proxy) %>%
  mutate(ancestry = factor(ancestry, levels = c("eur", "afr", "eas", "sas", "amr"))) %>%
  filter(
    analysis == "main" |
    # Include no-age cohorts that didn't have age data available
    (analysis == "noag" & cohort %in% c("mvplo", "biovu", "grace", "xstsa", "twing", "gothe"))
    ) %>%
  arrange(proxy_casecontrol, ancestry, cohort) %>%
  select(-sumstats_file) %>%
  left_join(abbr, by = c("cohort" = "abbr")) %>%
  rename(cohort_long = name) %>%
  relocate(cohort)

## Female-specific analysis table
# For proxy cohorts, use mother proxy phenotype (affected mother)
fema_proxy <- all_proxy %>%
  filter(parent == "mother" & analysis == "main") %>%
  mutate(analysis = parent) %>%
  select(-parent)

df_list[["fema"]] <- all_cc %>%
  filter(analysis == "fema") %>%
  rbind(fema_proxy) %>%
  mutate(ancestry = factor(ancestry, levels = c("eur", "afr", "eas", "sas", "amr"))) %>%
  arrange(proxy_casecontrol, ancestry, cohort) %>%
  select(-sumstats_file) %>%
  left_join(abbr, by = c("cohort" = "abbr")) %>%
  rename(cohort_long = name) %>%
  relocate(cohort)

## Male-specific analysis table
# For proxy cohorts, use father proxy phenotype (affected father)
male_proxy <- all_proxy %>%
  filter(parent == "father" & analysis == "main") %>%
  mutate(analysis = parent) %>%
  select(-parent)

df_list[["male"]] <- all_cc %>%
  filter(analysis == "male") %>%
  rbind(male_proxy) %>%
  mutate(ancestry = factor(ancestry, levels = c("eur", "afr", "eas", "sas", "amr"))) %>%
  arrange(proxy_casecontrol, ancestry, cohort) %>%
  select(-sumstats_file) %>%
  left_join(abbr, by = c("cohort" = "abbr")) %>%
  rename(cohort_long = name) %>%
  relocate(cohort)

## No-age restriction analysis table
# Analysis including all individuals regardless of age at diagnosis
df_list[["noag"]] <- rbind(all_cc, combined_proxy) %>%
  mutate(ancestry = factor(ancestry, levels = c("eur", "afr", "eas", "sas", "amr"))) %>%
  filter(analysis == "noag") %>%
  arrange(proxy_casecontrol, ancestry, cohort) %>%
  select(-sumstats_file) %>%
  left_join(abbr, by = c("cohort" = "abbr")) %>%
  rename(cohort_long = name) %>%
  relocate(cohort)

## APOE-stratified analysis table
# Analysis stratified by APOE genotype
df_list[["apoe"]] <- rbind(all_cc, combined_proxy) %>%
  mutate(ancestry = factor(ancestry, levels = c("eur", "afr", "eas", "sas", "amr"))) %>%
  filter(analysis == "apoe") %>%
  arrange(proxy_casecontrol, ancestry, cohort) %>%
  select(-sumstats_file) %>%
  left_join(abbr, by = c("cohort" = "abbr")) %>%
  rename(cohort_long = name) %>%
  relocate(cohort)


#### 5. Save Excel workbook with all analysis tables ####

print("Saving Excel file with cohort summary tables...")

# Write list of data frames to separate Excel sheets
# Each sheet corresponds to one analysis type (main, fema, male, noag, apoe)
write.xlsx(df_list, file = here("analysis/cohort_summary/cohort_summary.xlsx"))

print("Done! Output saved to: analysis/cohort_summary/cohort_summary.xlsx")
