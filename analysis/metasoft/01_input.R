#### 0. set-up ####
library(dplyr)
library(here)
library(data.table)
library(stringr)
library(purrr)


#### 1. make input ####
for (proxy_casec in c("case_control", "combined", "proxy")) {
  
  # initiate empty list
  sumstats <- NULL
  
  # List files in directory
  files <- list.files(here(paste0("data/sumstats/meta/main/", proxy_casec)))
  
  # Extract ancestry abbreviations (excluding "all")
  ancestries <- files %>%
    str_extract(paste0("(?<=main_", proxy_casec, "_)[a-z]+")) %>%
    unique() %>%
    .[. != "all"] %>%
    .[!is.na(.)]
  
  print("present ancestries: ")
  print(ancestries)
  
  for (anc in ancestries) {
    
    print("preparing input for: ")
    print(paste0(anc, " & ", proxy_casec))
    
    sumstats[[anc]] <- fread(paste0(here("data/sumstats/meta/main", proxy_casec), "/main_", proxy_casec, "_", anc, "_neff0.6_nsumstats1.txt.gz")) %>%
      select(variant_id, beta, standard_error)
    
  }
  
  combined <- sumstats |>
    imap(~ rename(.x, "{.y}_beta" := beta, "{.y}_se" := standard_error)) |>
    reduce(full_join, by = "variant_id")
  
  
  print("saving input dataframe")
  fwrite(
    x = combined, 
    file = paste0(here("analysis/metasoft/input/"), proxy_casec, "_neff0.6_nsumstats1.txt"), 
    quote = F, 
    col.names = F, 
    sep = " ", 
    na = "NA")
  
}





