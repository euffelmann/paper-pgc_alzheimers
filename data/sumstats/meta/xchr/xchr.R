#### 0. set-up ####
library(data.table)
library(dplyr)
library(here)


#### 1. load and filter main sumstats ####
for (casec in c("case_control", "proxy", "combined")) {
  for (anc in c("eur", "afr", "eas", "sas", "amr", "all")) {
    
    print(paste0("xchr_", casec, "_", anc))
    
    x_sumstats <- fread(here(
      "data/sumstats/meta/main",
      casec,
      paste0("main_", casec, "_", anc, ".txt.gz")
    )) %>%
      filter(chromosome == 23) %>%
      filter(neff >= 0.6 * max(.$neff))
    
    if (nrow(x_sumstats) > 0) {
      
      fwrite(
        x = x_sumstats,
        file = here(
          "data/sumstats/meta/xchr",
          casec,
          paste0("xchr_", casec, "_", anc, "_neff0.6_nsumstats1.txt.gz")
        ),
        append = F,
        quote = F,
        sep = "\t",
        na = NA
      )
      
    } else {
      
      print(paste0("no X-chromosome SNPs for ", "xchr_", casec, "_", anc))
      
    }
  }
}