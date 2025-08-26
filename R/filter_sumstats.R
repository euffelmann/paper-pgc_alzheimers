#### 00. set-up ####
if (!require(dplyr)) { install.packages("dplyr"); library(dplyr) }
if (!require(data.table)) { install.packages("data.table"); library(data.table) }  

### 01. function to filter sumstats ####
filter_sumstats <- function(input, neff_th = 0.6, n_sumstats_th = 1) {
  ## assumes GWAS catalog sumstats format
  
  sumstats <- fread(paste0(input, ".txt.gz")) %>%
    mutate(n_sumstats = nchar(direction) - stringr::str_count(direction, "\\?"))
  
  print(paste0("... filtering file: ", input, " ..."))
  
  output <- paste0(input, "_neff", neff_th, "_nsumstats", n_sumstats_th)
  
  sumstats %>%
    filter(neff >= neff_th * max(.$neff), n_sumstats >= n_sumstats_th) %>%
    fwrite(
      file = paste0(output, ".txt.gz"),
      quote = F,
      sep = "\t",
      na = "NA"
    )
}

#### 02. filter sumstats ####
input <- as.character(commandArgs(T)[1])
neff_th <- as.numeric(commandArgs(T)[2])
n_sumstats_th <- as.integer(commandArgs(T)[3])

filter_sumstats(
  input = input,
  neff_th = neff_th,
  n_sumstats_th = n_sumstats_th
)