#### 0. set-up ####
source("meta_analysis_f.R")
source("run_metal.R")

#### 1. load input ####
sumstats <- unlist(strsplit(as.character(commandArgs(T)[1]), ","))
out_dir <- as.character(commandArgs(T)[2])
outfile <- as.character(commandArgs(T)[3])
scheme <- as.character(commandArgs(T)[4])

#### 2. run meta_analysis_f.R ####
meta_analysis_f(
  sumstats = sumstats,
  out_dir = out_dir,
  outfile = outfile,
  scheme = scheme
)