#### 0. set-up ####
source("qc_plots.R")
library(data.table)
library(dplyr)
library(purrr)
library(here)
library(qqman)
library(stringr)

#### 1. input ####
sumstats_file <- as.character(commandArgs(T)[1])
outfile <- as.character(commandArgs(T)[2])
qc_dir <- as.character(commandArgs(T)[3])
ancestry <- as.character(commandArgs(T)[4])

#### 2. make qc plots ####
qc_plots(
  sumstats_file = sumstats_file,
  outfile = outfile,
  qc_dir = qc_dir,
  ancestry = ancestry
)