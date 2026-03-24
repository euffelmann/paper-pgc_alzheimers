#### 0. setup ####
library(here)
library(data.table)
library(dplyr)
library(openxlsx)
source(here("R/riskloc.R"))

#### 1. input ####
input <- as.character(commandArgs(T)[1])
output <- as.character(commandArgs(T)[2])
clump_kb <- as.integer(commandArgs(T)[3])
if (is.na(clump_kb)) {
  clump_kb <- 250
}

#### 2. define loci ####
sumstats <- fread(input)
loci <- riskloc(sumstats = sumstats, clump_kb = clump_kb)
fwrite(x = loci, file = output, quote = F, sep = "\t", na = NA)