#!/bin/bash

#SBATCH --nodes=1
#SBATCH -t 00:30:00
#SBATCH --job-name=02_format_raw_sumstats
#SBATCH --output=logs/02_format_raw_sumstats_out_%A.txt
#SBATCH --error=logs/02_format_raw_sumstats_error_%A.txt
#SBATCH --partition=rome
#SBATCH --mem=50G

## load modules
module load 2022
module load R/4.2.1-foss-2022a

Rscript 02_format_raw_sumstats.R
