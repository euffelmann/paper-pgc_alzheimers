#!/bin/bash
#SBATCH -t 00:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=rome
#SBATCH --output=logs/ldscMunge_%A_out.txt
#SBATCH --error=logs/ldscMunge_%A_error.txt

# =============================================================
#  ldscMunge.sh
#  Munge GWAS summary statistics for use with LDSC.
#  Output: <PHENO>_munge.sumstats.gz, used as input to
#          ldscH2.sh and ldscRg.sh.
# =============================================================

# -------------------------------------------------------------
# CONFIG — update these paths for your environment
# -------------------------------------------------------------
LDSC_PATH=""          # path to ldsc directory, e.g. /path/to/ldsc
LDSC_REF_PATH=""      # path to LDSC reference files, e.g. /path/to/LDSC_reffiles
MODULE_YEAR="2022"
MODULE_PYTHON="Python/2.7.18-GCCcore-11.3.0-bare"
# -------------------------------------------------------------

if [[ -z "$LDSC_PATH" || -z "$LDSC_REF_PATH" ]]; then
  echo "ERROR: LDSC_PATH and LDSC_REF_PATH must be set in the CONFIG block."
  exit 1
fi

module load "$MODULE_YEAR"
module load "$MODULE_PYTHON"

ANALYSIS_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"

# -------------------------------------------------------------
# Usage
# -------------------------------------------------------------
usage() {
  echo "
Usage: sbatch ldscMunge.sh -f <sumstats-file> -e <pheno-name> [-o <output-dir>] [-m <munge-options>]

Required:
  -f  GWAS summary statistics file. Must contain columns for SNP,
      BETA/OR, N, A1/A2, and P.
  -e  Phenotype name. Output files will be prefixed with this name.

Optional:
  -o  Output directory (default: submission directory).
  -m  Additional options passed directly to munge_sumstats.py,
      e.g. '--N-col OBS_CT --a2 REF'.
      Default: '--N-col N --a2 A2 --snp rsID --ignore SNP'
"
}

if [[ "$1" == "help" ]]; then
  usage
  exit 0
fi

# -------------------------------------------------------------
# Parse arguments
# -------------------------------------------------------------
while getopts f:o:e:m: option; do
  case "${option}" in
    f) SUMSTATFILE=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    e) PHENO=$OPTARG ;;
    m) MUNGEOPT=$OPTARG ;;
  esac
done

if [[ ! -f "${SUMSTATFILE}" ]]; then
  echo "ERROR: Sumstats file '${SUMSTATFILE}' does not exist."
  usage
  exit 1
fi

OUTDIR="${OUTDIR:-$ANALYSIS_DIR}"
MUNGEOPT="${MUNGEOPT:---N-col N --a2 A2 --snp rsID --ignore SNP}"

# -------------------------------------------------------------
# Run munge_sumstats
# -------------------------------------------------------------
echo "Munging sumstats: ${SUMSTATFILE}"

${LDSC_PATH}/munge_sumstats.py \
  --sumstats      ${SUMSTATFILE} \
  --merge-alleles ${LDSC_REF_PATH}/w_hm3.snplist \
  --out           ${OUTDIR}/${PHENO}_munge \
  ${MUNGEOPT}

echo "Done. Output: ${OUTDIR}/${PHENO}_munge.sumstats.gz"
