#!/bin/bash
#SBATCH -t 00:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=rome
#SBATCH --output=logs/ldsc_h2_%A_out.txt
#SBATCH --error=logs/ldsc_h2_%A_error.txt

# =============================================================
#  ldscH2.sh
#  Estimate SNP heritability (h2) with LDSC.
#  Expects munged sumstats as input (output of ldscMunge.sh).
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
Usage: sbatch ldscH2.sh -f <munged-sumstats> -e <pheno-name> [-o <output-dir>] [-h <h2-options>] [-l <ld-score-dir>]

Required:
  -f  Munged GWAS sumstats file (<PHENO>_munge.sumstats.gz from ldscMunge.sh).
  -e  Phenotype name. Output files will be prefixed with this name.

Optional:
  -o  Output directory (default: submission directory).
  -h  Additional options passed directly to ldsc.py --h2,
      e.g. '--pop-prev 0.01 --samp-prev 0.2' for case-control traits.
  -l  Directory containing LD score files.
      Default: EUR LD scores from the LDSC reference directory.
"
}

if [[ "$1" == "help" ]]; then
  usage
  exit 0
fi

# -------------------------------------------------------------
# Parse arguments
# -------------------------------------------------------------
while getopts f:o:e:h:l: option; do
  case "${option}" in
    f) SUMSTATFILE=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    e) PHENO=$OPTARG ;;
    h) H2OPT=$OPTARG ;;
    l) LDSCOREFILE=$OPTARG ;;
  esac
done

if [[ ! -f "${SUMSTATFILE}" ]]; then
  echo "ERROR: Sumstats file '${SUMSTATFILE}' does not exist."
  usage
  exit 1
fi

OUTDIR="${OUTDIR:-$ANALYSIS_DIR}"
LDSCOREFILE="${LDSCOREFILE:-${LDSC_REF_PATH}/eur_w_ld_chr/}"

# -------------------------------------------------------------
# Run h2 estimation
# -------------------------------------------------------------
echo "Estimating h2 for: ${PHENO}"

${LDSC_PATH}/ldsc.py \
  --h2    ${SUMSTATFILE} \
  --ref-ld ${LDSCOREFILE} \
  --w-ld   ${LDSCOREFILE} \
  --out    ${OUTDIR}/${PHENO}_h2 \
  ${H2OPT}

echo "Done. Output: ${OUTDIR}/${PHENO}_h2.log"
