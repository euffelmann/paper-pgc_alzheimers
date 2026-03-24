#!/bin/bash
#SBATCH -t 00:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=rome
#SBATCH --output=logs/ldsc_rg_%A_out.txt
#SBATCH --error=logs/ldsc_rg_%A_error.txt

# =============================================================
#  ldscRg.sh
#  Estimate genetic correlations (rg) between a target trait
#  and one or more comparison traits using LDSC.
#  Expects munged sumstats as input (output of ldscMunge.sh).
# =============================================================

# -------------------------------------------------------------
# CONFIG — update these paths for your environment
# -------------------------------------------------------------
LDSC_PATH=""          # path to ldsc directory, e.g. /path/to/ldsc
LDSC_REF_PATH=""      # path to LDSC reference files, e.g. /path/to/LDSC_reffiles
MODULE_YEAR="2023"
MODULE_PYTHON="Python/2.7.18-GCCcore-12.3.0"
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
Usage: sbatch ldscRg.sh -f <munged-sumstats> -e <pheno-name> [-o <output-dir>] [-r <comparison-sumstats>] [-l <ld-score-dir>]

Required:
  -f  Munged GWAS sumstats for the target trait
      (<PHENO>_munge.sumstats.gz from ldscMunge.sh).
  -e  Phenotype name. Output files will be prefixed with this name.

Optional:
  -o  Output directory (default: submission directory).
  -r  Comma-separated list of additional munged sumstats to correlate
      against, e.g. '/path/to/PHENO2.sumstats.gz,/path/to/PHENO3.sumstats.gz'.
      These are appended directly to the --rg argument.
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
while getopts f:o:e:r:l: option; do
  case "${option}" in
    f) SUMSTATFILE=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    e) PHENO=$OPTARG ;;
    r) RGOPT=$OPTARG ;;
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
# Run rg estimation
# -------------------------------------------------------------
echo "Estimating rg for: ${PHENO}"

${LDSC_PATH}/ldsc.py \
  --rg     ${SUMSTATFILE},${RGOPT} \
  --ref-ld ${LDSCOREFILE} \
  --w-ld   ${LDSCOREFILE} \
  --out    ${OUTDIR}/${PHENO}_rg

# -------------------------------------------------------------
# Format output table
# -------------------------------------------------------------
echo "Formatting output table"

echo "rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se" \
  > ${OUTDIR}/${PHENO}_rg.table

grep "$RGOPT" ${OUTDIR}/${PHENO}_rg.log \
  | grep -v "Reading" \
  | grep -v "ERROR" \
  | awk '{print $3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' \
  | awk -F "/" '{print $NF}' \
  | sed '1d' \
  >> ${OUTDIR}/${PHENO}_rg.table

echo "Done. Output: ${OUTDIR}/${PHENO}_rg.log and ${OUTDIR}/${PHENO}_rg.table"
