#!/bin/bash
#SBATCH -t 00:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=rome
#SBATCH --output=logs/ldsc_rg_%A_out.txt
#SBATCH --error=logs/ldsc_rg_%A_error.txt

# source /home/ctgukbio/ctgukbio_prjs1606/programs/ukbioCTGscriptsVcurrentbeta/scriptSettings.sh

module load 2022
module load Python/2.7.18-GCCcore-11.3.0-bare

export HELPER_SCRIPTDIR="/projects/0/ctgukbio/programs/ukbioCTGscriptsVcurrentbeta/helperfiles"
export LDSC_PATH=/projects/0/ctgukbio/programs/ldsc
export LDSC_REF_PATH=/projects/0/ctgukbio/programs/LDSC_reffiles
ANALYSIS_DIR="$SLURM_SUBMIT_DIR"
if [[ "$ANALYSIS_DIR" == "" ]]; then
  ANALYSIS_DIR=`pwd`
fi

echo "This script provides an example of running LDSC genetic correlation analysis with UKB sumstats and existing sumstats for a set of physical/psychological traits with reliable large-scale GWAS results available. Copy and modify as needed to suit your analysis"

helpString="
usage: ldscRg.sh -f <gwas-sumstat-file> -e <pheno-name> [-o <output-directory>] [-r <rg-extra-options>] [-l <ld-score-file>]\n
<gwas-sumstat-file> is a file containing SNP-based GWAS summary statistics, AFTER munging via the previous LDSC step (i.e. the <PHENO>_munge.sumstats.gz file) \n
<pheno-name> is the name of the phenotype for which GWAS results were produced. The results file will be saved with this same name.\n
[output-directory] is an optional directory to save results to \n
[rg-extra-options] optional additional phenotypes to include in ldsc rg script, e.g. ',/dir/to/sumstats/PHENO.sumstats.gz'. Make sure to include the leading comma. The text you enter will be pasted directly in the ldsc command\n
[ld-score-file] is an optional file containing LD scores; if not specified, default is to use pre-calculated EUR scores from https://github.com/bulik/ldsc
"

if [[ $1 == "help" ]]; then
        echo -e $helpString
        exit 0
fi

while getopts f:o:p:e:r:l: option
do
 case "${option}"
 in
 f) SUMSTATFILE=$OPTARG;;
 o) OUTDIR=$OPTARG;;
 e) PHENO=$OPTARG;;
 r) RGOPT=$OPTARG;;
 l) LDSCOREFILE=$OPTARG;;
 esac
done

if [[ ! -f "${SUMSTATFILE}" ]]; then
   echo "ERROR: Sumstat file ${SUMSTATFILE} does not exist"
   echo -e $helpString
   exit -1
fi

if [[ "${LDSCOREFILE}" == "" ]]; then
   LDSCOREFILE="${LDSC_REF_PATH}/eur_w_ld_chr/"
fi

if [[ "${OUTDIR}" == "" ]]; then
   OUTDIR="${ANALYSIS_DIR}"
fi

echo "Beginning rg calculation"

${LDSC_PATH}/ldsc.py \
  --rg ${SUMSTATFILE},\
${RGOPT} \
  --ref-ld ${LDSCOREFILE} \
  --w-ld ${LDSCOREFILE} \
  --out ${OUTDIR}/${PHENO}_rg


echo "Formatting output table"
echo 'rg se z p h2_obs h2_obs_se h2_int h2_int_se gcov_int gcov_int_se' > ${OUTDIR}/${PHENO}_rg.table
grep $RGOPT ${OUTDIR}/${PHENO}_rg.log | grep -v Reading | grep -v ERROR | awk '{print $3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' | awk -F "/" '{print $NF}' | sed '1d' >>  ${OUTDIR}/${PHENO}_rg.table

echo "job finished"
