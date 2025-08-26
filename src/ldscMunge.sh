#!/bin/bash
#SBATCH -t 00:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --partition=rome
#SBATCH --output=logs/ldscMunge_%A_out.txt
#SBATCH --error=logs/ldscMunge_%A_error.txt

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

echo "This script provides an example of running LDSC heritability and partitioned heritability analysis on UKB sumstats. Copy and modify as needed to suit your analysis"

helpString="
usage: sbatch ldscMunge.sh -f <gwas-sumstat-file> -e <pheno-name> [-o <output-directory>] [-m <munge-extra-options>] [-h <h2-extra-options>] [-l <ld-score-file>]\n
     <gwas-sumstat-file> is a file containing SNP-based GWAS summary statistics, with (at least) columns containing SNP, BETA/OR, N, A1/A2, and P \n
     <pheno-name> is the name of the phenotype for which GWAS results were produced. The results file will be saved with this same name.\n
     [output-directory] is an optional directory to save results to \n
     [munge-extra-options] optional additional parameters to submit to ldsc munging script, e.g. '--N-col=OBS_CT --a2 REF', to adjust how the sumstat file is read in. The text you enter will be pasted directly in the ldsc command. Default is to handle results file format produced by annotatePlinkResults.sh\n
     [h2-extra-options] optional additional parameters to submit to ldsc heritability script, e.g. '--pop-prev .01 --samp-prev 0.2', to adjust how the h2 calculation is run. The text you enter will be pasted directly in the ldsc command\n
     [1kG-population] is an optional 1000 Genomes population of UKB to use as a reference. Currently only EUR is supported; to use other populations you must provide your own file with ancestry-specific LD scores calculated\n
     [ld-score-file] is an optional file containing LD scores; if not specified, default is to use pre-calculated EUR scores from https://github.com/bulik/ldsc
"

if [[ $1 == "help" ]]; then
        echo -e $helpString
        exit 0
fi

while getopts f:o:e:m: option
do
 case "${option}"
 in
 f) SUMSTATFILE=$OPTARG;;
 o) OUTDIR=$OPTARG;;
 e) PHENO=$OPTARG;;
 m) MUNGEOPT=$OPTARG;;
 esac
done

if [[ ! -f "${SUMSTATFILE}" ]]; then
   echo "ERROR: Sumstat file ${SUMSTATFILE} does not exist"
   echo -e $helpString
   exit -1
fi

if [[ "${OUTDIR}" == "" ]]; then
   OUTDIR="${ANALYSIS_DIR}"
fi

if [[ "${MUNGEOPT}" == "" ]]; then
   MUNGEOPT="--N-col N --a2 A2 --snp rsID --ignore SNP"
fi

echo "Munging sumstats"

#Munge sumstats
${LDSC_PATH}/munge_sumstats.py \
  --sumstats ${SUMSTATFILE} \
  --merge-alleles ${LDSC_REF_PATH}/w_hm3.snplist \
  --out ${OUTDIR}/${PHENO}_munge ${MUNGEOPT}