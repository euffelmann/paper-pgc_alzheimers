import os

# Note that for 1kg data we used the output from the 1000g_eur and 1000g_all MAGMA implementation in FUMA
# In that case skip to run PoPS

########### run MAGMA ###########
# run MAGMA on reference population
# Define pop of interest
pop = 'pop_of_interest'

# define in/output directory
base = 'path/to/directory'

# path to the UKB reference files in plink .bim .bed .fam format
ref = f'dir_to_LD_reference_panels/ukb_{pop}'

if not os.path.exists(f'{base}/magma.genes.out'):
    cmd =  f'magma --bfile {ref} --pval {base}/sumstats.txt use=variant_id,p_value ncol=neff --gene-annot magma.genes.annot --out {base}/magma' 
    os.system(cmd)

#check if output exists, otherwise run
if not os.path.exists(f'{base}/magma.genes.out'):
    cmd =  f'magma --bfile {ref} --pval {base}/sumstats.txt use=variant_id,p_value ncol=neff --gene-annot magma.genes.annot --out {base}/magma' 
    os.system(cmd)



########### run MAGMA tissue type ###########
# run tissue type analysis
# define expression file, can be downloaded from https://zenodo.org/records/12635505
magma_exp = "path/to/gtex_v8_ts_avg_log2TPM.txt"

command = f'magma --gene-results {base}/magma.genes.raw --gene-covar {magma_exp} --model direction=greater condition-hide=Average --out {base}/magma_exp_gtex_v8_ts_avg_log2TPM'
os.system(command)


########### run PoPS ###########
# Be aware that this step uses considerable RAM 
# Features can be downloaded from https://zenodo.org/records/12635505
if not os.path.exists(f"{base}/PoPS_out_full.preds"):
        cmd = f"python path/to/PoPS/installation/pops.py \
        --gene_annot_path /path/to/Schipper_et_al/PoPS_features/PoPS_features_FUMA_compatible/gene_annots.txt \
        --feature_mat_prefix /path/to/Schipper_et_al/PoPS_features//PoPS_features_FUMA_compatible/pops_features_full_FUMA_compatible/features_munged/pops_features \
        --num_feature_chunks 116 \
        --magma_prefix {base}/magma \
        --control_features_path /path/to/Schipper_et_al/PoPS_features/pops_features_full_FUMA_compatible/control.features \
        --out_prefix {base}/PoPS_out_full"
