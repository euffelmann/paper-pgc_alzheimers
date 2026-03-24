import os

annot ='/Path/to/FLAMES/Annotation_data/'
indexfile = '/path/to/infexfile/listing/finemapped/locus/locations' #(see FLAMES github)
pops = '/path/to/pops/file/PoPS_out_full.preds'
magma = '/path/to/magma/file/magma.genes.out'
magma_tissue = '/path/to/magma/file/magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out'

# This section is for running command line VEP and CADD (see FLAMES github)
# This can be ignored by removing the -t, -cf, -cv and -vc flags, but this might result in slow annotation if annotating large numbers of loci.
tabix='/path/to/tabix/executable/tabix'
cadd='/path/to/cadd/scores/CADD_scores_hg19'
vep='/path/to/VEP/installation/ensembl-vep/vep'
vep_cache= '/path/to/VEP/cache/.vep'

# run FLAMES annotate
cmd = f'python /path/to/FLAMES/installation/FLAMES.py annotate \
    -a "{annot} \
    -p {pops} \
    -m {magma} \
    -mt {magma_tissue} \
    -id {indexfile} \
    -pc PIP \
    -sc SNP \
    -t {tabix} \
    -cf {cadd} \
    -cv {vep} \
    -vc {vep_cache}'

os.system(cmd)

# run final FLAMES scoring
cmd = f'python /path/to/FLAMES/installation/FLAMES.py FLAMES \
    -id {indexfile} \
    -o /path/to/output/dir'

os.system(cmd)