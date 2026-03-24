src/ — Third-party executables and tools
=========================================

liftOver/
  Executable for lifting over genomic coordinates between reference
  builds (e.g. hg38 to hg19). Used via snp_modifyBuild_local.R to
  remap SNP positions in summary statistics that were originally
  reported on a different genome build.
  Downloaded from:
    https://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/liftOver  (macOS)
    https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver   (Linux)

plink1.9/
  PLINK v1.9 (stable beta 7.2, 11 Dec 2023). Used for genotype QC
  in 01_process_raw_genotype.R, including SNP and individual
  filtering, sex checks, Hardy-Weinberg filtering, relatedness
  estimation, and MDS-based population structure plots.
  Downloaded from: https://www.cog-genomics.org/plink/

plink2/
  PLINK v2. Used alongside PLINK 1.9 in 01_process_raw_genotype.R
  for steps that require PLINK 2 functionality, specifically
  filtering multi-allelic and invariant variants.
  Downloaded from: https://www.cog-genomics.org/plink/2.0/

METAL-2020-05-05/
  METAL (release 2020-05-05). Software for fixed-effects GWAS
  meta-analysis. Called by run_metal.R via meta_analysis_f.R to
  perform within-ancestry and cross-ancestry meta-analyses of
  case-control, proxy, and combined summary statistics.
  Downloaded from: https://github.com/statgen/METAL/releases

generic-metal/
  Older METAL release. Retained for reproducibility of earlier
  analyses. The current pipeline uses METAL-2020-05-05 instead.
  Downloaded from: https://csg.sph.umich.edu/abecasis/metal/download/
