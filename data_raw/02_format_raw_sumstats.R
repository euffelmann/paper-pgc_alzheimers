# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                         #
#         script to format raw summary statistics         #
#         to a common GWAS catalog format                 #
#         (done manually)                                 #
#                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# - - - - - Table of Contents - - - - - - - - - - - - - - #
# - - - - - 00. set-up - - - - - - - - - - - - - - - - -  #
# - - - - - 01. Million Veteran Program (MVP) - - - - - - #
# - - - - - 02. Japanese Biobanks - - - - - - - - - - - - #
# - - - - - 03. Shigemizu et al. (2021) - - - - - - - - - #
# - - - - - 04. VUMC Amsterdam  - - - - - - - - - - - - - #
# - - - - - 05. 23andme - - - - - - - - - - - - - - - - - #
# - - - - - 06. IGAP  - - - - - - - - - - - - - - - - - - #
# - - - - - 07. Regeneron  - - - - - - - - - - - - - - -  #
# - - - - - 08. Genentech  - - - - - - - - - - - - - - -  #
# - - - - - 09. FinnGen  - - - - - - - - - - - - - - - -  #
# - - - - - 10. DeCODE - - - - - - - - - - - - - - - - -  #
# - - - - - 11. bioVU - - - - - - - - - - - - - - - - - - #
# - - - - - 12. Gothenburg - - - - - - - - - - - - - - -  #
# - - - - - 13. HUNT - - - - - - - - - - - - - - - - - -  #
# - - - - - 14. Estonian biobank - - - - - - - - - - - -  #
# - - - - - 15. GR@CE - - - - - - - - - - - - - - - - - - #
# - - - - - 16. STSA  - - - - - - - - - - - - - - - - - - #
# - - - - - 17. TwinGene  - - - - - - - - - - - - - - - - #
# - - - - - 18. UKB - - - - - - - - - - - - - - - - - - - #
# - - - - - 19. Lifelines - - - - - - - - - - - - - - - - #
# - - - - - 20. AllOfUs - - - - - - - - - - - - - - - - - #
# - - - - - 21. Eli Lilly, Johnson & Johnson, PROTECT - - #
# - - - - - 22. DemGene - - - - - - - - - - - - - - - - - #
# - - - - - 23. HUSK - - - - - - - - - - - - - - - - - -  #
# - - - - - 24. Copenhagen Hospital Bank (CHB) - - - - -  #
# - - - - - 25. EADB (Bellenguez 2021) - - - - - - - - -  #

#### 00. set-up ####

## load packages
library(data.table)
library(dplyr)
library(here)
library(purrr)
library(stringr)
library(bigsnpr)
source(here("R", "neff_f.R"))
source(here("R", "snp_modifyBuild_local.R"))

#### 01. Million Veteran Program (MVP) ####

## load summary statistics
for (sst in c("dadalz", "momalz", "60plus")) {
  
  ## set up variables
  if (sst == "dadalz") {
    outfile <- "mvplo_proxy_father_noag_afr"
    n_case <- 2256
    n_control <- 45970
    traitDescription <- "proxy-AD in fathers"
    proxy_casecontrol <- "proxy"
  } else if (sst == "momalz") {
    outfile <- "mvplo_proxy_mother_noag_afr"
    n_case <- 4385
    n_control <- 45970
    traitDescription <- "proxy-AD in mothers"
    proxy_casecontrol <- "proxy"
  } else if (sst == "60plus") {
    outfile <- "mvplo_case_control_noag_afr"
    n_case <- 4012
    n_control <- 18435
    traitDescription <-
      "AD case-control in individuals with age >= 60"
    proxy_casecontrol <- "case_control"
  }
  
  ## only run if sumstats don't exist
  if (!file.exists(here(
    "data_raw",
    "sumstats",
    proxy_casecontrol,
    paste0(outfile, ".txt.gz")
  ))
  ) {
    
    ## format sumstats
    print(paste0("load & format sumstats: ", outfile))
    # some of the columns have missing data, which leads to issues when 
    # loading the data. I only load lines that are complete.
    fread(
      cmd =
        paste0(
          "gunzip -c ",
          here("data_raw", "cohorts", "mvp_logue/"),
          "afr_",
          sst,
          "_2021_pgc.txt.gz | awk '{if(NF==18)print}'"
        ),
      header = T
    ) %>% 
      mutate(
        neff = neff_f(
          n_case = N_case,
          n_control = N_ctrl,
          proxy_casecontrol = proxy_casecontrol
        ),
        beta = log(OR)
      ) %>% 
      rename(
        chromosome = "#CHROM",
        base_pair_location = POS,
        effect_allele = A1,
        other_allele = AX,
        odds_ratio = OR,
        standard_error = "LOG(OR)_SE",
        effect_allele_frequency = A1_FREQ,
        p_value = P,
        info = MACH_R2,
        n_case = N_case,
        n_control = N_ctrl
      ) %>%
      mutate(
        variant_id = paste(
          chromosome, 
          base_pair_location, 
          other_allele, 
          effect_allele, 
          sep = ":")
      ) %>%
      select(
        variant_id,
        chromosome,
        base_pair_location,
        other_allele,
        effect_allele,
        odds_ratio,
        beta,
        standard_error,
        effect_allele_frequency,
        p_value,
        info,
        neff,
        n_case,
        n_control
      ) %>%
      arrange(chromosome, base_pair_location) %>%
      {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
      fwrite(
        file = here(
          "data_raw",
          "sumstats",
          proxy_casecontrol,
          paste0(outfile, ".txt.gz")
        ),
        quote = F,
        sep = "\t",
        na = "NA",
        row.names = F
      )
    
    ## save readme file
    readme <-
      file(here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, "_readme.txt")
      ))
    writeLines(c(
      paste0("Date: ", Sys.Date()),
      "- genomeAssembly: GRCh37/hg19",
      paste0("- traitDescription: ", traitDescription),
      paste0("- sampleSize: ", n_case + n_control),
      paste0("- caseCount: ", n_case),
      paste0("- controlCount: ", n_control),
      "- caseControlStudy: Yes",
      "- sampleAncestry: African",
      "- genotypingTechnology: MVP 1.0 custom Axiom array",
      "- analysisSoftware: logistic regression models implemented in PLINK 2.0",
      "- imputationPanel: agr",
      "- imputationSoftware: African Genome Rsources imputation panel",
      paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
      "- ancestryMethod: Subjects for this study included AFR MVP participants ",
      "     as determined by the genetically informed harmonized ancestry and ",
      "     race/ethnicity (HARE) method. This method combines self-report and ",
      "     genetically determined principal components in a machine learning model ",
      "     to assign participants to one of four groups: AFR, EUR, ASN (East Asian), ",
      "     and HIS (Hispanic). Individuals with indeterminate ancestry and ",
      "     individuals whose self-reported Race did not match genetically informed ",
      "     ancestry were excluded from classification",
      "- sortedByGenomicLocation: Yes",
      "- dosageOrHardcalls: Dosage"
    ), 
    readme)
    close(readme)
    
  } else {
    
    ## print if code is skipped because file already exists
    print(paste0(outfile, ": already exists")) 
    
  }
  
  ## clear workspace
  rm(list=setdiff(ls(), c("neff_f", "snp_modifyBuild_local")))
  
}


#### 02. Japanese Biobanks #####

## function to load, format, process, save, and plot jpb-sumstats
jpb_sumstats <- function(age_cutoff, covar, parent, proxy_casecontrol) {
  
  ## set file names
  infile <-
    here(
      "data_raw",
      "cohorts",
      "japanese_biobanks",
      paste0(proxy_casecontrol, age_cutoff),
      paste0("regenie_Age", age_cutoff, "_", parent, "_", covar),
      paste0("out_step2_chrALL_", parent, ".regenie")
    )
  
  if (covar == "withoutAge") {covar <- "noag"}
  if (covar == "AGE") {covar <- "main"}
  outfile <-
    paste0(
      "jpbyu_",
      proxy_casecontrol,
      "_",
      parent,
      "_",
      tolower(covar),
      "_eas"
    )
  
  ## only run if outfile does not exist
  if (!file.exists(here(
    "data_raw",
    "sumstats",
    proxy_casecontrol,
    paste0(outfile, ".txt.gz")
  ))) {
    
    ## load & format sumstats
    print(paste0("load & format sumstats: ", infile))
    fread(paste0(infile)) %>%
      {max(.$N_CASES, na.rm = T) ->> n_case;.} %>%
      {max(.$N_CONTROLS, na.rm = T) ->> n_control;.} %>%
      mutate(
        neff = neff_f(
          n_case = N_CASES,
          n_control = N_CONTROLS,
          proxy_casecontrol = proxy_casecontrol
        )
      ) %>% 
      rename(
        chromosome = CHROM,
        base_pair_location = GENPOS,
        other_allele = ALLELE0,
        effect_allele = ALLELE1,
        beta = BETA,
        standard_error = SE,
        effect_allele_frequency = A1FREQ_CONTROLS,
        p_value = P,
        info = INFO,
        neff = neff,
        n_case = N_CASES,
        n_control = N_CONTROLS
      )  %>%
      mutate(
        variant_id = paste(
          chromosome, 
          base_pair_location, 
          other_allele, 
          effect_allele, 
          sep = ":")
      ) %>%
      select(
        variant_id,
        chromosome,
        base_pair_location,
        other_allele,
        effect_allele,
        beta,
        standard_error,
        effect_allele_frequency,
        p_value,
        info,
        neff,
        n_case,
        n_control
      ) %>%
      arrange(chromosome, base_pair_location) %>%
      {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
      fwrite(
        file = here(
          "data_raw",
          "sumstats",
          proxy_casecontrol,
          paste0(outfile, ".txt.gz")
        ),
        quote = F,
        sep = "\t",
        na = "NA",
        row.names = F
      )
    
    ## save readme file
    readme <-
      file(here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, "_readme.txt")
      ))
    writeLines(c(
      paste0("Date: ", Sys.Date()),
      "- genomeAssembly: GRCh37/hg19",
      "- traitDescription: 
          Individuals who reported their parents to have 
          Alzheimer’s Disease / Dementia",
      paste0("- sampleSize: ", n_case + n_control),
      paste0("- caseCount: ", n_case),
      paste0("- controlCount: ", n_control),
      "- caseControlStudy: Yes",
      "- sampleAncestry: Japanese",
      "- genotypingTechnology: Japonica Array v.2 ",
      "- analysisSoftware: regenie",
      "- imputationPanel: 
            a combined reference panel of 1000 Genomes Project Phase 3 (n = 2,504) 
            and population-specific WGS data (that is, 3.5KJPNv2; n = 3,552)",
      "- imputationSoftware: 
            genotype data were prephased by using SHAPEIT2 software (r837), 
            and imputed by using IMPUTE4 software (r300.3)",
      paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
      "- ancestryMethod: removed non-East Asian outlier samples by genomic PCA",
      "- sortedByGenomicLocation: Yes",
      "- dosageOrHardcalls: Dosage"
    ), 
    readme)
    close(readme)
  } else {
    
    ## print if code is skipped because file already exists
    print(paste0(outfile, ": already exists")) 
    
  }
  
  ## clear workspace
  rm(list=setdiff(ls(), c("neff_f", "snp_modifyBuild_local")))
}

## parameters defining different sumstats
params <- expand.grid(
  age_cutoff = c("40"),
  covar = c("AGE", "APOE", "withoutAge"),
  parent = c("father", "mother"),
  proxy_casecontrol = "proxy",
  stringsAsFactors = F
)

params %>%
  pmap(~ jpb_sumstats(
    age_cutoff = ..1,
    covar = ..2,
    parent = ..3,
    proxy_casecontrol = ..4
  ))


#### 03. Shigemizu et al. (2021) ####

## only run if sumstats don't exist

## set up variables
outfile <- "indsh_case_control_main_eas" # ind for independent (i.e. no Biobank)
proxy_casecontrol <- "case_control"
n_case <- 3962
n_control <- 4074

if (!file.exists(here(
  "data_raw",
  "sumstats",
  "case_control",
  "indsh_case_control_main_eas.txt.gz"
))
) {
  
  ## format sumstats
  print(paste0("load & format sumstats: ", outfile))
  fread(here("data_raw", "cohorts", "shigemizu_2021", "NCGG_AD_GWAS2.txt"),
        header = T) %>%
    mutate(neff = neff_f(
             n_case = NMISS_A,
             n_control = NMISS_U,
             proxy_casecontrol = proxy_casecontrol),
           beta = log(OR)) %>%
    rowwise() %>%
    mutate(info = min(Info_NCGG, Info_Niigata, na.rm = TRUE)) %>%
    ungroup() %>%
    rename(
      rsid = SNP,
      chromosome = CHR,
      base_pair_location = BP,
      other_allele = A2,
      effect_allele = A1,
      odds_ratio = OR,
      standard_error = SE,
      effect_allele_frequency = MAF_U,
      p_value = P,
      n_case = NMISS_A,
      n_control = NMISS_U
    ) %>%
    mutate(
      variant_id = paste(
        chromosome, 
        base_pair_location, 
        other_allele, 
        effect_allele, 
        sep = ":")
    ) %>%
    select(
      variant_id,
      chromosome,
      base_pair_location,
      effect_allele,
      other_allele,
      odds_ratio,
      beta,
      standard_error,
      effect_allele_frequency,
      p_value,
      info,
      neff,
      n_case,
      n_control,
      rsid
    ) %>%
    arrange(chromosome, base_pair_location) %>%
    {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
    fwrite(
      file = here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, ".txt.gz")
      ),
      quote = F,
      sep = "\t",
      na = "NA",
      row.names = F
    )
  
  ## save readme file
  readme <-
    file(here(
      "data_raw",
      "sumstats",
      proxy_casecontrol,
      paste0(outfile, "_readme.txt")
    ))
  writeLines(c(
    paste0("Date: ", Sys.Date()),
    "- genomeAssembly: GRCh37/hg19",
    "- traitDescription: AD in East Asian (Japanese) individuals",
    paste0("- sampleSize: ", n_case + n_control),
    paste0("- caseCount: ", n_case),
    paste0("- controlCount: ", n_control),
    "- caseControlStudy: Yes",
    "- sampleAncestry: East Asians (Japanese)",
    "- genotypingTechnology: 
          Affymetrix Japonica Array22 for the NCGG subjects and 
          Affymetix GeneChip 6.0 microarrays for Niigata subjects",
    "- analysisSoftware: PLINK",
    "- imputationPanel: 
            3.5 K Japanese reference panel developed by 
            Tohoku Medical Megabank Organization",
    "- imputationSoftware: IMPUTE2",
    paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
    "- ancestryMethod: ?",
    "- sortedByGenomicLocation: Yes",
    "- dosageOrHardcalls: ?"
  ), 
  readme)
  close(readme)
} else {
  
  ## print if code is skipped because file already exists
  print(paste0(outfile, ": already exists")) 
  
}

## clear workspace
rm(list=setdiff(ls(), c("neff_f", "snp_modifyBuild_local")))


#### 04. VUMC Amsterdam ####

## function to load, format, process, save, and plot jpb-sumstats
vumc_sumstats <- function(analysis, proxy_casecontrol) {
  
  ## set file names
  infile <-
    Sys.glob(here(
      "data_raw",
      "cohorts",
      "vumc",
      paste0("*_", analysis, ".glm.logistic.hybrid.gz")
    ))
  if (analysis == "without_age") {analysis <- "noag"}
  if (analysis == "female") {analysis <- "fema"}
  if (analysis == "add_APOE") {analysis <- "apoe"}
  outfile <-
    paste0("vumcl_", proxy_casecontrol, "_", analysis, "_eur")
  
  ## only run if outfile does not exist
  if (!file.exists(here(
    "data_raw",
    "sumstats",
    proxy_casecontrol,
    paste0(outfile, ".txt.gz")
  ))) {
    
    ## load & format sumstats
    print(paste0("load, liftover (hg38 to GRCh37/hg19) & format sumstats: ", infile))
    fread(infile) %>% 
      mutate(
        variant_id_hg38 = paste(
          chromosome, 
          base_pair_location, 
          other_allele, 
          effect_allele, 
          sep = ":")
      ) %>%
      ## change column names for liftover
      rename(chr = chromosome,
             pos = base_pair_location) %>%
      snp_modifyBuild_local(liftOver = here("src", "liftOver"),
                      from = "hg38",
                      to = "hg19",
                      check_reverse = T) %>%
      rename(chromosome = chr,
             base_pair_location = pos) %>%
      {max(.$n_case) ->> n_case;.} %>%
      {max(.$n_control) ->> n_control;.} %>%
      mutate(
        neff = neff_f(n_case = n_case,
                      n_control = n_control,
                      proxy_casecontrol = proxy_casecontrol), 
        info = str_extract(info, "R2=([\\d\\.]+)") %>%
          str_replace_all("R2=", "") %>%
          as.numeric()
      ) %>% 
      mutate(
        variant_id = paste(
          chromosome, 
          base_pair_location, 
          other_allele, 
          effect_allele, 
          sep = ":")
      ) %>%
      select(
        variant_id,
        chromosome,
        base_pair_location,
        other_allele,
        effect_allele,
        beta,
        standard_error,
        effect_allele_frequency,
        p_value,
        info,
        neff,
        n_case,
        n_control,
        variant_id_hg38
      ) %>%
      arrange(chromosome, base_pair_location) %>%
      {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
      fwrite(
        file = here(
          "data_raw",
          "sumstats",
          proxy_casecontrol,
          paste0(outfile, ".txt.gz")
        ),
        quote = F,
        sep = "\t",
        na = "NA",
        row.names = F
      )
    
    ## save readme file
    readme <-
      file(here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, "_readme.txt")
      ))
    writeLines(c(
      paste0("Date: ", Sys.Date()),
      "- genomeAssembly: GRCh37/hg19",
      "- traitDescription: AD in European (Dutch) individuals",
      paste0("- sampleSize: ", n_case + n_control),
      paste0("- caseCount: ", n_case),
      paste0("- controlCount: ", n_control),
      "- caseControlStudy: Yes",
      "- sampleAncestry: European (Dutch)",
      "- genotypingTechnology: 
            Illumina Infinium Global Screening Array (GSA, GSAsharedCUSTOM_24+v1.0)",
      "- analysisSoftware: Plink2",
      "- imputationPanel: TOPMED",
      "- imputationSoftware: 
            Michigan Imputation server (https://imputation.biodatacatalyst.nhlbi.nih.gov/). 
            The server uses EAGLE (v2.4) to phase data and Minimac4 to perform 
            genotype imputation to the reference panel (version r2).",
      paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
      "- ancestryMethod: ?",
      "- sortedByGenomicLocation: Yes",
      "- dosageOrHardcalls: Dosage",
      "- notes: ~29,000 SNPs cannot be mapped with liftOver"
    ), 
    readme)
    close(readme)
  } else {
    
    ## print if code is skipped because file already exists
    print(paste0(outfile, ": already exists")) 
    
  }
}

## parameters defining different sumstats
params <- expand.grid(
  analysis = c("main", "male", "female", "without_age", "add_APOE"),
  proxy_casecontrol = "case_control",
  stringsAsFactors = F
)

params %>%
  pmap(~ vumc_sumstats(
    analysis = ..1,
    proxy_casecontrol = ..2
  ))

## clear workspace
rm(list=setdiff(ls(), c("neff_f", "snp_modifyBuild_local")))


#### 05. 23andme ####

## set file names
proxy_casecontrol <- "case_control"

infile <-
  Sys.glob(here(
    "data_raw",
    "cohorts",
    "23andme",
    "23andMe.txt.gz"
  ))

outfile <-
  paste0("twent_", proxy_casecontrol, "_main_eur")

## only run if outfile does not exist
if (!file.exists(here(
  "data_raw",
  "sumstats",
  proxy_casecontrol,
  paste0(outfile, ".txt.gz")
))) {
  
  ## load file with info scores
  print("load 23andme info file ...")
  info_imp <-
    fread(here("data_raw", "cohorts", "23andme", "im_snp_stat.txt.gz"))
  
  ## load 1000G reference data to add rsIDs
  if (!exists("g1000_eur")) {
    print("Load 1000G eur reference file ...")
    g1000_eur <- fread(here("data", "reference_data", "g1000_eur", "g1000_eur.bim"))
  }
  
  ## load & format sumstats
  print(paste0("load & format sumstats: ", infile))
  fread(infile) %>% 
    {max(.$N_cases) ->> n_case;.} %>%
    {max(.$N_controls) ->> n_control;.} %>%
    mutate(
      neff = neff_f(n_case = N_cases,
                    n_control = N_controls,
                    proxy_casecontrol = proxy_casecontrol)
    ) %>% 
    rename(other_allele = Non_Effect_allele,
           effect_allele = Effect_allele,
           chromosome = CHR,
           base_pair_location = POS_GRCh37,
           beta = BETA,
           standard_error = SE,
           effect_allele_frequency = Effect_allele_freq,
           p_value = P,
           n_case = N_cases,
           n_control = N_controls) %>%
    mutate(
      variant_id = paste(
        chromosome, 
        base_pair_location, 
        other_allele, 
        effect_allele, 
        sep = ":")
    ) %>%
    ## to get rsids to be able to join with info_imp
    left_join(g1000_eur, c("chromosome" = "V1", "base_pair_location" = "V4")) %>%
    left_join(info_imp, c("V2" = "assay.name")) %>%
    rename(rsid = V2,
           info = min.rsqr) %>%
    select(
      variant_id,
      chromosome,
      base_pair_location,
      other_allele,
      effect_allele,
      beta,
      standard_error,
      effect_allele_frequency,
      p_value,
      info,
      neff,
      n_case,
      n_control,
      rsid
    ) %>%
    arrange(chromosome, base_pair_location) %>%
    {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
    fwrite(
      file = here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, ".txt.gz")
      ),
      quote = F,
      sep = "\t",
      na = "NA",
      row.names = F
    )
  
  ## save readme file
  readme <-
    file(here(
      "data_raw",
      "sumstats",
      proxy_casecontrol,
      paste0(outfile, "_readme.txt")
    ))
  writeLines(c(
    paste0("Date: ", Sys.Date()),
    "- genomeAssembly: GRCh37/hg19",
    "- traitDescription: Proxy AD in individuals of European ancestry",
    paste0("- sampleSize: ", n_case + n_control),
    paste0("- caseCount: ", n_case),
    paste0("- controlCount: ", n_control),
    "- caseControlStudy: Yes",
    "- sampleAncestry: European",
    "- genotypingTechnology: 
          DNA extraction and genotyping were performed on saliva samples by 
          National Genetics Institute (NGI), a CLIA licensed clinical laboratory 
          and a subsidiary of Laboratory Corporation of America. Samples were 
          genotyped on one of five genotyping platforms. The v1 and v2 platforms 
          were variants of the Illumina HumanHap550+ BeadChip, including about 
          25,000 custom SNPs selected by 23andMe, with a total of about 
          560,000 SNPs. The v3 platform was based on the Illumina OmniExpress+ BeadChip, 
          with custom content to improve the overlap with our v2 array, with a 
          total of about 950,000 SNPs. The v4 platform was a fully customized array, 
          including a lower redundancy subset of v2 and v3 SNPs with additional 
          coverage of lower-frequency coding variation, and about 570,000 SNPs. 
          The v5 platform, in current use, is an Illumina Infinium Global Screening 
          Array (~640,000 SNPs) supplemented with ~50,000 SNPs of custom content. 
          This array was specifically designed to better capture global genetic 
          diversity and to help standardize the platform for genetic research. ",
    "- analysisSoftware: ?",
    "- imputationPanel: 
          Variants were imputed in two separated imputation reference panels.  
          For the first one, we combined the May 2015 release of the 1000 Genomes 
          Phase 3 haplotypes3 with the UK10K imputation reference panel4 to create 
          a single unified panel. To do this, multiallelic sites with N alternate 
          alleles were split into N separate biallelic sites. We then removed any 
          site whose minor allele appeared in only one sample. For each chromosome, 
          we used Minimac35 to impute the reference panels against each other, 
          reporting the best-guess genotype at each site. This gave us calls for 
          all samples over a single unified set of variants. We then joined these 
          together to get, for each chromosome, a single file with phased calls at 
          every site for 6,285 samples. Throughout, we treated structural variants 
          and small indels in the same way as SNPs. We used the Human Reference 
          Consortium (HRC) as the second imputation reference panel6. It consists 
          of 32,488 samples for 39,235,157 SNPs.",
    "- imputationSoftware: Minimac3",
    paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
    "- ancestryMethod: 
          Individuals were assigned ancestry by first partitioning the phased 
          genomic data into short windows of about 300 SNPs. Within each window,
          a support vector machine (SVM) classified individual haplotypes into 
          one of 31 reference populations (https://www.23andme.com/ancestrycomposition-guide/). 
          The SVM classifications are then fed into a hidden Markov model (HMM) 
          that accounts for switch errors and incorrect assignments, and gives 
          probabilities for each reference population in each window. Finally, 
          we used simulated admixed individuals to recalibrate the HMM 
          probabilities so that the reported assignments are consistent with the 
          simulated admixture proportions. Europeans were defined as those with 
          ancestry probabilities of European + Middle Eastern > 0.97 and 
          European > 0.90.",
    "- sortedByGenomicLocation: Yes",
    "- dosageOrHardcalls: Dosage"
  ), 
  readme)
  close(readme)
} else {
  
  ## print if code is skipped because file already exists
  print(paste0(outfile, ": already exists")) 
  
}

## clean up workspace
rm(list=setdiff(ls(), c("neff_f", "snp_modifyBuild_local")))


#### 06. IGAP ####

## set file names
proxy_casecontrol <- "case_control"

infile <-
  Sys.glob(here(
    "data_raw",
    "cohorts",
    "igap",
    "IGAP.txt.gz"
  ))

outfile <-
  paste0("xigap_", proxy_casecontrol, "_main_eur")

# see Kunkle_etal_2019_IGAP_summary_statistics_README.docx
n_case <- 21982 
n_control <- 41944 

## only run if outfile does not exist
if (!file.exists(here(
  "data_raw",
  "sumstats",
  proxy_casecontrol,
  paste0(outfile, ".txt.gz")
))) {
  
  ## load & format sumstats
  print(paste0("load & format sumstats: ", infile))
  fread(infile) %>% 
    mutate(
      n_case = n_case, 
      n_control = n_control) %>%
    mutate(
      neff = neff_f(n_case = n_case,
                    n_control = n_control,
                    proxy_casecontrol = proxy_casecontrol)
    ) %>% 
    rename(other_allele = Non_Effect_allele,
           effect_allele = Effect_allele,
           chromosome = CHR,
           base_pair_location = POS_GRCh37,
           beta = BETA,
           standard_error = SE,
           effect_allele_frequency = Effect_allele_freq,
           p_value = P,
           rsid = SNP) %>%
    mutate(
      variant_id = paste(
        chromosome, 
        base_pair_location, 
        other_allele, 
        effect_allele, 
        sep = ":"),
      rsid = ifelse(str_detect(rsid, pattern = "rs"), rsid, NA)
    ) %>%
    select(
      variant_id,
      chromosome,
      base_pair_location,
      other_allele,
      effect_allele,
      beta,
      standard_error,
      effect_allele_frequency,
      p_value,
      neff,
      n_case,
      n_control,
      rsid
    ) %>%
    arrange(chromosome, base_pair_location) %>%
    {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
    fwrite(
      file = here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, ".txt.gz")
      ),
      quote = F,
      sep = "\t",
      na = "NA",
      row.names = F
    )
  
  ## save readme file
  readme <-
    file(here(
      "data_raw",
      "sumstats",
      proxy_casecontrol,
      paste0(outfile, "_readme.txt")
    ))
  writeLines(c(
    paste0("Date: ", Sys.Date()),
    "- genomeAssembly: GRCh37/hg19",
    "- traitDescription: AD in individuals of European ancestry (stage 1)",
    paste0("- sampleSize: ", n_case + n_control),
    paste0("- caseCount: ", n_case),
    paste0("- controlCount: ", n_control),
    "- caseControlStudy: Yes",
    "- sampleAncestry: European",
    "- genotypingTechnology: Many, because meta-analysis",
    "- analysisSoftware: Mostly SNPTest",
    "- imputationPanel: 1000 Genomes",
    "- imputationSoftware: Mostly IMPUTE2, some cohorts Minimac3",
    paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
    "- ancestryMethod: 
          Individuals with non-European ancestry according to principal 
          components analysis of ancestry-informative markers were excluded 
          from the further analysis.",
    "- sortedByGenomicLocation: Yes",
    "- dosageOrHardcalls: Dosage"
  ), 
  readme)
  close(readme)
} else {
  
  ## print if code is skipped because file already exists
  print(paste0(outfile, ": already exists")) 

}


#### 07. Regeneron ####

cohorts <- c("COLORADO",
             "GHS",
             "INDIANA-CHALASANI",
             "MAYO-CLINIC",
             "SINAI",
             "UCLA",
             "UPENN-PMBB")

proxy_casecontrol <- "case_control"

for (cohort in cohorts) {
  
  if (cohort %in% c("UPENN-PMBB", "SINAI", "INDIANA-CHALASANI")) {
    geno_array <- "Illumina Global Screening Array genotyping chip"
  } else if (cohort == "GHS") {
    geno_array <- "GHS participants genotyped using either the Illumina Infinium 
                   OmniExpressExome or the Global Screening Array"
  } else {
    geno_array <- "targeted genomic sequencing using the twist diversity SNP 
                   panel, followed by the multipoint refinement using GLIMPSE"
  }
  
  if (cohort == "GHS") {name <- "rege1"}
  if (cohort == "INDIANA-CHALASANI") {name <- "rege2"}
  if (cohort == "MAYO-CLINIC") {name <- "rege3"}
  if (cohort == "SINAI") {name <- "rege4"}
  if (cohort == "UCLA") {name <- "rege5"}
  if (cohort == "UPENN-PMBB") {name <- "rege6"}
  if (cohort == "COLORADO") {name <- "rege7"}
  
  infile <-
    Sys.glob(here(
      "data_raw",
      "cohorts",
      "regeneron",
      paste0("*.", cohort, "*.gz")
    ))
  
  outfile <-
    paste0(name, "_", proxy_casecontrol, "_main_eur")
    
  ## only run if outfile does not exist
  if (!file.exists(here(
    "data_raw",
    "sumstats",
    proxy_casecontrol,
    paste0(outfile, ".txt.gz")
  ))) {
    
    ## load & format sumstats
    print(paste0("load & format sumstats: ", infile))
    fread(infile) %>%
      # liftOver requires chr23 to be coded as X
      mutate(Chr = ifelse(Chr == 23, "X", Chr),
             variant_id_hg38 = paste(
               Chr, 
               Pos, 
               Ref, 
               Alt, 
               sep = ":")) %>%
      rename(chr = Chr,
             pos = Pos) %>%
      snp_modifyBuild_local(liftOver = here("src", "liftOver"),
                      from = "hg38",
                      to = "hg19",
                      check_reverse = T) %>%
      mutate(
        beta = str_extract(Info, "REGENIE_BETA=(-?[\\d\\.]+)") %>%
          str_replace_all("REGENIE_BETA=", "") %>%
          as.numeric(),
        standard_error = str_extract(Info, "REGENIE_SE=([\\d\\.]+)") %>%
          str_replace_all("REGENIE_SE=", "") %>%
          as.numeric(),
        info = str_extract(Info, "INFO=([\\d\\.]+)") %>%
          str_replace_all("INFO=", "") %>%
          as.numeric(),
        neff = neff_f(n_case = Num_Cases, 
                      n_control = Num_Controls, 
                      proxy_casecontrol = proxy_casecontrol)
      ) %>% 
      rename(chromosome = chr,
             base_pair_location = pos,
             effect_allele = Alt,
             other_allele = Ref,
             odds_ratio = Effect,
             effect_allele_frequency = AAF,
             p_value = Pval,
             n_case = Num_Cases,
             n_control = Num_Controls) %>%
      mutate(variant_id = paste(
        chromosome,
        base_pair_location,
        other_allele,
        effect_allele,
        sep = ":"
      )) %>% 
      {max(.$n_case) ->> n_case_glob;.} %>%
      {max(.$n_control) ->> n_control_glob;.} %>%
      select(
        variant_id,
        chromosome,
        base_pair_location,
        other_allele,
        effect_allele,
        beta,
        odds_ratio,
        standard_error,
        effect_allele_frequency,
        p_value,
        info,
        neff,
        n_case,
        n_control,
        variant_id_hg38
      ) %>%
      arrange(chromosome, base_pair_location) %>%
      {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
      fwrite(
        file = here(
          "data_raw",
          "sumstats",
          proxy_casecontrol,
          paste0(outfile, ".txt.gz")
        ),
        quote = F,
        sep = "\t",
        na = "NA",
        row.names = F
      )
    
    ## save readme file
    readme <-
      file(here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, "_readme.txt")
      ))
    writeLines(c(
      paste0("Date: ", Sys.Date()),
      "- genomeAssembly: GRCh37/hg19",
      "- traitDescription: alzheimers_disease_broad in europeans",
      paste0("- sampleSize: ", n_case_glob + n_control_glob),
      paste0("- caseCount: ", n_case_glob),
      paste0("- controlCount: ", n_control_glob),
      "- caseControlStudy: Yes",
      "- sampleAncestry: European",
      paste0("- genotypingTechnology: ", geno_array),
      "- analysisSoftware: REGENIE",
      "- imputationPanel: TOPMed",
      "- imputationSoftware: MaCH",
      paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
      "- ancestryMethod: 
            We calculated PCs for the HapMap3 samples on which each of our samples 
            was projected. We trained a kernel density estimator (KDE) using the 
            HapMap3 PCs and used the KDEs to calculate the likelihood of a given 
            sample belonging to each of the five continental ancestry groups. 
            Individuals were assigned to European group if their likelihood of 
            belonging to European ancestry was > 0.3.",
      "- sortedByGenomicLocation: Yes",
      "- dosageOrHardcalls: Dosage",
      "- notes:
            ~38,000 variants could not be properly mapped using liftover."
    ), 
    readme)
    close(readme)
  } else {
    
    ## print if code is skipped because file already exists
    print(paste0(outfile, ": already exists")) 
    
  }
  
}

## clear workspace
rm(list=setdiff(ls(), c("neff_f", "snp_modifyBuild_local")))


#### 08. Genentech ####

analyses <- c("main", "male", "fema", "apoe", "without_age")

proxy_casecontrol <- "case_control"

for (analysis in analyses) {
  
  infile <-
    Sys.glob(here(
      "data_raw",
      "cohorts",
      "genentech",
      paste0(analysis, "*.gz")
    ))
  
  if (analysis == "without_age") {analysis <- "noag"}
  
  outfile <-
    paste0("genen_", proxy_casecontrol, "_", analysis, "_eur")
  
  ## only run if outfile does not exist
  if (!file.exists(here(
    "data_raw",
    "sumstats",
    proxy_casecontrol,
    paste0(outfile, ".txt.gz")
  ))) {
    
    ## load & format sumstats
    print(paste0("load & format sumstats: ", infile))
    fread(infile) %>%
      mutate(
        variant_id_hg38 = paste(
          chromosome, 
          base_pair_location, 
          other_allele, 
          effect_allele, 
          sep = ":")
      ) %>%
      rename(chr = chromosome,
             pos= base_pair_location) %>%
      snp_modifyBuild_local(liftOver = here("src", "liftOver"),
                      from = "hg38",
                      to = "hg19",
                      check_reverse = T) %>%
      rename(chromosome = chr,
             base_pair_location = pos) %>%
      {max(.$n_case) ->> n_case;.} %>%
      {max(.$n_control) ->> n_control;.} %>%
      mutate(
        neff = neff_f(n_case = n_case,
                      n_control = n_control,
                      proxy_casecontrol = proxy_casecontrol),
        variant_id = paste(
          chromosome, 
          base_pair_location, 
          other_allele, 
          effect_allele, 
          sep = ":")) %>%
      select(
        variant_id,
        chromosome,
        base_pair_location,
        other_allele,
        effect_allele,
        beta,
        standard_error,
        effect_allele_frequency,
        p_value,
        info,
        neff,
        n_case,
        n_control,
        variant_id_hg38
      ) %>%
      arrange(chromosome, base_pair_location) %>%
      {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
      fwrite(
        file = here(
          "data_raw",
          "sumstats",
          proxy_casecontrol,
          paste0(outfile, ".txt.gz")
        ),
        quote = F,
        sep = "\t",
        na = "NA",
        row.names = F
      )
    
    ## save readme file
    readme <-
      file(here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, "_readme.txt")
      ))
    writeLines(c(
      paste0("Date: ", Sys.Date()),
      "- genomeAssembly: GRCh37/hg19",
      "- traitDescription: 
        All cases met the criteria for Alzheimer's Disease according to the 
        NINCDS-ADRDA, with biomarker evidence of amyloid pathology",
      paste0("- sampleSize: ", n_case + n_control),
      paste0("- caseCount: ", n_case),
      paste0("- controlCount: ", n_control),
      "- caseControlStudy: Yes",
      "- sampleAncestry: European",
      "- genotypingTechnology: 
            Genotyping of the samples was conducted using Whole Genome Sequencing (WGS) 
            or genotyping arrays (Illumina Infinium Multi-Ethnic Global 8 and GSAv2.0)",
      "- analysisSoftware: REGENIE",
      "- imputationPanel: GNE reference panel with >25k multi-ethnic individuals",
      "- imputationSoftware: 
        We used Beagle 5.1 for imputation. The imputations were conducted with 
        hg38 HapMap genetic maps and a Genentech reference panel, which includes 
        over 25,000 multi-ethnic individuals.",
      paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
      "- ancestryMethod: 
        European ancestry was determined through an analysis using ADMIXTURE 1.3. 
        We set an ancestry threshold of 0.7 and used samples from Phase 3 of 
        the 1000 Genomes Project as reference, classifying them by their 
        superpopulation labels (AFR = African, AMR = AdMixed American, 
        EAS = East Asian, EUR = European)",
      "- sortedByGenomicLocation: Yes",
      "- dosageOrHardcalls: Hard calls",
      "- notes: ~75 SNPs failed with liftOver"
    ), 
    readme)
    close(readme)
  } else {
    
    ## print if code is skipped because file already exists
    print(paste0(outfile, ": already exists")) 
    
  }
  
}

## clear workspace
rm(list=setdiff(ls(), c("neff_f", "snp_modifyBuild_local")))


#### 09. FinnGen ####

## set file names
proxy_casecontrol <- "case_control"

infile <-
  here("data_raw",
       "cohorts",
       "finngen",
       "finngen_R10_G6_AD_WIDE.gz")

outfile <-
  paste0("finng_", proxy_casecontrol, "_main_eur")

## sample size
manifest <-
  fread(here(
    "data_raw",
    "cohorts",
    "finngen",
    "R10_manifest.tsv"
  ))
n_case <- manifest[manifest$phenocode == "G6_AD_WIDE",]$num_cases
n_control <- manifest[manifest$phenocode == "G6_AD_WIDE",]$num_controls

## only run if outfile does not exist
if (!file.exists(here(
  "data_raw",
  "sumstats",
  proxy_casecontrol,
  paste0(outfile, ".txt.gz")
))) {
  
  
  g1000_eur <- fread(here("data", "reference_data", "g1000_eur", "g1000_eur.frq"))
  
  ## load & format sumstats
  print(paste0("load & format sumstats: ", infile))
  fread(infile) %>% 
    mutate(
      variant_id_hg38 = paste(
        `#chrom`, 
        pos, 
        ref, 
        alt, 
        sep = ":"),
      n_case = n_case,
      n_control = n_control
    ) %>%
    # liftOver requires chr23 to be coded as X
    mutate(`#chrom` = ifelse(`#chrom` == 23, "X", `#chrom`)) %>%
    rename(chr = `#chrom`,
           pos = pos) %>%
    snp_modifyBuild_local(liftOver = here("src", "liftOver"),
                    from = "hg38",
                    to = "hg19",
                    check_reverse = T) %>%
    mutate(chr = as.integer(ifelse(chr == "X", 23, chr))) %>%
    rename(
      chromosome = chr,
      base_pair_location = pos,
      other_allele = ref,
      effect_allele = alt,
      beta = beta,
      standard_error = sebeta,
      effect_allele_frequency = af_alt,
      p_value = pval,
      rsid = rsids
    ) %>%
    # get base_pair_location for SNPs where liftOver failed
    left_join(g1000_eur, by = c("rsid" = "SNP")) %>%
    mutate(base_pair_location = ifelse(is.na(base_pair_location), POS, base_pair_location),
           chromosome = ifelse(is.na(chromosome), CHR, chromosome)) %>%
    mutate(
      variant_id = paste(
        chromosome, 
        base_pair_location, 
        other_allele, 
        effect_allele, 
        sep = ":"),
      neff = neff_f(n_case = n_case,
                    n_control = n_control,
                    proxy_casecontrol = proxy_casecontrol)
    ) %>%
    select(
      variant_id,
      chromosome,
      base_pair_location,
      other_allele,
      effect_allele,
      beta,
      standard_error,
      effect_allele_frequency,
      p_value,
      neff,
      n_case,
      n_control,
      rsid,
      variant_id_hg38
    ) %>%
    arrange(chromosome, base_pair_location) %>%
    {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
    fwrite(
      file = here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, ".txt.gz")
      ),
      quote = F,
      sep = "\t",
      na = "NA",
      row.names = F
    )
  
  ## save readme file
  readme <-
    file(here(
      "data_raw",
      "sumstats",
      proxy_casecontrol,
      paste0(outfile, "_readme.txt")
    ))
  writeLines(c(
    paste0("Date: ", Sys.Date()),
    "- genomeAssembly: GRCh37/hg19",
    "- traitDescription: Alzheimer’s disease, wide definition",
    paste0("- sampleSize: ", n_case + n_control),
    paste0("- caseCount: ", n_case),
    paste0("- controlCount: ", n_control),
    "- caseControlStudy: Yes",
    "- sampleAncestry: European (Finnish)",
    "- genotypingTechnology: Illumina and Affymetrix chip arrays",
    "- analysisSoftware: REGENIE",
    "- imputationPanel: 
          Chip genotype data were imputed using the population-specific 
          SISu v4.2 imputation reference panel of 8,554 whole genomes.",
    "- imputationSoftware: Beagle 4.1 (version 27Jan18.7e1)",
    paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
    "- ancestryMethod: 
          FinnGen data was merged with the 1k genome project (1kgp) data.
          A round of PCA was performed and a bayesian algorithm was used to 
          spot outliers. This process got rid of 14,547 FinnGen samples. 
          See for more details https://finngen.gitbook.io/documentation/methods/phewas/quality-checks",
    "- sortedByGenomicLocation: Yes",
    "- dosageOrHardcalls: ?",
    "- notes: 
          info column not available but this is listed under 'Sample QC and PCA' 
          (https://finngen.gitbook.io/documentation/methods/phewas/quality-checks): 
          only variants with a minimum info score of 0.9 in all batches are kept.
    
          169,070 variants have not been mapped with liftOver (mostly on chr23, which
          is also coded for chrY in finngen)"
  ), 
  readme)
  close(readme)
} else {
  
  ## print if code is skipped because file already exists
  print(paste0(outfile, ": already exists")) 
  
}

## clear workspace
rm(list=setdiff(ls(), c("neff_f", "snp_modifyBuild_local")))


#### 10. deCODE ####

## set file names
proxy_casecontrol <- "case_control"

infile <- here(
  "data_raw",
  "cohorts",
  "decode",
  "deCODE_Alzheimers_ICD10_F00_G30_sumstat_july2024",
  "deCODE_Alzheimers_ICD10_F00_G30_summarystats_july2024.txt.gz"
)

outfile <- paste0("decod_", proxy_casecontrol, "_main_eur")

# sample size reported "deCODE_Alzheimers_ICD10_F00_G30_summarystats_info_july2024.txt"
n_case <- 9559
n_control <- 256998

## only run if outfile does not exist
if (!file.exists(here(
  "data_raw",
  "sumstats",
  proxy_casecontrol,
  paste0(outfile, ".txt.gz")
))) {
  
  ## load & format sumstats
  print(paste0("load & format sumstats: ", infile))
  fread(infile) %>% 
    mutate(
      Chr = ifelse(Chr == "chrX", "chr23", Chr)
    ) %>%
    mutate(
      Chr = as.integer(str_remove(Chr, "chr")),
      n_case = !!n_case,
      n_control = !!n_control,
      EAFrq = EAFrq / 100
      ) %>%
    mutate(
      variant_id_hg38 = paste(Chr, Pos, OA, EA, sep = ":")
    ) %>%
    ## change column names for liftover
    rename(chr = Chr,
           pos = Pos) %>%
    snp_modifyBuild_local(
      liftOver = here("src", "liftOver"),
      from = "hg38",
      to = "hg19",
      check_reverse = T
    ) %>% 
    rename(
      chromosome = chr,
      base_pair_location = pos,
      other_allele = OA,
      effect_allele = EA, # minor allele
      beta = Beta,
      standard_error = SE,
      effect_allele_frequency = EAFrq,
      p_value = P,
      rsid = rsID
    ) %>% 
    mutate(
      variant_id = paste(
        chromosome, 
        base_pair_location, 
        other_allele, 
        effect_allele, 
        sep = ":"),
      neff = neff_f(
        n_case = n_case, 
        n_control = n_control,
        proxy_casecontrol = !!proxy_casecontrol)
    ) %>%
    select(
      variant_id,
      chromosome,
      base_pair_location,
      other_allele,
      effect_allele,
      beta,
      standard_error,
      effect_allele_frequency,
      p_value,
      info,
      neff,
      n_case,
      n_control,
      variant_id_hg38,
      rsid
    ) %>%
    arrange(chromosome, base_pair_location) %>%
    {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
    fwrite(
      file = here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, ".txt.gz")
      ),
      quote = F,
      sep = "\t",
      na = "NA",
      row.names = F
    )
  
  ## save readme file
  readme <-
    file(here(
      "data_raw",
      "sumstats",
      proxy_casecontrol,
      paste0(outfile, "_readme.txt")
    ))
  writeLines(c(
    paste0("Date: ", Sys.Date()),
    "- genomeAssembly: GRCh37/hg19",
    "- traitDescription: AD in icelandic population",
    paste0("- sampleSize: ", n_case + n_control),
    paste0("- caseCount: ", n_case),
    paste0("- controlCount: ", n_control),
    "- caseControlStudy: Yes",
    "- sampleAncestry: European (Icelandish)",
    "- genotypingTechnology: ?",
    "- analysisSoftware: ?",
    "- imputationPanel: ?",
    "- imputationSoftware: ?",
    paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
    "- ancestryMethod: ?",
    "- sortedByGenomicLocation: Yes",
    "- dosageOrHardcalls: ?",
    "- covariates: ?",
    "- notes: 
      - 534,583 could not be mapped with liftOver"
  ), 
  readme)
  close(readme)
} else {
  
  ## print if code is skipped because file already exists
  print(paste0(outfile, ": already exists")) 
  
}

## clear workspace
rm(list=setdiff(ls(), c("neff_f", "snp_modifyBuild_local")))


#### 11. bioVU ####

## set file names
proxy_casecontrol <- "case_control"

infile <-
  here("data_raw",
       "cohorts",
       "biovu",
       "BioVU.txt.gz")

outfile <-
  paste0("biovu_", proxy_casecontrol, "_noag_eur")

# sample size reported in Wightman et al. (2021) and Alz_BioVU_readme.docx
n_case <- 600
n_control <- 36059

## only run if outfile does not exist
if (!file.exists(here(
  "data_raw",
  "sumstats",
  proxy_casecontrol,
  paste0(outfile, ".txt.gz")
))) {
  
  ## load & format sumstats
  print(paste0("load & format sumstats: ", infile))
  fread(infile) %>% 
    rename(
      chromosome = CHR,
      base_pair_location = POS_GRCh37,
      other_allele = Non_Effect_allele,
      effect_allele = Effect_allele,
      beta = BETA,
      standard_error = SE,
      p_value = P,
      rsid = SNP
      ) %>% 
    mutate(
      n_case = n_case,
      n_control = n_control
    ) %>%
    mutate(
      variant_id = paste(
        chromosome,
        base_pair_location,
        other_allele,
        effect_allele,
        sep = ":"
      ),
      neff = neff_f(n_case = n_case,
                    n_control = n_control,
                    proxy_casecontrol = proxy_casecontrol),
      effect_allele_frequency = (
        ((n_case / (n_case + n_control)) * Effect_allele_freq_cases) + 
          ((n_control / (n_case + n_control)) * Effect_allele_freq_controls)
        )
      ) %>% 
    select(
      variant_id,
      chromosome,
      base_pair_location,
      other_allele,
      effect_allele,
      beta,
      standard_error,
      effect_allele_frequency,
      p_value,
      neff,
      n_case,
      n_control,
      rsid
    ) %>%
    arrange(chromosome, base_pair_location) %>%
    {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
    fwrite(
      file = here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, ".txt.gz")
      ),
      quote = F,
      sep = "\t",
      na = "NA",
      row.names = F
    )
  
  ## save readme file
  readme <-
    file(here(
      "data_raw",
      "sumstats",
      proxy_casecontrol,
      paste0(outfile, "_readme.txt")
    ))
  writeLines(c(
    paste0("Date: ", Sys.Date()),
    "- genomeAssembly: GRCh37/hg19",
    "- traitDescription: 
          Cases were defined as individuals diagnosed with ICD-10 G30 and 
          ICD-9 331.0. Controls were individuals without any of the following 
          ICD-10 diagnoses; G30, F01, F02, F03, F10.27, F10.97, F13.27, F13.97, 
          F18.17, F18.27, F18.97, F19.17, F19.27, F19.97, G31.0, G31.83 and the 
          following ICD-9 diagnoses; 331.0 ,290, 291.2, 292.82, 294.1, 294.10, 
          294.11, 294.2, 294.20, 294.21, 331.19, 331.82. Individuals with a 
          family history of dementia in their electronic health records were 
          also excluded from the control sample.",
    paste0("- sampleSize: ", n_case + n_control),
    paste0("- caseCount: ", n_case),
    paste0("- controlCount: ", n_control),
    "- caseControlStudy: Yes",
    "- sampleAncestry: European",
    "- genotypingTechnology: Illumina MEGAEX array",
    "- analysisSoftware: SAIGE",
    "- imputationPanel: HRC panel",
    "- imputationSoftware: Michigan Imputation Server",
    paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
    "- ancestryMethod: 
          Principal component analysis (PCA) was used to determine BioVU 
          individuals of European genetic ancestry. First, we performed PCA 
          using FlashPCA55 on BioVU combined with CEU, YRI, and CHB reference 
          sets from 1000 Genomes Project Phase 343. Principal components were 
          scaled so that the axes could be interpreted as proportions of 
          genetic ancestry. We selected BioVU individuals who were within 40% 
          of the CEU cluster along the CEU-CHB axis and within 30% of the CEU 
          cluster on the CEU-YRI axis, generating a oncePCA filtered European 
          set. To ensure subsequent steps would remove SNPs associated with 
          reduced quality rather than cryptic population substructure, we 
          filtered the previously identified BioVU European cluster to identify 
          individuals falling within the CEU, TSI, and GIH 1000 genomes 
          populations, producing a twice-filtered European set",
    "- sortedByGenomicLocation: Yes",
    "- dosageOrHardcalls: Hard calls",
    "- covariates: Sex and top 10 PCs",
    "- notes: 
          Sumstats do not include 'info' column, but readme states they filtered
          for SNPs with 'info' > 0.3. A very loose threshold, considering they
          used hard calls and not dosages." 
  ), 
  readme)
  close(readme)
} else {
  
  ## print if code is skipped because file already exists
  print(paste0(outfile, ": already exists")) 
  
}

## clear workspace
rm(list=setdiff(ls(), c("neff_f", "snp_modifyBuild_local")))



#### 12. Gothenburg ####

## function to load, format, process, save, and plot jpb-sumstats
gothenburg_sumstats <- function(analysis, proxy_casecontrol) {
  
  ## set file names
  infile <-
    Sys.glob(here(
      "data_raw/cohorts/gothenburg/sumstats",
      paste0("*_", analysis, ".txt.gz")
    ))
  
  if (analysis == "sex") {analysis <- "noag"}
  if (analysis == "sex_APOE4") {analysis <- "apoe"}
  if (analysis == "males") {analysis <- "male"}
  if (analysis == "females") {analysis <- "fema"}
  
  outfile <- paste0("gothe_", proxy_casecontrol, "_", analysis, "_eur")
  
  ## only run if outfile does not exist
  if (!file.exists(here(
    "data_raw",
    "sumstats",
    proxy_casecontrol,
    paste0(outfile, ".txt.gz")
  ))) {
    
    ## load & format sumstats
    print(paste0("load, liftover (hg38 to GRCh37/hg19) & format sumstats: ", infile))
    fread(infile) %>% 
      mutate(
        variant_id_hg38 = paste(
          chromosome, 
          base_pair_location, 
          other_allele, 
          effect_allele, 
          sep = ":")
      ) %>%
      ## change column names for liftover
      rename(chr = chromosome,
             pos = base_pair_location) %>%
      snp_modifyBuild_local(liftOver = here("src", "liftOver"),
                            from = "hg38",
                            to = "hg19",
                            check_reverse = T) %>%
      rename(chromosome = chr,
             base_pair_location = pos) %>%
      {max(.$n_case) ->> n_case;.} %>%
      {max(.$n_control) ->> n_control;.} %>%
      mutate(
        neff = neff_f(
          n_case = n_case,
          n_control = n_control,
          proxy_casecontrol = !!proxy_casecontrol)
      ) %>% 
      mutate(
        variant_id = paste(
          chromosome, 
          base_pair_location, 
          other_allele, 
          effect_allele, 
          sep = ":")
      ) %>%
      select(
        variant_id,
        chromosome,
        base_pair_location,
        other_allele,
        effect_allele,
        beta,
        standard_error,
        effect_allele_frequency,
        p_value,
        info,
        neff,
        n_case,
        n_control,
        variant_id_hg38
      ) %>%
      arrange(chromosome, base_pair_location) %>%
      {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
      fwrite(
        file = here(
          "data_raw",
          "sumstats",
          proxy_casecontrol,
          paste0(outfile, ".txt.gz")
        ),
        quote = F,
        sep = "\t",
        na = "NA",
        row.names = F
      )
    
    print(paste0(
      "saved in: ",
      here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, ".txt.gz")
      )
    ))
    
    if (analysis == "noag") {covars <- "sex,PCs1-10"}
    if (analysis == "apoe") {covars <- "sex,APOE4,PC1-10"}
    if (analysis == "male") {covars <- "PCs1-10"}
    if (analysis == "fema") {covars <- "PCs1-10"}
    
    ## save readme file
    readme <-
      file(here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, "_readme.txt")
      ))
    writeLines(c(
      paste0("Date: ", Sys.Date()),
      "- genomeAssembly: GRCh37/hg19",
      "- traitDescription: 
          The Gothenburg H70 Birth Cohort Studies and Clinical AD from Sweden 
          (Gothenburg) AD cases originate from Sweden and were either collected 
          in memory clinics (in different parts of Sweden) or as a part of two 
          population-based epidemiological studies in Gothenburg; the Prospective 
          Population Study of Women (PPSW) and the Gothenburg Birth Cohort 
          Studies (H70, H85 and 95+), described in detail previously 57–60. 
          Controls originate from the Gothenburg Birth Cohort Studies and PPSW.",
      paste0("- sampleSize: ", n_case + n_control),
      paste0("- caseCount: ", n_case),
      paste0("- controlCount: ", n_control),
      "- caseControlStudy: Yes",
      "- sampleAncestry: European (Swedish)",
      "- genotypingTechnology: Illumina Neurochip array",
      "- analysisSoftware: ?",
      "- imputationPanel: ?",
      "- imputationSoftware: ?",
      paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
      "- ancestryMethod: 
          Individuals of non-European descent were excluded as part of the QC 
          of the GWAS-data",
      "- sortedByGenomicLocation: Yes",
      "- dosageOrHardcalls: ?",
      paste0("- covariates: ", covars),
      "- notes: " 
    ), 
    readme)
    close(readme)
  } else {
    
    ## print if code is skipped because file already exists
    print(paste0(outfile, ": already exists")) 
    
  }
}

## parameters defining different sumstats
params <- expand.grid(
  analysis = c("males", "females", "sex_APOE4", "sex"),
  proxy_casecontrol = "case_control",
  stringsAsFactors = F
)

params %>%
  pmap(~ gothenburg_sumstats(
    analysis = ..1,
    proxy_casecontrol = ..2
  ))

## clear workspace
rm(list=setdiff(ls(), c("neff_f", "snp_modifyBuild_local")))


#### 13. HUNT ####

## function to load, format, process, save, and plot jpb-sumstats
hunt_sumstats <- function(proxy_casecontrol, analysis, analysis_nr, ancestry) { 
  
  infile <- here(
    "data_raw/cohorts/hunt/2024_08_28",
    paste0("HUNT.MORENO.CASE_CONTROL_GWAS_AD_DEF", analysis_nr, ".txt.gz")
    )
  
  outfile <- paste("xhunt", proxy_casecontrol, analysis, "eur", sep = "_")
  
  ## only run if outfile does not exist
  if (!file.exists(here(
    "data_raw",
    "sumstats",
    proxy_casecontrol,
    paste0(outfile, ".txt.gz")
  ))) {
    
    ## load & format sumstats
    print(paste0("load & format sumstats: ", infile))
    fread(infile) %>% 
      ## delimiters between columns differ, hence the following extraction
      mutate(
        rsid = stringr::str_extract(`rsid n_control n_case`, "^[^ ]+"),  # Extract the first part before the first space
        n_control = stringr::str_extract(`rsid n_control n_case`, "(?<= )[0-9]+(?= )"),  # Extract the middle part (control)
        n_case = stringr::str_extract(`rsid n_control n_case`, "(?<= )[0-9]+$")  # Extract the last part (case)
      ) %>% 
      mutate(
        n_control = as.integer(n_control),
        n_case = as.integer(n_case)     
        ) %>%
      select(-`rsid n_control n_case`) %>%
      rename(
        standard_error = SEBETA,
        p_value = PVALUE,
        info = INFO      
        ) %>% 
      mutate(
        variant_id = paste(
          chromosome,
          base_pair_location,
          other_allele,
          effect_allele,
          sep = ":"
        ),
        neff = neff_f(
          n_case = n_case,
          n_control = n_control,
          proxy_casecontrol = proxy_casecontrol
          )
      ) %>% 
      {max(.$n_control, na.rm = T) ->> n_control;.} %>%
      {max(.$n_case, na.rm = T) ->> n_case;.} %>%
      select(
        variant_id,
        chromosome,
        base_pair_location,
        other_allele,
        effect_allele,
        beta,
        standard_error,
        effect_allele_frequency,
        p_value,
        info,
        neff,
        n_case,
        n_control,
        rsid
      ) %>%
      arrange(chromosome, base_pair_location) %>%
      {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
      fwrite(
        file = here(
          "data_raw",
          "sumstats",
          proxy_casecontrol,
          paste0(outfile, ".txt.gz")
        ),
        quote = F,
        sep = "\t",
        na = "NA",
        row.names = F
      )
    print(paste0(
      "saved in: ",
      here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, ".txt.gz")
      )
    ))
    
    ## save readme file
    readme <-
      file(here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, "_readme.txt")
      ))
    writeLines(c(
      paste0("Date: ", Sys.Date()),
      "- genomeAssembly: GRCh37/hg19",
      "- traitDescription: 
          followed analysis plan v2.6",
      paste0("- sampleSize: ", n_case + n_control),
      paste0("- caseCount: ", n_case),
      paste0("- controlCount: ", n_control),
      "- caseControlStudy: Yes",
      "- sampleAncestry: European (Norwegian)",
      "- genotypingTechnology: 
          Illumina HumanCoreExome arrays 
          (HumanCoreExome12 v1.0, HumanCoreExome12 v1.1, or UM HUNT Biobank v1.0)",
      "- analysisSoftware: (SAIGE)",
      "- imputationPanel: HRC",
      "- imputationSoftware: - ",
      paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
      "- ancestryMethod: 
          -",
      "- sortedByGenomicLocation: Yes",
      "- dosageOrHardcalls: dosages",
      "- covariates: analysis nr. + batch ",
      "- notes: 
          -" 
    ), 
    readme)
    close(readme)
  } else {
    
    ## print if code is skipped because file already exists
    print(paste0(outfile, ": already exists")) 
    
  }
}

## parameters defining different sumstats
case_control_params <- data.table(
  proxy_casecontrol = rep("case_control", 5),
  analysis = c("main", "male", "fema", "noag", "apoe"),
  analysis_nr = c(1, 2, 3, 4, 5),
  ancestry = rep("eur", 5)
)

case_control_params %>%
  pmap(~ hunt_sumstats(
    proxy_casecontrol = ..1,
    analysis = ..2,
    analysis_nr = ..3,
    ancestry = ..4
  ))



## clear workspace
rm(list=setdiff(ls(), c("neff_f", "snp_modifyBuild_local")))



#### 14. Estonian biobank ####

## set file names
proxy_casecontrol <- "case_control"

infile <-
  here("data_raw",
       "cohorts",
       "estonia",
       "alzheimer_all_EstBB_GWAS_results.txt.gz")

outfile <-
  paste0("eston_", proxy_casecontrol, "_main_eur")

# sample size reported in README.txt
n_case <- 533
n_control <- 194928

## only run if outfile does not exist
if (!file.exists(here(
  "data_raw",
  "sumstats",
  proxy_casecontrol,
  paste0(outfile, ".txt.gz")
))) {
  
  ## load & format sumstats
  print(paste0("load & format sumstats: ", infile))
  fread(infile) %>% 
    rename(
      chromosome = CHR,
      base_pair_location = POS,
      other_allele = Allele1,
      effect_allele = Allele2,
      beta = BETA,
      standard_error = SE,
      p_value = p.value    
    ) %>% 
    mutate(
      n_case = n_case,
      n_control = n_control
    ) %>%
    mutate(
      variant_id = paste(
        chromosome,
        base_pair_location,
        other_allele,
        effect_allele,
        sep = ":"
      ),
      neff = neff_f(n_case = n_case,
                    n_control = n_control,
                    proxy_casecontrol = proxy_casecontrol),
      effect_allele_frequency = (
        ((n_case / (n_case + n_control)) * AF.Cases) + 
          ((n_control / (n_case + n_control)) * AF.Controls)
      )
    ) %>% 
    select(
      variant_id,
      chromosome,
      base_pair_location,
      other_allele,
      effect_allele,
      beta,
      standard_error,
      effect_allele_frequency,
      p_value,
      neff,
      n_case,
      n_control
    ) %>%
    arrange(chromosome, base_pair_location) %>%
    {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
    fwrite(
      file = here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, ".txt.gz")
      ),
      quote = F,
      sep = "\t",
      na = "NA",
      row.names = F
    )
  
  print(paste0(
    "saved in: ",
    here(
      "data_raw",
      "sumstats",
      proxy_casecontrol,
      paste0(outfile, ".txt.gz")
    )
  ))
  
  ## save readme file
  readme <-
    file(here(
      "data_raw",
      "sumstats",
      proxy_casecontrol,
      paste0(outfile, "_readme.txt")
    ))
  writeLines(c(
    paste0("Date: ", Sys.Date()),
    "- genomeAssembly: GRCh37/hg19",
    "- traitDescription: 
          Cases of Alzheimer's disease are individuals with F00* or G30* diagnosis.
          All individuals with AD diagnosis before the age of 30 were removed among cases (n=25)
          Controls are all individuals without the following ICD-10 diagnoses: 
          F01*,F02*,F03*,F05.1,F10.6,F10.73,G31*,I67.3",
    paste0("- sampleSize: ", n_case + n_control),
    paste0("- caseCount: ", n_case),
    paste0("- controlCount: ", n_control),
    "- caseControlStudy: Yes",
    "- sampleAncestry: European (Estonian)",
    "- genotypingTechnology: ?",
    "- analysisSoftware: SAIGE",
    "- imputationPanel: ?",
    "- imputationSoftware: ?",
    paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
    "- ancestryMethod: ?",
    "- sortedByGenomicLocation: Yes",
    "- dosageOrHardcalls: ?",
    "- covariates: Gender, birth year, 10PC",
    "- notes: info column is all NA" 
  ), 
  readme)
  close(readme)
} else {
  
  ## print if code is skipped because file already exists
  print(paste0(outfile, ": already exists")) 
  
}

## clear workspace
rm(list=setdiff(ls(), c("neff_f", "snp_modifyBuild_local")))



#### 15. GR@CE ####

## set file names
proxy_casecontrol <- "case_control"

infile <-
  here("data_raw",
       "cohorts",
       "grace",
       "GRACE_StageI.txt")

outfile <-
  paste0("grace_", proxy_casecontrol, "_noag_eur")

# sample size reported in readme.txt
n_case <- 4120
n_control <- 3289

## only run if outfile does not exist
if (!file.exists(here(
  "data_raw",
  "sumstats",
  proxy_casecontrol,
  paste0(outfile, ".txt.gz")
))) {
  
  ## load & format sumstats
  print(paste0("load & format sumstats: ", infile))
  fread(infile) %>% 
    rename(
      chromosome = CHR,
      base_pair_location = BP,
      other_allele = NonEffect_Allele,
      effect_allele = Effect_Allele,
      odds_ratio = OR,
      standard_error = SE,
      p_value = P,
      effect_allele_frequency = FREQ_Effect_Allele,
      rsid = rsID
    ) %>% 
    mutate(
      n_case = n_case,
      n_control = n_control
    ) %>%
    mutate(
      variant_id = paste(
        chromosome,
        base_pair_location,
        other_allele,
        effect_allele,
        sep = ":"
      ),
      neff = neff_f(n_case = n_case,
                    n_control = n_control,
                    proxy_casecontrol = proxy_casecontrol),
      beta = log(odds_ratio)
    ) %>% 
    select(
      variant_id,
      chromosome,
      base_pair_location,
      other_allele,
      effect_allele,
      beta,
      odds_ratio,
      standard_error,
      effect_allele_frequency,
      p_value,
      neff,
      n_case,
      n_control,
      rsid
    ) %>%
    arrange(chromosome, base_pair_location) %>%
    {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
    fwrite(
      file = here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, ".txt.gz")
      ),
      quote = F,
      sep = "\t",
      na = "NA",
      row.names = F
    )
  
  print(paste0(
    "saved in: ",
    here(
      "data_raw",
      "sumstats",
      proxy_casecontrol,
      paste0(outfile, ".txt.gz")
    )
  ))
  
  ## save readme file
  readme <-
    file(here(
      "data_raw",
      "sumstats",
      proxy_casecontrol,
      paste0(outfile, "_readme.txt")
    ))
  writeLines(c(
    paste0("Date: ", Sys.Date()),
    "- genomeAssembly: GRCh37/hg19",
    "- traitDescription: 
          All AD patients included in the GR@ACE study received a thorough 
          structured neurological evaluation that included: history, examination, 
          Mini-Mental State Examination (MMSE), Blessed Dementia Rating Scale (BDRS),
          Neuropsychiatric Inventory-questionnaire (NPI-Q), Tinnetti scale for 
          gait and balance, Clinical Dementia Rating (CDR), Global Deterioration Scale (GDR)
          scoring and Hachinski Ischemia Scale. The Fundacio ACE neuropsychological 
          battery (NBACE) was administered to all patients. The NBACE includes 
          measures of cognitive information processing speed, orientation, attention, 
          verbal learning and memory, language, visuoperception, praxis and 
          executive functions. Family members or caregivers are interviewed by a 
          social worker. 
          
          See Supplementary Methods in Moreno-Grau et al. (2019)",
    paste0("- sampleSize: ", n_case + n_control),
    paste0("- caseCount: ", n_case),
    paste0("- controlCount: ", n_control),
    "- caseControlStudy: Yes",
    "- sampleAncestry: European (Spanish)",
    "- genotypingTechnology: Axiom 815K Spanish Biobank array",
    "- analysisSoftware: PLINK 1.9",
    "- imputationPanel: Haplotype Reference Consortium panel",
    "- imputationSoftware: Michigan Imputation Server",
    paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
    "- ancestryMethod: 
          To detect population outliers of non-European ancestry (.6 SD from 
          European popula- tion mean), principal component analysis was conducted 
          using SMARTPCA from EIGENSOFT 6.1.4",
    "- sortedByGenomicLocation: Yes",
    "- dosageOrHardcalls: Dosage",
    "- covariates: 4 PCs",
    "- notes: 
          - Info column is not available and covariates do not include age and sex
          - Sumstats are the Stage I (only GR@CE subjects) results
          - Paper: Moreno-Grau et al. (2019). Genome-wide association analysis of 
                   dementia and its clinical endophenotypes reveal novel loci 
                   associated with Alzheimer’s disease and three causality 
                   networks: The GR@ACE project. Alzheimer’s & Dementia, 15(10), 
                   1333–1347. https://doi.org/10.1016/j.jalz.2019.06.4950
          - Downloaded from GWAS Catalogue: https://www.ebi.ac.uk/gwas/publications/31473137" 
  ), 
  readme)
  close(readme)
} else {
  
  ## print if code is skipped because file already exists
  print(paste0(outfile, ": already exists")) 
  
}

## clear workspace
rm(list=setdiff(ls(), c("neff_f", "snp_modifyBuild_local")))



#### 16. STSA ####

## function to load, format, process, save, and plot jpb-sumstats
stsa_sumstats <- function(analysis, proxy_casecontrol) {
  
  ## set file names
  infile <-
    Sys.glob(here(
      "data_raw",
      "cohorts",
      "stsa",
      paste0("*_", analysis, ".txt.gz")
    ))
  
  if (analysis == "sex") {analysis <- "noag"}
  if (analysis == "sex_APOE4") {analysis <- "apoe"}
  if (analysis == "males") {analysis <- "male"}
  if (analysis == "females") {analysis <- "fema"}
  
  outfile <-
    paste0("xstsa_", proxy_casecontrol, "_", analysis, "_eur")
  
  ## only run if outfile does not exist
  if (!file.exists(here(
    "data_raw",
    "sumstats",
    proxy_casecontrol,
    paste0(outfile, ".txt.gz")
  ))) {
    
    ## load & format sumstats
    print(paste0("load, liftover (hg38 to GRCh37/hg19) & format sumstats: ", infile))
    fread(infile) %>% 
      mutate(
        variant_id_hg38 = paste(
          chromosome, 
          base_pair_location, 
          other_allele, 
          effect_allele, 
          sep = ":")
      ) %>%
      ## change column names for liftover
      rename(chr = chromosome,
             pos = base_pair_location) %>%
      snp_modifyBuild_local(liftOver = here("src", "liftOver"),
                            from = "hg38",
                            to = "hg19",
                            check_reverse = T) %>%
      rename(chromosome = chr,
             base_pair_location = pos) %>%
      {max(.$n_case) ->> n_case;.} %>%
      {max(.$n_control) ->> n_control;.} %>%
      mutate(
        neff = neff_f(n_case = n_case,
                      n_control = n_control,
                      proxy_casecontrol = proxy_casecontrol)
      ) %>% 
      mutate(
        variant_id = paste(
          chromosome, 
          base_pair_location, 
          other_allele, 
          effect_allele, 
          sep = ":")
      ) %>%
      select(
        variant_id,
        chromosome,
        base_pair_location,
        other_allele,
        effect_allele,
        beta,
        standard_error,
        effect_allele_frequency,
        p_value,
        info,
        neff,
        n_case,
        n_control,
        variant_id_hg38
      ) %>%
      arrange(chromosome, base_pair_location) %>%
      {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
      fwrite(
        file = here(
          "data_raw",
          "sumstats",
          proxy_casecontrol,
          paste0(outfile, ".txt.gz")
        ),
        quote = F,
        sep = "\t",
        na = "NA",
        row.names = F
      )
    
    print(paste0(
      "saved in: ",
      here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, ".txt.gz")
      )
    ))
    
    if (analysis == "noag") {covars <- "sex,PCs1-10"}
    if (analysis == "apoe") {covars <- "sex,APOE4,PC1-10"}
    if (analysis == "male") {covars <- "PCs1-10"}
    if (analysis == "fema") {covars <- "PCs1-10"}
    
    ## save readme file
    readme <-
      file(here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, "_readme.txt")
      ))
    writeLines(c(
      paste0("Date: ", Sys.Date()),
      "- genomeAssembly: GRCh37/hg19",
      "- traitDescription: 
          Alzheimer’s disease patients were diagnosed as part of the studies 
          according to the NINCDS/ADRDA criteria. In addition, information on 
          disease after last study participation was retrieved from three 
          population- based health care registers: The National Patient Register, 
          the Causes of Death Register, and the Prescribed Drug Register.",
      paste0("- sampleSize: ", n_case + n_control),
      paste0("- caseCount: ", n_case),
      paste0("- controlCount: ", n_control),
      "- caseControlStudy: Yes",
      "- sampleAncestry: European (Swedish)",
      "- genotypingTechnology: Illumina Infinium PsychArray (see jansen 2019 supplement for more details)",
      "- analysisSoftware: regenie v3.4.1",
      "- imputationPanel: TOPMed",
      "- imputationSoftware: -",
      paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
      "- ancestryMethod: 
          non-European individuals within the datasets were removed based on 
          principle component analysis (PCA), using the 1KG Phase 3 dataset as a 
          reference. The PCA pipeline was repeated including all European 
          individuals in all genotype level datasets to identify individuals 
          across the datasets with a pihat > 0.2 for exclusion from the analysis.",
      "- sortedByGenomicLocation: Yes",
      "- dosageOrHardcalls: dosages",
      paste0("- covariates: ", covars),
      "- notes: 
          - Age was not available as a covariate, but cases and controls were matched.
            Info in Supplement of Wightman et al. (2021)" 
    ), 
    readme)
    close(readme)
  } else {
    
    ## print if code is skipped because file already exists
    print(paste0(outfile, ": already exists")) 
    
  }
}

## parameters defining different sumstats
params <- expand.grid(
  analysis = c("males", "females", "sex_APOE4", "sex"),
  proxy_casecontrol = "case_control",
  stringsAsFactors = F
)

params %>%
  pmap(~ stsa_sumstats(
    analysis = ..1,
    proxy_casecontrol = ..2
  ))

## clear workspace
rm(list=setdiff(ls(), c("neff_f", "snp_modifyBuild_local")))


#### 17. TwinGene ####

## function to load, format, process, save, and plot jpb-sumstats
twingene_sumstats <- function(analysis, proxy_casecontrol) {
  
  ## set file names
  infile <-
    Sys.glob(here(
      "data_raw/cohorts/twingene/sumstats",
      paste0("*_", analysis, ".txt.gz")
    ))
  
  if (analysis == "sex") {analysis <- "noag"}
  if (analysis == "sex_APOE4") {analysis <- "apoe"}
  if (analysis == "males") {analysis <- "male"}
  if (analysis == "females") {analysis <- "fema"}
  
  outfile <-
    paste0("twing_", proxy_casecontrol, "_", analysis, "_eur")
  
  ## only run if outfile does not exist
  if (!file.exists(here(
    "data_raw",
    "sumstats",
    proxy_casecontrol,
    paste0(outfile, ".txt.gz")
  ))) {
    
    ## load & format sumstats
    print(paste0("load, liftover (hg38 to GRCh37/hg19) & format sumstats: ", infile))
    fread(infile) %>% 
      mutate(
        variant_id_hg38 = paste(
          chromosome, 
          base_pair_location, 
          other_allele, 
          effect_allele, 
          sep = ":")
      ) %>%
      ## change column names for liftover
      rename(chr = chromosome,
             pos = base_pair_location) %>%
      snp_modifyBuild_local(liftOver = here("src", "liftOver"),
                            from = "hg38",
                            to = "hg19",
                            check_reverse = T) %>%
      rename(chromosome = chr,
             base_pair_location = pos) %>%
      {max(.$n_case) ->> n_case;.} %>%
      {max(.$n_control) ->> n_control;.} %>%
      mutate(
        neff = neff_f(
          n_case = n_case,
          n_control = n_control,
          proxy_casecontrol = !!proxy_casecontrol)
      ) %>% 
      mutate(
        variant_id = paste(
          chromosome, 
          base_pair_location, 
          other_allele, 
          effect_allele, 
          sep = ":")
      ) %>%
      select(
        variant_id,
        chromosome,
        base_pair_location,
        other_allele,
        effect_allele,
        beta,
        standard_error,
        effect_allele_frequency,
        p_value,
        info,
        neff,
        n_case,
        n_control,
        variant_id_hg38
      ) %>%
      arrange(chromosome, base_pair_location) %>%
      {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
      fwrite(
        file = here(
          "data_raw",
          "sumstats",
          proxy_casecontrol,
          paste0(outfile, ".txt.gz")
        ),
        quote = F,
        sep = "\t",
        na = "NA",
        row.names = F
      )
    
    print(paste0(
      "saved in: ",
      here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, ".txt.gz")
      )
    ))
    
    if (analysis == "noag") {covars <- "sex,PCs1-10"}
    if (analysis == "apoe") {covars <- "sex,APOE4,PC1-10"}
    if (analysis == "male") {covars <- "PCs1-10"}
    if (analysis == "fema") {covars <- "PCs1-10"}
    
    ## save readme file
    readme <-
      file(here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, "_readme.txt")
      ))
    writeLines(c(
      paste0("Date: ", Sys.Date()),
      "- genomeAssembly: GRCh37/hg19",
      "- traitDescription: 
          TwinGene1 is a population-based study of older twins drawn from the 
          Swedish Twin  Registry. Written informed consent was obtained from 
          all participants and the study was  approved by the Regional Ethics 
          Board in Stockholm. DNA was extracted from blood  samples and 
          genotyped using Illumina Human OmniExpress for 1,791 individuals.  
          Information about Alzheimer’s disease (n cases = 343, n controls = 9070) 
          was extracted from the National Patient Register, the Causes of Death 
          Register, and the Prescribed Drug  Register, all of which are 
          population-based health care registers with nationwide coverage.",
      paste0("- sampleSize: ", n_case + n_control),
      paste0("- caseCount: ", n_case),
      paste0("- controlCount: ", n_control),
      "- caseControlStudy: Yes",
      "- sampleAncestry: European (Swedish)",
      "- genotypingTechnology: Illumina Human OmniExpress",
      "- analysisSoftware: regenie v3.4.1",
      "- imputationPanel: TOPMed",
      "- imputationSoftware: -",
      paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
      "- ancestryMethod: 
          -",
      "- sortedByGenomicLocation: Yes",
      "- dosageOrHardcalls: dosages",
      paste0("- covariates: ", covars),
      "- notes: 
          - Age was not available as a covariate, but cases and controls were matched.
            Info in Supplement of Wightman et al. (2021)" 
    ), 
    readme)
    close(readme)
  } else {
    
    ## print if code is skipped because file already exists
    print(paste0(outfile, ": already exists")) 
    
  }
}

## parameters defining different sumstats
params <- expand.grid(
  analysis = c("males", "females", "sex_APOE4", "sex"),
  proxy_casecontrol = "case_control",
  stringsAsFactors = F
)

params %>%
  pmap(~ twingene_sumstats(
    analysis = ..1,
    proxy_casecontrol = ..2
  ))

## clear workspace
rm(list=setdiff(ls(), c("neff_f", "snp_modifyBuild_local")))



#### 18. UKB ####

## function to load, format, process, save, and plot jpb-sumstats
ukb_sumstats <- function(proxy_casecontrol, analysis, analysis_out, ancestry) {
  
  infile <-
    paste0(here(
      "data_raw",
      "cohorts",
      "ukb",
      "PGCALZ3_UKB_For_Emil",
      paste0("MAF_",
             ancestry,
             "_",
             analysis,
             ".regenie.gz")
    ))
  
  outfile <- paste0("xukbb_", proxy_casecontrol, "_", analysis_out, "_", tolower(ancestry))
  
  ## only run if outfile does not exist
  if (!file.exists(here(
    "data_raw",
    "sumstats",
    proxy_casecontrol,
    paste0(outfile, ".txt.gz")
  ))) {
    
    ## load & format sumstats
    print(paste0("load & format sumstats: ", infile))
    fread(infile) %>% 
      {max(.$n_control, na.rm = T) ->> n_control;.} %>%
      {max(.$n_case, na.rm = T) ->> n_case;.} %>%
      mutate(
        neff = neff_f(n_case = n_case,
                      n_control = n_control,
                      proxy_casecontrol = proxy_casecontrol),
        variant_id = paste(
          chromosome, 
          base_pair_location, 
          other_allele, 
          effect_allele, 
          sep = ":")
      ) %>%
      select(
        variant_id,
        chromosome,
        base_pair_location,
        other_allele,
        effect_allele,
        beta,
        standard_error,
        effect_allele_frequency,
        p_value,
        info,
        neff,
        n_case,
        n_control
      ) %>%
      arrange(chromosome, base_pair_location) %>%
      {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
      fwrite(
        file = here(
          "data_raw",
          "sumstats",
          proxy_casecontrol,
          paste0(outfile, ".txt.gz")
        ),
        quote = F,
        sep = "\t",
        na = "NA",
        row.names = F
      )
    
    print(paste0(
      "saved in: ",
      here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, ".txt.gz")
      )
    ))
    
    if (proxy_casecontrol == "case_control") {
      caseControlStudy <- "yes"
    } else {
      caseControlStudy <- "no"
    }
    
    ## save readme file
    readme <-
      file(here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, "_readme.txt")
      ))
    writeLines(c(
      paste0("Date: ", Sys.Date()),
      "- genomeAssembly: GRCh37/hg19",
      "- traitDescription: 
          True case definition:
          Proxy case definition:",
      paste0("- sampleSize: ", n_case + n_control),
      paste0("- caseCount: ", n_case),
      paste0("- controlCount: ", n_control),
      paste0("- caseControlStudy: ", caseControlStudy),
      paste0("- sampleAncestry: ", tolower(ancestry)),
      "- genotypingTechnology: ",
      "- analysisSoftware: ",
      "- imputationPanel: -",
      "- imputationSoftware: -",
      paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
      "- ancestryMethod: PCA",
      "- dosageOrHardcalls: Hard calls",
      paste0("- covariates: see /home/emilu/projects/pgc_alzheimers/data_raw/cohorts/ukb/PGCALZ3_UKB_For_Emil/methods/PGCALZ3_UKB_methds.docx"),
      "- notes:  "
    ), 
    readme)
    close(readme)
  } else {
    
    ## print if code is skipped because file already exists
    print(paste0(outfile, ": already exists")) 
    
  }
}

## parameters defining different sumstats
params_casecontrol <- expand.grid(
  proxy_casecontrol = "case_control",
  analysis = c(
    "AD_age_age2_females_AD",
    "AD_age_age2_males_AD",
    "AD_sex__AD",
    "AD_sex_age_age2_agesex__AD",
    "AD_sex_age_age2_agesex_APOE4__AD"
  ),
  ancestry = c("EUR", "AFR", "SAS")
) %>%
  mutate(
    analysis_out =
      case_when(
        analysis == "AD_age_age2_females_AD" ~ "fema",
        analysis == "AD_age_age2_males_AD" ~ "male",
        analysis == "AD_sex__AD" ~ "noag",
        analysis == "AD_sex_age_age2_agesex__AD" ~ "main",
        analysis == "AD_sex_age_age2_agesex_APOE4__AD" ~ "apoe"
      )
  )

params_proxy <- expand.grid(
  proxy_casecontrol = "proxy",
  analysis = c(
    "Maternal_sex_age_age2_MotherAge_MotherAge2__Maternal",
    "Maternal_sex__Maternal",
    "Maternal_sex_age_age2_MotherAge_MotherAge2_APOE4__Maternal",
    "Paternal_sex_age_age2_FatherAge_FatherAge2__Paternal",
    "Paternal_sex__Paternal",
    "Paternal_sex_age_age2_FatherAge_FatherAge2_APOE4__Paternal"
  ),
  ancestry = c("EUR", "AFR", "SAS")
) %>%
  mutate(
    analysis_out =
      case_when(
        analysis == "Maternal_sex_age_age2_MotherAge_MotherAge2__Maternal" ~ "mother_main",
        analysis == "Maternal_sex__Maternal" ~ "mother_noag",
        analysis == "Maternal_sex_age_age2_MotherAge_MotherAge2_APOE4__Maternal" ~ "mother_apoe",
        analysis == "Paternal_sex_age_age2_FatherAge_FatherAge2__Paternal" ~ "father_main",
        analysis == "Paternal_sex__Paternal" ~ "father_noag",
        analysis == "Paternal_sex_age_age2_FatherAge_FatherAge2_APOE4__Paternal" ~ "father_apoe",
      )
  )

params <- rbind(
  params_casecontrol,
  params_proxy
)

params %>%
  pmap(~ ukb_sumstats(
    proxy_casecontrol = ..1,
    analysis = ..2,
    ancestry = ..3,
    analysis_out = ..4
  ))

## clear workspace
rm(list=setdiff(ls(), c("neff_f", "snp_modifyBuild_local")))



#### 19. Lifelines ####

## function to load, format, process, save, and plot jpb-sumstats
lifelines_sumstats <- function(analysis, cohort, parent, proxy_casecontrol, Ns) {
  
  ## set file names
  infile <-
    Sys.glob(here(
      "data_raw",
      "cohorts",
      "lifelines",
      "sumstats",
      paste0("*_", cohort, tolower(parent), "_", analysis, "_", parent, ".txt.gz")
    ))
  
  ## rename analysis
  covars_readme <- analysis
  if (analysis == "SexPCs") {analysis <- "noag"}
  if (analysis == "SexPCsMother_ageMother_age2AgeAge2APOE4" |
      analysis == "SexPCsFather_ageFather_age2AgeAge2APOE4") {
    analysis <- "apoe"
  }
  if (analysis == "SexPCsMother_ageMother_age2AgeAge2" |
      analysis == "SexPCsFather_ageFather_age2AgeAge2") {
    analysis <- "main"
  }
  
  ## rename parent
  if (parent == "Paternal") {parent <- "father"}
  if (parent == "Maternal") {parent <- "mother"}
  
  ## Genotyping array
  if (cohort == "GWAS") {
    geno_array <- "Illumina CytoSNP-12v2 array"
    name <- "life1"
    }
  if (cohort == "UGLI1") {
    geno_array <- "Infinium Global Screening Array (GSA) MultiEthnic Disease Version"
    name <- "life2"
    }
  if (cohort == "UGLI2") {
    geno_array <- "FinnGen Thermo Fisher Axiom custom array"
    name <- "life3"
    }
  
  outfile <-
    paste0(name, "_", proxy_casecontrol, "_", parent, "_", analysis, "_eur")
  
  ## only run if outfile does not exist
  if (!file.exists(here(
    "data_raw",
    "sumstats",
    proxy_casecontrol,
    paste0(outfile, ".txt.gz")
  ))) {
    
    n_case <- Ns$N[Ns$cohort == cohort & Ns$parent == parent & Ns$case_control == "case"]
    n_control <- Ns$N[Ns$cohort == cohort & Ns$parent == parent & Ns$case_control == "control"]
    
    ## load & format sumstats
    print(paste0("load & format sumstats: ", infile))
    fread(infile) %>% 
      ## check that sample size is correct
      {if (max(.$n, na.rm = T) != (n_case + n_control)) stop("n_total != (n_case + n_control)") else print("n_total == (n_case + n_control)");.} %>%
      mutate(
        n_case = n_case,
        n_control = n_control
      ) %>%
      mutate(
        neff = neff_f(n_case = n_case,
                      n_control = n_control,
                      proxy_casecontrol = proxy_casecontrol),
        variant_id = paste(
          chromosome, 
          base_pair_location, 
          other_allele, 
          effect_allele, 
          sep = ":")
      ) %>%
      rename(rsid = rs_id) %>%
      select(
        variant_id,
        chromosome,
        base_pair_location,
        other_allele,
        effect_allele,
        beta,
        standard_error,
        effect_allele_frequency,
        p_value,
        info,
        neff,
        n_case,
        n_control
      ) %>%
      arrange(chromosome, base_pair_location) %>%
      {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
      fwrite(
        file = here(
          "data_raw",
          "sumstats",
          proxy_casecontrol,
          paste0(outfile, ".txt.gz")
        ),
        quote = F,
        sep = "\t",
        na = "NA",
        row.names = F
      )
    
    print(paste0(
      "saved in: ",
      here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, ".txt.gz")
      )
    ))
    
    ## save readme file
    readme <-
      file(here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, "_readme.txt")
      ))
    writeLines(c(
      paste0("Date: ", Sys.Date()),
      "- genomeAssembly: GRCh37/hg19",
      "- traitDescription: 
            Proxy analysis: We defined the Maternal dementia phenotype by assigning an any individual 
            reporting that their mother suffered from dementia as a case (report a 
            1 in 'dementia_mother_fam_q_1_a' or 'dementia_mother_fam_q_1_b'). 
            We repeated this for the Paternal phenotype using 'dementia_father_fam_q_1_a' 
            and 'dementia_father_fam_q_1_b. To avoid overlapping individuals between 
            the phenotypes, all individuals who were cases in both the Maternal 
            and Paternal phenotypes were removed from the Maternal phenotype, 
            individuals who were a case in one phenotype and a control in another 
            were removed from the phenotype where they were a control, and 
            individuals who were controls in both phenotypes were split in half 
            and randomly assigned to one of the phenotypes.",
      paste0("- sampleSize: ", n_case + n_control),
      paste0("- caseCount: ", n_case),
      paste0("- controlCount: ", n_control),
      "- caseControlStudy: Yes",
      "- sampleAncestry: European (Dutch)",
      paste0("- genotypingTechnology: ", geno_array),
      "- analysisSoftware: regenie v3.4",
      "- imputationPanel: Genome of The Netherlands (GoNL) release 5 and the 1000 Genomes phase1 v3",
      "- imputationSoftware: Beagle 3.1.0 using Minimac version 2012.10.3",
      paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
      "- ancestryMethod: 
            To estimate ancestry, the datasets were merged with the 1000 Genomes 
            (1KG) reference data and individuals outside of 4 standard deviations 
            of the mean 1KG European value of PC1 and PC2 were removed 
            (first 5 PCs in the UGLI2 dataset). ",
      "- sortedByGenomicLocation: Yes",
      "- dosageOrHardcalls: Dosage (except for X chromosome)",
      paste0("- covariates: ", covars_readme),
      "- notes: 
            Because X chromomosomes are coded differently between males and
            females, it was analysed as hard calls. Autosomes were analysed
            as dosages."
    ), 
    readme)
    close(readme)
  } else {
    
    ## print if code is skipped because file already exists
    print(paste0(outfile, ": already exists")) 
    
  }
}

## sample sizes based on QCtrackindividuals.xlsx
Ns <- expand.grid(
  cohort = c("GWAS", "UGLI1", "UGLI2"),
  parent = c("mother", "father"),
  case_control = c("case", "control"),
  stringsAsFactors = F
) %>%
  mutate(N = case_when(
    cohort == "GWAS" & parent == "mother" & case_control == "case" ~ 444,
    cohort == "UGLI1" & parent == "mother" & case_control == "case" ~ 907,
    cohort == "UGLI2" & parent == "mother" & case_control == "case" ~ 896,
    cohort == "GWAS" & parent == "father" & case_control == "case" ~ 284,
    cohort == "UGLI1" & parent == "father" & case_control == "case" ~ 556,
    cohort == "UGLI2" & parent == "father" & case_control == "case" ~ 625,
    cohort == "GWAS" & parent == "mother" & case_control == "control" ~ 2877,
    cohort == "UGLI1" & parent == "mother" & case_control == "control" ~ 7187,
    cohort == "UGLI2" & parent == "mother" & case_control == "control" ~ 7812,
    cohort == "GWAS" & parent == "father" & case_control == "control" ~ 2803,
    cohort == "UGLI1" & parent == "father" & case_control == "control" ~ 7079,
    cohort == "UGLI2" & parent == "father" & case_control == "control" ~ 7720,
  ))

## parameters defining different sumstats
params <- expand.grid(
  analysis = c(
    "SexPCs",
    "SexPCsMother_ageMother_age2AgeAge2APOE4",
    "SexPCsMother_ageMother_age2AgeAge2"
  ),
  cohort = c("GWAS", "UGLI1", "UGLI2"),
  parent = c("Maternal"),
  stringsAsFactors = F
) %>% 
  rbind(
    expand.grid(
      analysis = c(
        "SexPCs",
        "SexPCsFather_ageFather_age2AgeAge2APOE4",
        "SexPCsFather_ageFather_age2AgeAge2"
      ),
      cohort = c("GWAS", "UGLI1", "UGLI2"),
      parent = c("Paternal"),
      stringsAsFactors = F
    )
  )

params %>%
  pmap(~ lifelines_sumstats(
    analysis = ..1,
    cohort = ..2,
    parent = ..3,
    proxy_casecontrol = "proxy",
    Ns = Ns
  ))

## clear workspace
rm(list=setdiff(ls(), c("neff_f", "snp_modifyBuild_local")))


#### 20. AllOfUs ####

## function to load, format, process, save, and plot jpb-sumstats
allofus_sumstats <- function(proxy_casecontrol, analysis, ancestry) {
  
  ## set file names
  if (proxy_casecontrol == "proxy") {
    infile <-
      Sys.glob(here(
        "data_raw",
        "cohorts",
        "allofus",
        proxy_casecontrol,
        paste0("proxy_", analysis, "_", ancestry, "*.glm.logistic.gz")
      ))
  } else {
    infile <-
      Sys.glob(here(
        "data_raw",
        "cohorts",
        "allofus",
        proxy_casecontrol,
        paste0(analysis, "_", ancestry, "*.glm.logistic.gz")
      ))
  }

  
  covars_readme <- system(paste0(
    "grep 'Covariates Used:' ",
    here(
      "data_raw",
      "cohorts",
      "allofus",
      proxy_casecontrol,
      paste0("*", analysis, "_", ancestry, "*.glm.logistic.meta")
    )
  ), intern = T) %>% str_remove("Covariates Used:\t")
  
  ## name of file to be saved
  # case_control
  if (analysis == "bothSex") {analysis <- "main"}
  if (analysis == "female") {analysis <- "fema"}
  if (analysis == "bothSex_withoutAge") {analysis <- "noag"}
  # proxy
  if (analysis == "father") {analysis <- "father_main"}
  if (analysis == "father_withoutAge") {analysis <- "father_noag"}
  if (analysis == "mother") {analysis <- "mother_main"}
  if (analysis == "mother_withoutAge") {analysis <- "mother_noag"}
  outfile <- paste("allof", proxy_casecontrol, analysis, ancestry, sep = "_")
  
  ## only run if outfile does not exist
  if (!file.exists(here(
    "data_raw",
    "sumstats",
    proxy_casecontrol,
    paste0(outfile, ".txt.gz")
  ))) {
    
    ## load & format sumstats
    print(paste0("load & format sumstats: ", infile))
    fread(infile) %>% 
      mutate(
        variant_id_hg38 = paste(
          chromosome, 
          base_pair_location, 
          other_allele, 
          effect_allele, 
          sep = ":")
      ) %>%
      ## change column names for liftover
      rename(chr = chromosome,
             pos = base_pair_location) %>%
      snp_modifyBuild_local(liftOver = here("src", "liftOver"),
                      from = "hg38",
                      to = "hg19",
                      check_reverse = T) %>%
      rename(chromosome = chr,
             base_pair_location = pos) %>%
      {max(.$n_control, na.rm = T) ->> n_control;.} %>%
      {max(.$n_case, na.rm = T) ->> n_case;.} %>%
      mutate(
        neff = neff_f(n_case = n_case,
                      n_control = n_control,
                      proxy_casecontrol = proxy_casecontrol),
        variant_id = paste(
          chromosome, 
          base_pair_location, 
          other_allele, 
          effect_allele, 
          sep = ":")
      ) %>%
      select(
        variant_id,
        chromosome,
        base_pair_location,
        other_allele,
        effect_allele,
        beta,
        standard_error,
        effect_allele_frequency,
        p_value,
        info,
        neff,
        n_case,
        n_control,
        variant_id_hg38
      ) %>%
      arrange(chromosome, base_pair_location) %>%
      {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
      fwrite(
        file = here(
          "data_raw",
          "sumstats",
          proxy_casecontrol,
          paste0(outfile, ".txt.gz")
        ),
        quote = F,
        sep = "\t",
        na = "NA",
        row.names = F
      )
    
    print(paste0(
      "saved in: ",
      here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, ".txt.gz")
      )
    ))
    
    if (proxy_casecontrol == "case_control") {
      caseControlStudy <- "yes"
    } else {
      caseControlStudy <- "no"
    }
    
    ## save readme file
    readme <-
      file(here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, "_readme.txt")
      ))
    writeLines(c(
      paste0("Date: ", Sys.Date()),
      "- genomeAssembly: GRCh37/hg19",
      "- traitDescription: 
          True case definition:
          Cases were mainly defined by the self-reported questionnaire and the 
          predetermined Alzheimer's related code from All of Us. We used the 
          following 'concept codes' to capture AD cases:
            378419, 43530664, 4220313 and 4278830
          each of these codes can be linked to a corresponding SNOMED code (a US 
          Gov standard) and OMOP concept, details can be found in the following site:
          https://databrowser.researchallofus.org/ehr/conditions
      
          As for self-reported questionnaire, we extracted results from the following 
          questions (concept ID included):
            - Has a doctor or health care provider ever told you that you have? (select 
              all that apply) (1384526) [ One box for: Dementia (including AD, 
              Vascular Dementia..) ] 
            - About how old were you when you were first told you had Dementia (43530525)
            - Are you currently prescribed medications and/or receiving treatment 
              for Dementia (43528852)
            - Are you still seeing a doctor or health care provider for Dementia (43530367)
      
          We further perform the following filtering to refine to AD:
            - Remove all reports related to 'Parkinson', 'Lewy body' or 'Mixed dementia'
            - For self-report, the answers do not contain the word Alzheimer's disease, 
              but include traits such as Neuropathy, Traumatic Brain Injury etc. 
              As such we extract any answers related to 'Dementia' and/or 'Memory Loss'.
      
          Proxy-defintion:
          Those were all defined by self-reported questionnaire about family history 
          (https://www.researchallofus.org/data-tools/survey-explorer/personal-family-health-history-survey/ & 
          https://databrowser.researchallofus.org/survey/personal-and-family-health-history/dementia). 
          They were asked 'Have you or anyone in your family ever been diagnosed with 
          the following brain and nervous system conditions? Think only of the people 
          you are related to by blood. - Dementia (includes Alzheimer's, vascular, etc.)'
          ",
      paste0("- sampleSize: ", n_case + n_control),
      paste0("- caseCount: ", n_case),
      paste0("- controlCount: ", n_control),
      paste0("- caseControlStudy: ", caseControlStudy),
      paste0("- sampleAncestry: ", ancestry),
      "- genotypingTechnology: Illumina Global Diversity Array (GDA)",
      "- analysisSoftware: plink2",
      "- imputationPanel: -",
      "- imputationSoftware: -",
      paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
      "- ancestryMethod: PCA",
      "- dosageOrHardcalls: Hard Calls",
      paste0("- covariates: ", covars_readme),
      "- notes: 
          - no imputation because of high densitiy genotyping array
          - ~16 variants cannot be mapped"
    ), 
    readme)
    close(readme)
  } else {
    
    ## print if code is skipped because file already exists
    print(paste0(outfile, ": already exists")) 
    
  }
}

## parameters defining different sumstats
case_control_params <- expand.grid(
  proxy_casecontrol = "case_control",
  analysis = c("bothSex", "bothSex_withoutAge", "female", "male"),
  ancestry = c("eur", "afr", "amr")
)

proxy_params <- expand.grid(
  proxy_casecontrol = "proxy",
  analysis = c("father", "father_withoutAge", "mother", "mother_withoutAge"),
  ancestry = c("eur", "afr", "amr")
)

params <- rbind(
  case_control_params,
  proxy_params
)

params %>%
  pmap(~ allofus_sumstats(
    proxy_casecontrol = ..1,
    analysis = ..2,
    ancestry = ..3
  ))

## clear workspace
rm(list=setdiff(ls(), c("neff_f", "snp_modifyBuild_local")))


#### 21. Eli Lilly, Johnson & Johnson, PROTECT UK #### 

## function to load, format, process, save, and plot jpb-sumstats
lilly_sumstats <- function(proxy_casecontrol, analysis, analysis_nr, ancestry) {
  
  ## set file names
    infile <-
      Sys.glob(here(
        "data_raw",
        "cohorts",
        "lilly_jj_protect",
        "transfer_93587_files_cf928460",
        paste0("pharma_protect.pgc_setup_", analysis_nr,".ALZ.glm.logistic.hybrid.gz")
      ))
  
  ## name of file to be saved
  outfile <- paste("lilly", proxy_casecontrol, analysis, ancestry, sep = "_")
  
  ## only run if outfile does not exist
  if (!file.exists(here(
    "data_raw",
    "sumstats",
    proxy_casecontrol,
    paste0(outfile, ".txt.gz")
  ))) {
    
    ## load & format sumstats
    print(paste0("load & format sumstats: ", infile))
    fread(infile) %>% 
      mutate(
        variant_id_hg38 = paste(
          `#CHROM`, 
          POS, 
          OMITTED, 
          A1, 
          sep = ":")
      ) %>%
      ## change column names (for liftover)
      rename(
        chr = "#CHROM",
        pos = "POS",
        rsid = RSID,
        effect_allele = A1,
        other_allele = OMITTED,
        effect_allele_frequency = A1_FREQ,
        p_value = P,
        n_case = NCASE,
        n_control = NCONTROL,
        standard_error = "LOG(OR)_SE") %>%
      snp_modifyBuild_local(liftOver = here("src", "liftOver"),
                            from = "hg38",
                            to = "hg19",
                            check_reverse = T) %>%
      rename(chromosome = chr, base_pair_location = pos) %>%
      {max(.$n_control, na.rm = T) ->> n_control;.} %>%
      {max(.$n_case, na.rm = T) ->> n_case;.} %>%
      mutate(
        neff = neff_f(n_case = n_case,
                      n_control = n_control,
                      proxy_casecontrol = proxy_casecontrol),
        variant_id = paste(
          chromosome, 
          base_pair_location, 
          other_allele, 
          effect_allele, 
          sep = ":"),
        beta = log(OR)
      ) %>%
      select(
        variant_id,
        chromosome,
        base_pair_location,
        other_allele,
        effect_allele,
        beta,
        standard_error,
        effect_allele_frequency,
        p_value,
        neff,
        n_case,
        n_control,
        variant_id_hg38, 
        rsid
      ) %>%
      arrange(chromosome, base_pair_location) %>%
      {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
      fwrite(
        file = here(
          "data_raw",
          "sumstats",
          proxy_casecontrol,
          paste0(outfile, ".txt.gz")
        ),
        quote = F,
        sep = "\t",
        na = "NA",
        row.names = F
      )
    
    print(paste0(
      "saved in: ",
      here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, ".txt.gz")
      )
    ))
    
    if (proxy_casecontrol == "case_control") {
      caseControlStudy <- "yes"
    } else {
      caseControlStudy <- "no"
    }
    
    ## save readme file
    readme <-
      file(here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, "_readme.txt")
      ))
    writeLines(c(
      paste0("Date: ", Sys.Date()),
      "- genomeAssembly: GRCh37/hg19",
      "- traitDescription: ?",
      paste0("- sampleSize: ", n_case + n_control),
      paste0("- caseCount: ", n_case),
      paste0("- controlCount: ", n_control),
      paste0("- caseControlStudy: ", caseControlStudy),
      paste0("- sampleAncestry: ", ancestry),
      "- genotypingTechnology: Different for different genotyping batches as listed below.
	      Eli Lilly Expedition3 trial (ALZ cases): GSA+Multi disease panel version 2 (one batch)
	      Eli Lilly old trials (ALZ cases): Illumina 5M panel version 3 (three batches)
	      J&J (ALZ cases): GDAneuroBooster array. Custom modification of Illumina ",
      "- analysisSoftware: plink2 (v2.00a6LM AVX2 AMD (9 Jun 2024))",
      "- imputationPanel: HGDP + 1KG (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9900804/)",
      "- imputationSoftware: -",
      paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
      "- ancestryMethod: PCA",
      "- dosageOrHardcalls: Dosages",
      paste0("- covariates: -"),
      "- notes: 
          - "
    ), 
    readme)
    close(readme)
  } else {
    
    ## print if code is skipped because file already exists
    print(paste0(outfile, ": already exists")) 
    
  }
}

## parameters defining different sumstats
case_control_params <- data.table(
  proxy_casecontrol = rep("case_control", 5),
  analysis = c("main", "male", "fema", "noag", "apoe"),
  analysis_nr = c(1, 2 ,3 ,4 , 5),
  ancestry = rep("eur", 5)
)

case_control_params %>%
  pmap(~ lilly_sumstats(
    proxy_casecontrol = ..1,
    analysis = ..2,
    analysis_nr = ..3,
    ancestry = ..4
  ))


#### 22. demGene ####

## function to load, format, process, save, and plot sumstats
demgene_sumstats <- function(proxy_casecontrol, analysis, ancestry) {
  
  ## set file names
  infile <-
    Sys.glob(here(
      paste0(
        "data_raw/cohorts/demgene/sumstats/demgene_OmniExpress_noHUNT_",
        analysis,
        "MAF.txt.gz"
      )
    ))
  
  if (analysis == "sex") {analysis <- "noag"}
  if (analysis == "APOE4") {analysis <- "apoe"}
  if (analysis == "males") {analysis <- "male"}
  if (analysis == "females") {analysis <- "fema"}
  
  ## name of file to be saved
  outfile <- paste("demge", proxy_casecontrol, analysis, ancestry, sep = "_")
  
  ## only run if outfile does not exist
  if (!file.exists(here(
    "data_raw",
    "sumstats",
    proxy_casecontrol,
    paste0(outfile, ".txt.gz")
  ))) {
    
    ## load & format sumstats
    print(paste0("load & format sumstats: ", infile))
    fread(infile) %>% 
      mutate(
        variant_id_hg38 = paste(
          chromosome, 
          base_pair_location, 
          other_allele, 
          effect_allele, 
          sep = ":")
      ) %>%
      ## change column names for liftover
      rename(
        chr = "chromosome",
        pos = "base_pair_location"
        ) %>%
      snp_modifyBuild_local(
        liftOver = here("src", "liftOver"),
        from = "hg38",
        to = "hg19",
        check_reverse = T
      ) %>%
      rename(chromosome = chr, base_pair_location = pos) %>% 
      {max(.$n_control, na.rm = T) ->> n_control;.} %>%
      {max(.$n_case, na.rm = T) ->> n_case;.} %>%
      mutate(
        neff = neff_f(
          n_case = n_case,
          n_control = n_control,
          proxy_casecontrol = !!proxy_casecontrol
        ),
        variant_id = paste(
          chromosome,
          base_pair_location,
          other_allele,
          effect_allele,
          sep = ":"
        )
      ) %>% 
      select(
        variant_id,
        chromosome,
        base_pair_location,
        other_allele,
        effect_allele,
        beta,
        standard_error,
        effect_allele_frequency,
        p_value,
        info,
        neff,
        n_case,
        n_control,
        variant_id_hg38     
        ) %>%
      arrange(chromosome, base_pair_location) %>%
      {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
      fwrite(
        file = here(
          "data_raw",
          "sumstats",
          proxy_casecontrol,
          paste0(outfile, ".txt.gz")
        ),
        quote = F,
        sep = "\t",
        na = "NA",
        row.names = F
      )
    
    print(paste0(
      "saved in: ",
      here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, ".txt.gz")
      )
    ))
    
    if (proxy_casecontrol == "case_control") {
      caseControlStudy <- "yes"
    } else {
      caseControlStudy <- "no"
    }
    
    if (analysis == "main") {covars <- "PCs1-20,age,age2Msex,age*sex"}
    if (analysis == "noag") {covars <- "PCs1-20,sex"}
    if (analysis == "apoe") {covars <- "PCs1-20,age,age2,sex,age*sex,APOE4"}
    if (analysis == "male") {covars <- "PCs1-20,age,age2"}
    if (analysis == "fema") {covars <- "PCs1-20,age,age2"}
    
    ## save readme file
    readme <-
      file(here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, "_readme.txt")
      ))
    writeLines(c(
      paste0("Date: ", Sys.Date()),
      "- genomeAssembly: GRCh37/hg19",
      "- traitDescription: -",
      paste0("- sampleSize: ", n_case + n_control),
      paste0("- caseCount: ", n_case),
      paste0("- controlCount: ", n_control),
      paste0("- caseControlStudy: ", caseControlStudy),
      paste0("- sampleAncestry: ", ancestry),
      "- genotypingTechnology: 
      HumanOmniExpress, HumanOmniExpress-12v1-1B/HumanOmniExpress-12v1_H, 
      HumanOmniExpress-12v1-1B/HumanOmniExpress-24v1-0A, HumanOmniExpress-24v1-0A, 
      HumanOmniExpress-24v1-1_A, HumanOmniExpress-24v1-2_A1, InfiniumOmniExpress-24v1-2_A1",
      "- imputationPanel: TOPMed",
      "- imputationSoftware: -",
      paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
      "- ancestryMethod: -",
      "- dosageOrHardcalls: dosages",
      paste0("- covariates: ", covars), 
      "- notes: 
          - demgene_1706_batch11 & demgene_DeCodeGenetics_GSA_no1706_batch11 
            had QC problems, so we decided to exclude them and only analyse
            demgene_OmniExpress"
    ), 
    readme)
    close(readme)
  } else {
    
    ## print if code is skipped because file already exists
    print(paste0(outfile, ": already exists")) 
    
  }
}

## parameters defining different sumstats
case_control_params <- expand.grid(
  proxy_casecontrol = "case_control",
  analysis = c("main", "males", "females", "sex", "APOE4"),
  ancestry = "eur"
)

case_control_params %>%
  pmap(~ demgene_sumstats(
    proxy_casecontrol = ..1,
    analysis = ..2,
    ancestry = ..3
  ))


#### 23. HUSK ####

## set file names
proxy_casecontrol <- "case_control"

infile <- here(
  "data_raw",
  "cohorts",
  "husk/NOR_HUSK_Alzheimers_ICD10_F00_G30_ICD9_331_june2024",
  "Norway_HUSK_Alzheimers_ICD10_F00_G30_ICD9_331_vs_Controls_age70plus_AdjYOB_summarystats_27062024.txt.gz"
)

outfile <- paste0("xhusk_", proxy_casecontrol, "_main_eur")

# sample size reported in:
# "Norway_HUSK_Alzheimers_ICD10_F00_G30_ICD9_331_vs_Controls_age70plus_summarystats_info_27062024.txt"
n_case <- 624
n_control <- 15329

## only run if outfile does not exist
if (!file.exists(here(
  "data_raw",
  "sumstats",
  proxy_casecontrol,
  paste0(outfile, ".txt.gz")
))) {
  
  ## load & format sumstats
  print(paste0("load & format sumstats: ", infile))
  fread(infile) %>% 
    mutate(
      Chr = ifelse(Chr == "chrX", "chr23", Chr)
    ) %>%
    mutate(
      Chr = as.integer(str_remove(Chr, "chr")),
      n_case = !!n_case,
      n_control = !!n_control,
      EAFrq = EAFrq / 100
    ) %>%
    mutate(
      variant_id_hg38 = paste(Chr, Pos, OA, EA, sep = ":")
    ) %>%
    ## change column names for liftover
    rename(chr = Chr,
           pos = Pos) %>%
    snp_modifyBuild_local(
      liftOver = here("src", "liftOver"),
      from = "hg38",
      to = "hg19",
      check_reverse = T
    ) %>% 
    rename(
      chromosome = chr,
      base_pair_location = pos,
      other_allele = OA,
      effect_allele = EA, # minor allele
      beta = Beta,
      standard_error = SE,
      effect_allele_frequency = EAFrq,
      p_value = P,
      rsid = rsID
    ) %>% 
    mutate(
      variant_id = paste(
        chromosome, 
        base_pair_location, 
        other_allele, 
        effect_allele, 
        sep = ":"),
      neff = neff_f(
        n_case = n_case,
        n_control = n_control,
        proxy_casecontrol = !!proxy_casecontrol
        )
    ) %>%
    select(
      variant_id,
      chromosome,
      base_pair_location,
      other_allele,
      effect_allele,
      beta,
      standard_error,
      effect_allele_frequency,
      p_value,
      info,
      neff,
      n_case,
      n_control,
      variant_id_hg38,
      rsid
    ) %>%
    arrange(chromosome, base_pair_location) %>%
    {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
    fwrite(
      file = here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, ".txt.gz")
      ),
      quote = F,
      sep = "\t",
      na = "NA",
      row.names = F
    )
  
  ## save readme file
  readme <-
    file(here(
      "data_raw",
      "sumstats",
      proxy_casecontrol,
      paste0(outfile, "_readme.txt")
    ))
  writeLines(c(
    paste0("Date: ", Sys.Date()),
    "- genomeAssembly: GRCh37/hg19",
    "- traitDescription: Alzheimers-ICD10-F00,G30 and/or ICD9-331 vs Controls that have reached age 70 or over",
    paste0("- sampleSize: ", n_case + n_control),
    paste0("- caseCount: ", n_case),
    paste0("- controlCount: ", n_control),
    "- caseControlStudy: Yes",
    "- sampleAncestry: European",
    "- genotypingTechnology: ?",
    "- analysisSoftware: ?",
    "- imputationPanel: ?",
    "- imputationSoftware: ?",
    paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
    "- ancestryMethod: ?",
    "- sortedByGenomicLocation: Yes",
    "- dosageOrHardcalls: ?",
    "- covariates: ?",
    "- notes: 
      -  Male and female combined
	                male	female	mean year-of-birth
        Case	    214	  410	    1931.6
        Control	  7240	8089	  1944.2
    
      - 401,714 variants have not been mapped with liftOver"
  ), 
  readme)
  close(readme)
} else {
  
  ## print if code is skipped because file already exists
  print(paste0(outfile, ": already exists")) 
  
}

## clear workspace
rm(list=setdiff(ls(), c("neff_f", "snp_modifyBuild_local")))



#### 24. Copenhagen Hospital Biobank (xchbb) ####

## function to load, format, process, save, and plot sumstats
xchbb_sumstats <- function(analysis, xchbb_dir) {
  
  ## set file names
  proxy_casecontrol <- "case_control"
  
  infile <- here(
    paste0(
      "data_raw/cohorts/copenhagen_hospital_biobank/",
      xchbb_dir,
      "/",
      "PGC_GWAS_CHB_",
      analysis,
      "_",
      gsub(
        pattern = "_",
        replacement = "",
        x = xchbb_dir
      ),
      # remove underscore
      ".gz"
    )
  )
  
  if (analysis == "model1_main") {analysis <- "main"}
  if (analysis == "model2_males") {analysis <- "male"}
  if (analysis == "model3_females") {analysis <- "fema"}
  if (analysis == "model4_woage") {analysis <- "noag"}
  if (analysis == "model5_apoe") {analysis <- "apoe"}
  
  
  outfile <- paste0("xchbb_case_control_", analysis, "_eur")
  
  ## only run if outfile does not exist
  if (!file.exists(here(
    "data_raw",
    "sumstats",
    proxy_casecontrol,
    paste0(outfile, ".txt.gz")
  ))) {
    
    ## load & format sumstats
    print(paste0("load & format sumstats: ", infile))
    fread(infile) %>% 
      {max(.$n_control, na.rm = T) ->> n_control;.} %>%
      {max(.$n_case, na.rm = T) ->> n_case;.} %>%
      mutate(
        variant_id_hg38 = paste(
          chromosome, 
          base_pair_location, 
          other_allele, 
          effect_allele, 
          sep = ":")
      ) %>%
      ## change column names for liftover
      rename(chr = chromosome,
             pos = base_pair_location) %>%
      snp_modifyBuild_local(
        liftOver = here("src", "liftOver"),
        from = "hg38",
        to = "hg19",
        check_reverse = T
      ) %>% 
      rename(
        chromosome = chr,
        base_pair_location = pos
      ) %>% 
      mutate(
        variant_id = paste(
          chromosome, 
          base_pair_location, 
          other_allele, 
          effect_allele, 
          sep = ":"),
        neff = neff_f(
          n_case = n_case,
          n_control = n_control,
          proxy_casecontrol = !!proxy_casecontrol
        )
      ) %>%
      ## there are some SNPs with neff ~ 0. I will remove these SNPs manually
      filter(neff > 100) %>%
      select(
        variant_id,
        chromosome,
        base_pair_location,
        other_allele,
        effect_allele,
        beta,
        standard_error,
        effect_allele_frequency,
        p_value,
        info,
        neff,
        n_case,
        n_control,
        variant_id_hg38
      ) %>%
      arrange(chromosome, base_pair_location) %>%
      {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
      fwrite(
        file = here(
          "data_raw",
          "sumstats",
          proxy_casecontrol,
          paste0(outfile, ".txt.gz")
        ),
        quote = F,
        sep = "\t",
        na = "NA",
        row.names = F
      )
    
    ## save readme file
    readme <-
      file(here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, "_readme.txt")
      ))
    writeLines(c(
      paste0("Date: ", Sys.Date()),
      "- genomeAssembly: GRCh37/hg19",
      "- traitDescription: closely followed 2024_05_13_pgcalz3_analysis_plan_v2.6",
      paste0("- sampleSize: ", n_case + n_control),
      paste0("- caseCount: ", n_case),
      paste0("- controlCount: ", n_control),
      "- caseControlStudy: Yes",
      "- sampleAncestry: European (Danish)",
      "- genotypingTechnology: Illumina GSA chip",
      "- analysisSoftware: Regenie/3.3",
      "- imputationPanel:",
      "- imputationSoftware:",
      paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
      "- ancestryMethod:",
      "- sortedByGenomicLocation:",
      "- dosageOrHardcalls: Dosage",
      "- covariates: chip batch, PC1-20s + per-analysis covariates",
      "- notes
      ~327,726 SNPs could not be mapped: 
"
    ), 
    readme)
    close(readme)
  } else {
    
    ## print if code is skipped because file already exists
    print(paste0(outfile, ": already exists")) 
    
  }
}

## parameters defining different sumstats
xchbb_params <- data.frame(
  analysis = c(
    "model1_main",
    "model2_males",
    "model3_females",
    "model4_woage",
    "model5_apoe"
  ),
  xchbb_dir = c(
    "2024_09_24",
    "2024_10_24",
    "2024_10_24",
    "2024_10_24",
    "2024_10_24"
  )
)

xchbb_params %>%
  pmap(~ xchbb_sumstats(
    analysis = ..1,
    xchbb_dir = ..2
  ))


## clear workspace
rm(list=setdiff(ls(), c("neff_f", "snp_modifyBuild_local")))



#### 25. EADB bellenguez 2022 (xeadb) ####

## set file names
proxy_casecontrol <- "combined"

infile <- here(
  "data_raw/published_sumstats/ad_bellenguez",
  "GCST90027158_buildGRCh38.tsv.gz"
)

outfile <- paste0("xeadb_", proxy_casecontrol, "_noag_eur")

## only run if outfile does not exist
if (!file.exists(here(paste0(
  "data_raw/sumstats/combined/", outfile, ".txt.gz"
)))) {
  
  ## load & format sumstats
  print(paste0("load & format sumstats: ", infile))
  fread(infile) %>% 
    filter(effect_allele_frequency > 0.01) %>%
    mutate(
      variant_id_hg38 = paste(
        chromosome, 
        base_pair_location, 
        other_allele, 
        effect_allele, 
        sep = ":")
    ) %>%
    ## change column names for liftover
    rename(
      chr = chromosome,
      pos = base_pair_location,
      n_case = n_cases,
      n_control = n_controls
    ) %>%
    snp_modifyBuild_local(
      liftOver = here("src", "liftOver"),
      from = "hg38",
      to = "hg19",
      check_reverse = T
    ) %>% 
    rename(
      chromosome = chr,
      base_pair_location = pos
    ) %>% 
    
    # # # based on supplementary table 2 in Bellenguez 2022 # # #
    mutate(
      n_proxy_case = case_when(
        n_case > 40000 ~ 49275,
        .default = 0
      ),
      n_case = case_when(
        n_case > 40000 ~ n_case - 49275,
        .default = n_case
      ),
      n_proxy_control = case_when(
        n_control > 100000 ~ 338440,
        .default = 0
      ),
      n_control = case_when(
        n_control > 100000 ~ n_control - 338440,
        .default = n_control
      ),
    # # # based on supplementary table 2 in Bellenguez 2022 # # #
    
      variant_id = paste(
        chromosome, 
        base_pair_location, 
        other_allele, 
        effect_allele, 
        sep = ":"),
      neff = round( (4 / ((1 / n_case) + (1 / n_control))) + ((4 / ((1 / n_proxy_case) + (1 / n_proxy_control))) / 4) )
    ) %>%
    {max(.$n_control, na.rm = T) ->> n_control;.} %>%
    {max(.$n_case, na.rm = T) ->> n_case;.} %>%
    {max(.$n_proxy_control, na.rm = T) ->> n_proxy_control;.} %>%
    {max(.$n_proxy_case, na.rm = T) ->> n_proxy_case;.} %>%
    ## there are some SNPs with neff ~ 0. I will remove these SNPs manually
    filter(neff > 100) %>%
    select(
      variant_id,
      chromosome,
      base_pair_location,
      other_allele,
      effect_allele,
      beta,
      standard_error,
      effect_allele_frequency,
      p_value,
      neff,
      n_case,
      n_control,
      n_proxy_case,
      n_proxy_control,
      variant_id_hg38
    ) %>%
    arrange(chromosome, base_pair_location) %>%
    {min(.$effect_allele_frequency, na.rm = T) ->> min_effect_allele_frequency;.} %>%
    fwrite(
      file = here(
        "data_raw",
        "sumstats",
        proxy_casecontrol,
        paste0(outfile, ".txt.gz")
      ),
      quote = F,
      sep = "\t",
      na = "NA",
      row.names = F
    )
  
  ## save readme file
  readme <-
    file(here(
      "data_raw",
      "sumstats",
      proxy_casecontrol,
      paste0(outfile, "_readme.txt")
    ))
  writeLines(c(
    paste0("Date: ", Sys.Date()),
    "- genomeAssembly: GRCh37/hg19",
    "- traitDescription: closely followed 2024_05_13_pgcalz3_analysis_plan_v2.6",
    paste0("- sampleSize: ", n_case + n_control + n_proxy_control + n_proxy_case),
    paste0("- caseCount: ", n_case),
    paste0("- controlCount: ", n_control),
    paste0("- proxycaseCount: ", n_proxy_case),
    paste0("- proxycontrolCount: ", n_proxy_control),
    "- caseControlStudy: case-control + proxy",
    "- sampleAncestry: European",
    "- genotypingTechnology: mixed",
    "- analysisSoftware: plink1.9, SNPTEST, SAIGE ",
    "- imputationPanel:",
    "- imputationSoftware:",
    paste0("- effectAlleleFreqLowerLimit: ", min_effect_allele_frequency),
    "- ancestryMethod:",
    "- sortedByGenomicLocation:",
    "- dosageOrHardcalls: Dosage",
    "- covariates: see Supplementary Table 2. variable number of PCs, no age,
        no sex, sometimes batch, sometimes no covariats at all.",
    "- notes:
      proxies are simply coded as 0 and 1, betas and SEs are multiplied by 2
      ~32,942 SNPs could not be mapped: "
  ), 
  readme)
  close(readme)
} else {
  
  ## print if code is skipped because file already exists
  print(paste0(outfile, ": already exists")) 
  
}

## clear workspace
rm(list=setdiff(ls(), c("neff_f", "snp_modifyBuild_local")))