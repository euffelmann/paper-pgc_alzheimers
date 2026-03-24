# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                         #
#         script to process raw genotype data             #
#         run remotely                                    #
#                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# - - - - - Table of Contents - - - - - - - - - - - - - - #
# - - - - - 00. set-up - - - - - - - - - - - - - - - - -  #
# - - - - - 01. DemGene - - - - - - - - - - - - - - - - - #
  # -------- a. make covariate file --------------------- #
  # -------- b. clean DemGene batches ------------------- #
  # -------- c. combine individual batches -------------- #
  # -------- d. PCA plots ------------------------------- #
  # -------- e. save data ------------------------------- #

#### 00. set-up ####

library(here)
library(dplyr)
library(data.table)
library(stringr)
plink1.9 <- here("src", "plink1.9", "plink1.9")
plink2 <- here("src", "plink2", "plink2")

#### 01. DemGene ####

##### a. make covariate file #####

if (!file.exists(here("data_raw", "cohorts", "demgene", "geno", "demgene.cov"))) {
  
  ## load all covariates + AD status
  cases <- fread(
    here(
      "data_raw",
      "cohorts",
      "demgene",
      "demgene4pgc_export_04032024",
      "demgene_redcap_ad.txt"
    )
  ) %>%
    select(
      FID,
      IID,
      demo_diagnosis_group,
      demo_age,
      demo_gender,
    ) %>%
    rename(
      diagnosis = demo_diagnosis_group,
      age = demo_age,
      sex = demo_gender
    ) %>%
    mutate(
      sample = "demgene"
    )
  
  cases_mci <- fread(
    here(
      "data_raw",
      "cohorts",
      "demgene",
      "demgene4pgc_export_04032024",
      "demgene_redcap_mci.txt"
    )
  ) %>%
    select(
      FID,
      IID,
      demo_diagnosis_group,
      demo_age,
      demo_gender,
    ) %>%
    rename(
      diagnosis = demo_diagnosis_group,
      age = demo_age,
      sex = demo_gender
    ) %>%
    mutate(
      sample = "demgene"
    )
  
  controls_demgene <- fread(
    here(
      "data_raw",
      "cohorts",
      "demgene",
      "demgene4pgc_export_04032024",
      "demgene_redcap_control.txt"
    )
  ) %>%
    select(
      FID,
      IID,
      demo_diagnosis_group,
      demo_age,
      demo_gender,
    ) %>%
    rename(
      diagnosis = demo_diagnosis_group,
      age = demo_age,
      sex = demo_gender
    ) %>%
    mutate(
      sample = "demgene"
    )
  
  controls_top <- fread(
    here(
      "data_raw",
      "cohorts",
      "demgene",
      "demgene4pgc_export_04032024",
      "top_controls_redcap_age_sex_patient_status.txt"
    )
  ) %>%
    select(
      FID,
      IID,
      STATUS,
      AGE,
      SEX
    ) %>%
    rename(
      diagnosis = STATUS,
      age = AGE,
      sex = SEX
    ) %>%
    mutate(
      sex = case_when(
        sex == "male" ~ 1,
        sex == "female" ~ 2,
        .default = 0
      ),
      sample = "top"
    )
  
  ## make covariate file + batch + AD status
  temp <- rbind(cases, cases_mci, controls_demgene, controls_top) %>%
    mutate(batch = sapply(str_split(IID, "__"), `[`, 1))
  
  ## add genotyping platform
  genoplat <- fread(
    here(
      "data_raw",
      "cohorts",
      "demgene",
      "demgene4pgc_export_04032024",
      "batch_genotyping_array_info.txt"
    ),
    header = F
  )
  
  covars <- temp %>%
    inner_join(genoplat, by = c("batch" = "V1")) %>%
    rename(array = V2)
  
  
  fwrite(
    x = covars,
    file = here(
      "data_raw",
      "cohorts",
      "demgene",
      "geno",
      "demgene.cov"
    ),
    quote = F,
    sep = "\t",
    na = "NA",
    row.names = F
  )
}


##### b. clean DemGene batches #####

geno_dir <- here("data_raw",
                 "cohorts",
                 "demgene",
                 "demgene4pgc_export_04032024",
                 "geno/")

## demgene batches
files <- list.files(path = geno_dir,
                    pattern = glob2rx("*.bed"),
                    full.names = F)
files <- gsub(".bed", "", files)

for (batch in files) {
  
  if (!file.exists(paste0(
    here("data_raw", "cohorts", "demgene", "qc/"), batch, "/", batch, "_qc11.bed"
  ))) {
   
    ## make batch-specific qc directory
    system(paste0("mkdir ", here("data_raw", "cohorts", "demgene", "qc/"), batch))
    qc_dir <- here("data_raw", "cohorts", "demgene", "qc", batch)
    
    ###### step 1 (SNP filtering) ###### 
    
    ## include only SNPs with call rate ≥ 0.95 and remove X-chromosome
    system(
      paste0(
        plink1.9, 
        " --bfile ", geno_dir, batch, 
        " --geno 0.05", 
        " --make-bed", 
        " --out ", qc_dir, "/", batch, "_qc1"
      )
    )
    
    ##### change SNP IDs in .bim files  #####  
    # different SNP IDs refer to the same SNP
    # change to CHR:BP format. because multi-allelic SNPs have been removed already
    # I don't need to add the alleles. Moreover, because of strand issues, adding 
    # the alleles causes merge problems (e.g. 2:2947571:A:G becomes 2:2947571:T:C)
    # when strands don't align.
    
    fread(paste0(qc_dir, "/", batch, "_qc1.bim"), header = F) %>%
      mutate(V2 = paste0(V1, ":", V4)) %>%
      fwrite(
        file = paste0(qc_dir, "/", batch, "_qc1.bim"),
        quote = F,
        sep = "\t",
        na = NA,
        col.names = F,
        row.names = F
      )

    
    ## exclude palindromic SNPs (A/T C/G) with MAF > 0.3
    # get allele frequencies
    system(
      paste0(
        plink1.9, 
        " --bfile ", qc_dir, "/", batch, "_qc1", 
        " --freq", 
        " --out ", qc_dir, "/", batch, "_qc1"
      )
    )
    
    # write file with SNPs to remove
    system(
      paste0(
        "awk '{if ($3==\"A\" && $4==\"T\" && $5>0.3) print $2}' ", 
        qc_dir, "/", batch, "_qc1.frq > ", qc_dir, "/", batch, "_qc1_badsnps.txt"
      )
    )
    system(
      paste0(
        "awk '{if ($3==\"T\" && $4==\"A\" && $5>0.3) print $2}' ", 
        qc_dir, "/", batch, "_qc1.frq >> ", qc_dir, "/", batch, "_qc1_badsnps.txt"
      )
    )
    system(
      paste0(
        "awk '{if ($3==\"C\" && $4==\"G\" && $5>0.3) print $2}' ", 
        qc_dir, "/", batch, "_qc1.frq >> ", qc_dir, "/", batch, "_qc1_badsnps.txt"
      )
    )
    system(
      paste0(
        "awk '{if ($3==\"G\" && $4==\"C\" && $5>0.3) print $2}' ", 
        qc_dir, "/", batch, "_qc1.frq >> ", qc_dir, "/", batch, "_qc1_badsnps.txt"
      )
    )
    
    # remove badsnps
    system(
      paste0(
        plink1.9, 
        " --bfile ", qc_dir, "/", batch, "_qc1", 
        " --exclude ", qc_dir, "/", batch, "_qc1_badsnps.txt", 
        " --make-bed", 
        " --out ", qc_dir, "/", batch, "_qc2"
      )
    )
    
    ##### step 2 (individual filtering) #####
    
    ## include only individuals with call rate ≥ 0.98
    system(
      paste0(
        plink1.9, 
        " --bfile ", qc_dir, "/", batch, "_qc2", 
        " --mind 0.02", 
        " --make-bed", 
        " --out ", qc_dir, "/", batch, "_qc3"
      )
    )
    
    ## exclude individuals when FHET outside +/- 0.20
    system(
      paste0(
        plink1.9, 
        " --bfile ", qc_dir, "/", batch, "_qc3", 
        " --het", 
        " --out ", qc_dir, "/", batch, "_qc3"
      )
    )
    
    system(
      paste0(
        "awk '{if (($6 > 0.2) || ($6 < -0.2)) print $1,$2}' ", 
        qc_dir, "/", batch, "_qc3.het > ", qc_dir, "/", batch, "_qc3_fhet_badinds.txt"
      )
    )
    
    system(
      paste0(
        plink1.9, 
        " --bfile ", qc_dir, "/", batch, "_qc3", 
        " --remove ", qc_dir, "/", batch, "_qc3_fhet_badinds.txt", 
        " --make-bed", 
        " --out ", qc_dir, "/", batch, "_qc4"
      )
    )
    
    ## exclude sex violations - genetic sex does not match pedigree sex
    
    # add sex information if it's missing
    covars <- fread(here("data_raw", "cohorts", "demgene", "geno", "demgene.cov"))
    fam <- fread(paste0(qc_dir, "/", batch, "_qc4.fam"))
    
    fam_covars <- left_join(fam, covars, c("V1" = "FID", "V2" = "IID")) %>%
      mutate(V5 = case_when(
        V5 == 0 ~ sex,
        .default = V5
      )) %>%
      select(V1, V2, V3, V4, V5, V6)
    
    # overwrite the fam file
    fwrite(
      fam_covars,
      file = paste0(qc_dir, "/", batch, "_qc4.fam"),
      quote = F,
      sep = " ",
      col.names = F,
      row.names = F
    )
    
    system(
      paste0(
        plink1.9, 
        " --bfile ", qc_dir, "/", batch, "_qc4", 
        " --check-sex", 
        " --out ", qc_dir, "/", batch, "_qc4_sex"
      )
    )
    
    system(
      paste0(
        "grep PROBLEM ", qc_dir, "/", batch, "_qc4_sex.sexcheck > ", qc_dir, "/", batch, "_qc4_sex_badinds.txt"
      )
    )
    
    system(
      paste0(
        plink1.9, 
        " --bfile ", qc_dir, "/", batch, "_qc4", 
        " --remove ", qc_dir, "/", batch, "_qc4_sex_badinds.txt",
        " --make-bed",
        " --out ", qc_dir, "/", batch, "_qc5"
      )
    )
    
    ##### step 3 (SNP filtering) ##### 
    
    ## include only SNPs with call rate ≥ 0.98 
    ## (note: after filtering out “bad” individuals, we now apply a stricter call rate threshold)
    system(
      paste0(
        plink1.9, 
        " --bfile ", qc_dir, "/", batch, "_qc5", 
        " --geno 0.02", 
        " --make-bed", 
        " --out ", qc_dir, "/", batch, "_qc6"
      )
    )
    
    
    ## exclude invariant and multi-allelic SNPs
    system(
      paste0(
        plink2, 
        " --bfile ", qc_dir, "/", batch, "_qc6", 
        " --max-alleles 2", 
        " --make-bed", 
        " --out ", qc_dir, "/", batch, "_qc7"
      )
    )
    
    system(
      paste0(
        plink2, 
        " --bfile ", qc_dir, "/", batch, "_qc7", 
        " --min-alleles 2", 
        " --make-bed", 
        " --out ", qc_dir, "/", batch, "_qc8"
      )
    )
    
    ## include only SNPs with MAF ≥ 0.01
    system(
      paste0(
        plink1.9, 
        " --bfile ", qc_dir, "/", batch, "_qc8", 
        " --freq", 
        " --out ", qc_dir, "/", batch, "_qc8"
      )
    )
    
    system(
      paste0(
        "awk '{if ($5<0.01) print $2}' ", 
        qc_dir, "/", batch, "_qc8.frq > ", qc_dir, "/", batch, "_qc8_badsnps.txt"
      )
    )
    
    # remove badsnps
    system(
      paste0(
        plink1.9, 
        " --bfile ", qc_dir, "/", batch, "_qc8", 
        " --exclude ", qc_dir, "/", batch, "_qc8_badsnps.txt", 
        " --make-bed", 
        " --out ", qc_dir, "/", batch, "_qc9"
      )
    )
    
    ## include only SNPs where Hardy-Weinberg equilibrium (HWE) in controls p value ≥ 1e-06
    
    # subset to controls
    system(
      paste0(
        "awk 'NR > 1 {print $11,$12}' ", 
        here("data_raw/cohorts/demgene/demgene4pgc_export_04032024/demgene_redcap_control.txt"),
        " > ", qc_dir, "/", batch, "_qc9_controls.txt"
      )
    )
    
    system(
      paste0(
        "awk 'NR > 1 {print $3,$4}' ", 
        here("data_raw/cohorts/demgene/demgene4pgc_export_04032024/top_controls_redcap_age_sex_patient_status.txt"),
        " >> ", qc_dir, "/", batch, "_qc9_controls.txt"
      )
    )
    
    # keep only controls
    system(
      paste0(
        plink1.9, 
        " --bfile ", qc_dir, "/", batch, "_qc9", 
        " --keep ", qc_dir, "/", batch, "_qc9_controls.txt", 
        " --make-bed", 
        " --out ", qc_dir, "/", batch, "_qc9_controls"
      )
    )
    
    system(
      paste0(
        plink1.9, 
        " --bfile ", qc_dir, "/", batch, "_qc9_controls", 
        " --hardy", 
        " --out ", qc_dir, "/", batch, "_qc9_controls"
      )
    )
    
    # keep only cases
    system(
      paste0(
        plink1.9, 
        " --bfile ", qc_dir, "/", batch, "_qc9", 
        " --remove ", qc_dir, "/", batch, "_qc9_controls.txt", 
        " --make-bed", 
        " --out ", qc_dir, "/", batch, "_qc9_cases"
      )
    )
    
    system(
      paste0(
        plink1.9, 
        " --bfile ", qc_dir, "/", batch, "_qc9_cases", 
        " --hardy", 
        " --out ", qc_dir, "/", batch, "_qc9_cases"
      )
    )
    
    # bad SNPs in controls
    system(
      paste0(
        "awk 'NR > 1 {if ($9+0 < 0.000001) print $2}' ", 
        qc_dir, "/", batch, "_qc9_controls.hwe",
        " > ", qc_dir, "/", batch, "_qc9_badsnps.txt"
      )
    )
    
    # bad SNPs in cases
    system(
      paste0(
        "awk 'NR > 1 {if ($9+0 < 0.0000000001) print $2}' ", 
        qc_dir, "/", batch, "_qc9_cases.hwe",
        " >> ", qc_dir, "/", batch, "_qc9_badsnps.txt"
      )
    )
    
    # remove badsnps
    system(
      paste0(
        plink1.9, 
        " --bfile ", qc_dir, "/", batch, "_qc9", 
        " --exclude ", qc_dir, "/", batch, "_qc9_badsnps.txt", 
        " --make-bed", 
        " --out ", qc_dir, "/", batch, "_qc10"
      )
    )
    
#    ## remove X-chromosome SNPs
#    system(
#      paste0(
#        "awk 'NR > 1 {if ($1 == 23) print $2}' ", 
#        qc_dir, "/", batch, "_qc10.bim",
#        " > ", qc_dir, "/", batch, "qc10_chr23_badsnps.txt"
#      )
#    )
#    
#    system(
#      paste0(
#        plink1.9, 
#        " --bfile ", qc_dir, "/", batch, "_qc10", 
#        " --exclude ", qc_dir, "/", batch, "qc10_chr23_badsnps.txt", 
#        " --make-bed", 
#        " --out ", qc_dir, "/", batch, "_qc11"
#      )
#    )
    
    ## check for error messages
    log_files <- list.files(qc_dir, pattern = glob2rx("*log"), full.names = TRUE)
    results <- data.frame(file = character(), line = character(), stringsAsFactors = FALSE)
    for (file in log_files) {
      lines <- readLines(file)
      error_lines <- lines[str_detect(lines, regex("error|warning", ignore_case = TRUE))]
      
      if (length(error_lines) > 0) {
        temp_df <- data.frame(file = basename(file), line = error_lines, stringsAsFactors = FALSE)
        results <- bind_rows(results, temp_df)
      }
    }
    write.table(
      results,
      file = here("data_raw", "cohorts", "demgene", "qc", paste0(batch, "_error_warning.txt")),
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
    
    print(results)
     
  }
}



##### c. combine individual batches #####
if (FALSE) {
  
  ###### save batch and array info ###### 

  covars <- fread(here("data_raw", "cohorts", "demgene", "geno", "demgene.cov")) %>%
    filter(diagnosis %in% c("AD", "control"))

  summary_arrays <- covars %>%
    group_by(batch, array, diagnosis) %>%
    summarise(count = n(), .groups = 'drop') %>%
    tidyr::pivot_wider(names_from = diagnosis, values_from = count, values_fill = list(count = 0))
  
  # write the summary
  fwrite(
    summary_arrays,
    file = here("data_raw/cohorts/demgene/summary_info/batch_array_diagnosis.txt"),
    sep = "\t",
    quote = F,
    na = NA,
    row.names = F
  )
  
  #####  split by genotyping array  #####
  ## the SNPs genotyped differs between arrays, which causes merging to
  ## drastically reduce the genotyping rate
  # 1. OmniExpress
  # 2. DeCodeGenetics_V1_20012591_A1
  # 3. GSAv3
  
  batches_omni <- summary_arrays %>%
    filter(
      str_detect(array, "OmniExpress"),
      batch != "1706_batch14", # all individuals removed in qc
      batch != "2009_lbd8" # only 133 controls that causes a lot of merge problems
      ) %>%
    pull(batch)
  
  batches_decode <- summary_arrays %>%
    filter(
      str_detect(array, "DeCode"),
      batch != "1706_batch14" # all individuals removed in qc
      ) %>%
    pull(batch)
  
  batches_gsa <- summary_arrays %>%
    filter(
      str_detect(array, "GSA"),
      batch != "1706_batch14") %>%
    pull(batch)
  
  #####  1. OmniExpress #####  
  ## batches to merge
  batches_dirs <- covars %>%
    select(batch) %>%
    filter(batch %in% batches_omni) %>% # has no individuals after qc
    unique() %>%
    mutate(
      dir = paste0(here("data_raw", "cohorts", "demgene", "qc/"), batch, "/", batch, "_qc10")
    )
    
  ## save names of batches to be merged
  writeLines(
    sort(batches_dirs$dir), 
    here("data_raw/cohorts/demgene/qc/merged_files/demgene_OmniExpress_files_to_merge.txt")
    )
  
  ## merge batches
  files_to_merge <- here("data_raw/cohorts/demgene/qc/merged_files/demgene_OmniExpress_files_to_merge.txt")
  outfile <- here("data_raw/cohorts/demgene/qc/merged_files/demgene_OmniExpress")
  system(paste0(
    plink1.9,
    " --merge-list ", files_to_merge,
    " --make-bed --out ", outfile
  ))
  
  #####  2. DeCodeGenetics & GSA   #####  
  
  ## batches to merge
  batches_dirs <- covars %>%
    select(batch) %>%
    filter(batch %in% batches_decode) %>% # has no individuals after qc
    unique() %>%
    mutate(
      dir = paste0(here("data_raw", "cohorts", "demgene", "qc/"), batch, "/", batch, "_qc10")
    )
    
  writeLines(
    sort(batches_dirs$dir), 
    here("data_raw/cohorts/demgene/qc/merged_files/demgene_DeCodeGenetics_files_to_merge.txt")
  )
  
  dirs <- batches_dirs %>%
    pull(dir) %>%
    sort()
  
  batches <- batches_dirs %>%
    pull(batch) %>%
    sort()
  
  ## 1. merge (2 batches)

  # looks like strand issues
  i <- 1
  system(paste0(
    plink1.9,
    " --bfile ", dirs[i],
    " --bmerge ", dirs[i+1],
    " --make-bed", 
    " --out ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), "merged_decode", i+1, "_", batches[i+1])
  )
  
  # flip strands
  # duplicated SNPs which I need to remove first
  system(paste0(
    plink1.9,
    " --bfile ", dirs[i+1],
    " --flip ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), "merged_decode", i+1, "_", batches[i+1], "-merge.missnp",
    " --make-bed", 
    " --out ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), batches[i+1], "_flip"
    ))
  
  ## find and remove duplicates
  for (j in 1:length(dirs)) {
    bim <- fread(paste0(dirs[j], ".bim"))
    duplicates <- bim %>%
      group_by(across(all_of("V2"))) %>%
      filter(n() > 1) %>%
      ungroup() %>%
      unique()
    
    writeLines(
      sort(duplicates$V2), 
      paste0(here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), batches[j], "_duplicate_snp.txt")
    )
    
    system(paste0(
      plink1.9,
      " --bfile ", dirs[j],
      " --exclude ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), batches[j], "_duplicate_snp.txt",
      " --make-bed", 
      " --out ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), batches[j], "_no_duplicates"
    ))
  }

  # flip strands
  system(paste0(
    plink1.9,
    " --bfile ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), batches[i+1], "_no_duplicates",
    " --flip ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), "merged_decode", i+1, "_", batches[i+1], "-merge.missnp",
    " --make-bed", 
    " --out ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), batches[i+1], "_flip"
  ))
  
  # try merge again. no more errors
  system(paste0(
    plink1.9,
    " --bfile ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), batches[i], "_no_duplicates",
    " --bmerge ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), batches[i+1], "_flip",
    " --make-bed", 
    " --out ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), "merged_decode", i+1, "_", batches[i+1])
  )
  
  ## merge the other 4 batches
  for (i in 2:length(dirs)-1) {
    
    # strand issues again
    system(paste0(
      plink1.9,
      " --bfile ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), "merged_decode", i, "_", batches[i],
      " --bmerge ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), batches[i+1], "_no_duplicates",
      " --make-bed", 
      " --out ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), "merged_decode", i+1, "_", batches[i+1])
    )
    
    # flip strands
    system(paste0(
      plink1.9,
      " --bfile ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), batches[i+1], "_no_duplicates",
      " --flip ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), "merged_decode", i+1, "_", batches[i+1], "-merge.missnp",
      " --make-bed", 
      " --out ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), batches[i+1], "_flip"
    ))
    
    system(paste0(
      plink1.9,
      " --bfile ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), "merged_decode", i, "_", batches[i],
      " --bmerge ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), batches[i+1], "_flip",
      " --make-bed", 
      " --out ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), "merged_decode", i + 1, "_", batches[i + 1])
    )  
    
  }
  
  ## merge gsa batches
  gsa_batches_dirs <- covars %>%
    select(batch) %>%
    filter(batch %in% batches_gsa) %>% # has no individuals after qc
    unique() %>%
    mutate(
      dir = paste0(here("data_raw", "cohorts", "demgene", "qc/"), batch, "/", batch, "_qc10")
    )
  
  writeLines(
    sort(batches_dirs$dir), 
    here("data_raw/cohorts/demgene/qc/merged_files/demgene_GSA_files_to_merge.txt")
  )
  
  gsa_dirs <- gsa_batches_dirs %>%
    pull(dir) %>%
    sort()
  
  gsa_batches <- gsa_batches_dirs %>%
    pull(batch) %>%
    sort()

  ## find and remove duplicate variants
  for (j in 1:length(gsa_dirs)) {
    bim <- fread(paste0(gsa_dirs[j], ".bim"))
    duplicates <- bim %>%
      group_by(across(all_of("V2"))) %>%
      filter(n() > 1) %>%
      ungroup() %>%
      unique()
    
    writeLines(
      sort(duplicates$V2), 
      paste0(here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), gsa_batches[j], "_duplicate_snp.txt")
    )
    
    system(paste0(
      plink1.9,
      " --bfile ", gsa_dirs[j],
      " --exclude ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), gsa_batches[j], "_duplicate_snp.txt",
      " --make-bed", 
      " --out ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), gsa_batches[j], "_no_duplicates"
    ))
  }
  
  ## add 1. GSA batch
  system(paste0(
    plink1.9,
    " --bfile ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), "merged_decode", 6, "_", batches[6],
    " --bmerge ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), gsa_batches[1], "_no_duplicates",
    " --make-bed", 
    " --out ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), "merged_decode_gsa", 1, "_", gsa_batches[1])
  )
  
  # flip strands
  system(paste0(
    plink1.9,
    " --bfile ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), gsa_batches[1], "_no_duplicates",
    " --flip ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), "merged_decode_gsa", 1, "_", gsa_batches[1], "-merge.missnp",
    " --make-bed", 
    " --out ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), gsa_batches[1], "_flip"
  ))
  
  system(paste0(
    plink1.9,
    " --bfile ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), "merged_decode", 6, "_", batches[6],
    " --bmerge ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), gsa_batches[1], "_flip",
    " --make-bed", 
    " --out ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), "merged_decode_gsa", 1, "_", gsa_batches[1])
  )
  
  ## add 2. GSA batch
  system(paste0(
    plink1.9,
    " --bfile ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), "merged_decode_gsa", 1, "_", gsa_batches[1],
    " --bmerge ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), gsa_batches[2], "_no_duplicates",
    " --make-bed", 
    " --out ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), "merged_decode_gsa", 2, "_", gsa_batches[2])
  )
  
  # flip strands
  system(paste0(
    plink1.9,
    " --bfile ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), gsa_batches[2], "_no_duplicates",
    " --flip ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), "merged_decode_gsa", 2, "_", gsa_batches[2], "-merge.missnp",
    " --make-bed", 
    " --out ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), gsa_batches[2], "_flip"
  ))
  
  system(paste0(
    plink1.9,
    " --bfile ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), "merged_decode_gsa", 1, "_", gsa_batches[1],
    " --bmerge ", here("data_raw", "cohorts", "demgene", "qc", "temp_merge_files/"), gsa_batches[2], "_flip",
    " --make-bed", 
    " --out ", here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_DeCodeGenetics_GSA")
  )

  #### d. PCA plots ####
  
  library(ggplot2)
  library(cowplot)
  
  AD_controls <- fread(here("data_raw/cohorts/demgene/geno/demgene.cov")) %>%
    filter(diagnosis %in% c("AD", "control")) %>%
    select(FID, IID)
  
  fwrite(
    x = AD_controls,
    file = here(
      "data_raw/cohorts/demgene/qc/merged_files",
      "AD_control_ids.txt"
    ),
    quote = F,
    sep = "\t",
    col.names = F,
    row.names = F
  )
  
  for (batch in c("DeCodeGenetics_GSA", "OmniExpress")) {
    
    ## only keep cases and controls
    system(
      paste0(
        plink1.9, 
        " --bfile ", here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_", batch, 
        " --keep ", here("data_raw/cohorts/demgene/qc/merged_files", "AD_control_ids.txt"), 
        " --make-bed", 
        " --out ", here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_AD_control_", batch
      )
    )
    
    ## only keep SNPs with genotyping rate > 0.98
    system(
      paste0(
        plink1.9, 
        " --bfile ", here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_AD_control_", batch, 
        " --geno 0.02", 
        " --make-bed", 
        " --out ", here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_geno0.02_AD_control_", batch
      )
    )
    
    system(paste0(
      plink1.9,
      " --bfile ", here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_geno0.02_AD_control_", batch,
      " --genome",
      " --out ", here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_geno0.02_AD_control_", batch
    ))
    
    system(paste0(
      plink1.9,
      " --bfile ", here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_geno0.02_AD_control_", batch,
      " --cluster",
      " --mds-plot 20",
      " --out ", here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_geno0.02_AD_control_", batch
    ))
    
    
    ## make plot
    covars <- fread(here("data_raw/cohorts/demgene/geno/demgene.cov")) %>%
      select(FID, IID, batch, array)
    mds <- fread(paste0(here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_geno0.02_AD_control_", batch, ".mds")) %>%
      left_join(covars, by = c("FID", "IID"))
    
    p_batch <- mds %>%
      ggplot(aes(x = C1, y = C2, colour = batch)) +
      geom_point()
    
    p_array <- mds %>%
      ggplot(aes(x = C1, y = C2, colour = array)) +
      geom_point()
    
    p <- plot_grid(
      p_batch,
      p_array,
      labels = c('batches', 'arrays'),
      ncol = 1,
      align = "v"
    )
    
    ggsave(
      plot = p,
      paste0(
        here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "pca_demgene_geno0.02_AD_control_", batch, ".png"
      ),
      width = 22,
      height = 20,
      units = "cm"
    )
    
  }
  
  ## after inspection of the PCA plots I decided to 
  ### 1. isolate 1706_batch11
  
  ## 1706_batch11 (to be analysed later)
  covars <- fread(here("data_raw", "cohorts", "demgene", "geno", "demgene.cov"))
  batches_to_remove <- covars %>%
    filter(batch == "1706_batch11") %>%
    select(FID, IID)
  
  fwrite(
    x = batches_to_remove,
    file = here(
      "data_raw/cohorts/demgene/qc/merged_files",
      "1706_batch11_ids.txt"
    ),
    quote = F,
    sep = "\t",
    col.names = F,
    row.names = F
  )
  
  system(
    paste0(
      plink1.9, 
      " --bfile ", here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_geno0.02_AD_control_DeCodeGenetics_GSA", 
      " --remove ", here("data_raw/cohorts/demgene/qc/merged_files","1706_batch11_ids.txt"), 
      " --make-bed", 
      " --out ", here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_geno0.02_AD_control_DeCodeGenetics_GSA_no1706_batch11"
    )
  )
  
  system(paste0(
    plink1.9,
    " --bfile ", here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_geno0.02_AD_control_DeCodeGenetics_GSA_no1706_batch11",
    " --genome",
    " --out ", here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_geno0.02_AD_control_DeCodeGenetics_GSA_no1706_batch11"
  ))
  
  system(paste0(
    plink1.9,
    " --bfile ", here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_geno0.02_AD_control_DeCodeGenetics_GSA_no1706_batch11",
    " --cluster",
    " --mds-plot 20",
    " --out ", here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_geno0.02_AD_control_DeCodeGenetics_GSA_no1706_batch11"
  ))
  
  ## make plot
  mds <- fread(paste0(here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_geno0.02_AD_control_DeCodeGenetics_GSA_no1706_batch11", ".mds")) %>%
    left_join(covars, by = c("FID", "IID"))
  
  p_batch <- mds %>%
    ggplot(aes(x = C1, y = C2, colour = batch)) +
    geom_point()
  
  p_array <- mds %>%
    ggplot(aes(x = C1, y = C2, colour = array)) +
    geom_point()
  
  p <- plot_grid(
    p_batch,
    p_array,
    labels = c('batches', 'arrays'),
    ncol = 1,
    align = "v"
  )
  
  ggsave(
    plot = p,
    paste0(
      here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "pca_demgene_geno0.02_AD_control_DeCodeGenetics_GSA_no1706_batch11", ".png"
    ),
    width = 22,
    height = 20,
    units = "cm"
  )
  
  ## inspection of the plots identified 3 additional outliers which will also be removed
  
  mds <- fread(
    paste0(
      here("data_raw", "cohorts", "demgene", "qc/merged_files/"),
      "demgene_geno0.02_AD_control_DeCodeGenetics_GSA_no1706_batch11",
      ".mds"
    )
  )
  
  remove_outliers <- mds %>%
    filter(C2 > 0.1) %>%
    select(FID, IID)
  
  fwrite(
    x = remove_outliers,
    file = paste0(
      here("data_raw", "cohorts", "demgene", "qc/merged_files/"),
      "demgene_geno0.02_AD_control_DeCodeGenetics_GSA_no1706_batch11_outliers.txt"
    ),
    quote = F,
    sep = " ",
    col.names = F
  )
  
  system(
    paste0(
      plink1.9, 
      " --bfile ", here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_geno0.02_AD_control_DeCodeGenetics_GSA_no1706_batch11", 
      " --remove ", paste0(here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_geno0.02_AD_control_DeCodeGenetics_GSA_no1706_batch11_outliers.txt"), 
      " --make-bed", 
      " --out ", here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_geno0.02_nooutlier_AD_control_DeCodeGenetics_GSA_no1706_batch11"
    )
  )
  
  system(paste0(
    plink1.9,
    " --bfile ", here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_geno0.02_nooutlier_AD_control_DeCodeGenetics_GSA_no1706_batch11",
    " --genome",
    " --out ", here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_geno0.02_nooutlier_AD_control_DeCodeGenetics_GSA_no1706_batch11"
  ))
  
  system(paste0(
    plink1.9,
    " --bfile ", here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_geno0.02_nooutlier_AD_control_DeCodeGenetics_GSA_no1706_batch11",
    " --cluster",
    " --mds-plot 20",
    " --out ", here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_geno0.02_nooutlier_AD_control_DeCodeGenetics_GSA_no1706_batch11"
  ))
  
  ## make plot
  mds <- fread(paste0(here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "demgene_geno0.02_nooutlier_AD_control_DeCodeGenetics_GSA_no1706_batch11", ".mds")) %>%
    left_join(covars, by = c("FID", "IID"))
  
  p_batch <- mds %>%
    ggplot(aes(x = C1, y = C2, colour = batch)) +
    geom_point()
  
  p_array <- mds %>%
    ggplot(aes(x = C1, y = C2, colour = array)) +
    geom_point()
  
  p <- plot_grid(
    p_batch,
    p_array,
    labels = c('batches', 'arrays'),
    ncol = 1,
    align = "v"
  )
  
  ggsave(
    plot = p,
    paste0(
      here("data_raw", "cohorts", "demgene", "qc/merged_files/"), "pca_demgene_geno0.02_nooutlier_AD_control_DeCodeGenetics_GSA_no1706_batch11", ".png"
    ),
    width = 22,
    height = 20,
    units = "cm"
  )
  
  ##### e. save data #####
  
  ###### OmniExpress ###### 
  covars <- fread(here("data_raw", "cohorts", "demgene", "geno", "demgene.cov"))
  omni_fam <- fread(here("data_raw/cohorts/demgene/qc/merged_files", "demgene_geno0.02_AD_control_OmniExpress.fam")) %>%
    left_join(covars, by = c("V1" = "FID", "V2" = "IID")) %>%
    mutate(V6 = case_when(
      diagnosis == "AD" ~ 2,
      diagnosis == "control" ~ 1,
      .default = -9
    ))
  
  ## save .fam file with phenotype
  omni_fam %>%
    select(V1:V6) %>%
    fwrite(
      file = here(
        "data_raw/cohorts/demgene/geno",
        "demgene_OmniExpress.fam"
      ),
      quote = F,
      col.names = F,
      sep = " "
    )
  
  ## save .cov file
  omni_fam %>%
    select(V1, V2, diagnosis, age, sex, batch, array) %>%
    rename(
      FID = V1,
      IID = V2
    ) %>%
    fwrite(
      file = here(
        "data_raw/cohorts/demgene/geno",
        "demgene_OmniExpress.cov"
      ),
      quote = F,
      col.names = T,
      sep = " "
    )
  
  ## copy .bed and .bim file
  system(paste0(
    "cp ",
    here(
      "data_raw/cohorts/demgene/qc/merged_files",
      "demgene_geno0.02_AD_control_OmniExpress.bim ",
      here("data_raw/cohorts/demgene/geno/demgene_OmniExpress.bim")
    )
  ))
  system(paste0(
    "cp ",
    here(
      "data_raw/cohorts/demgene/qc/merged_files",
      "demgene_geno0.02_AD_control_OmniExpress.bed ",
      here("data_raw/cohorts/demgene/geno/demgene_OmniExpress.bed")
    )
  ))
  system(paste0(
    "cp ",
    here(
      "data_raw/cohorts/demgene/qc/merged_files",
      "demgene_geno0.02_AD_control_OmniExpress.hh ",
      here("data_raw/cohorts/demgene/geno/demgene_OmniExpress.hh")
    )
  ))
  
  ###### DeCodeGenetics_GSA_no1706_batch11 ###### 
  covars <- fread(here("data_raw", "cohorts", "demgene", "geno", "demgene.cov"))
  decode_fam <- fread(here("data_raw/cohorts/demgene/qc/merged_files", "demgene_geno0.02_nooutlier_AD_control_DeCodeGenetics_GSA_no1706_batch11.fam")) %>%
    left_join(covars, by = c("V1" = "FID", "V2" = "IID")) %>%
    mutate(V6 = case_when(
      diagnosis == "AD" ~ 2,
      diagnosis == "control" ~ 1,
      .default = -9
    ))
  
  ## save .fam file with phenotype
  decode_fam %>%
    select(V1:V6) %>%
    fwrite(
      file = here(
        "data_raw/cohorts/demgene/geno",
        "demgene_DeCodeGenetics_GSA_no1706_batch11.fam"
      ),
      quote = F,
      col.names = F,
      sep = " "
    )
  
  ## save .cov file
  decode_fam %>%
    select(V1, V2, diagnosis, age, sex, batch, array) %>%
    rename(
      FID = V1,
      IID = V2
    ) %>%
    fwrite(
      file = here(
        "data_raw/cohorts/demgene/geno",
        "demgene_DeCodeGenetics_GSA_no1706_batch11.cov"
      ),
      quote = F,
      col.names = T,
      sep = " "
    )
  
  ## copy .bed and .bim file
  system(paste0(
    "cp ",
    here(
      "data_raw/cohorts/demgene/qc/merged_files",
      "demgene_geno0.02_nooutlier_AD_control_DeCodeGenetics_GSA_no1706_batch11.bim ",
      here("data_raw/cohorts/demgene/geno/demgene_DeCodeGenetics_GSA_no1706_batch11.bim")
    )
  ))
  system(paste0(
    "cp ",
    here(
      "data_raw/cohorts/demgene/qc/merged_files",
      "demgene_geno0.02_nooutlier_AD_control_DeCodeGenetics_GSA_no1706_batch11.bed ",
      here("data_raw/cohorts/demgene/geno/demgene_DeCodeGenetics_GSA_no1706_batch11.bed")
    )
  ))
  system(paste0(
    "cp ",
    here(
      "data_raw/cohorts/demgene/qc/merged_files",
      "demgene_geno0.02_nooutlier_AD_control_DeCodeGenetics_GSA_no1706_batch11.hh ",
      here("data_raw/cohorts/demgene/geno/demgene_DeCodeGenetics_GSA_no1706_batch11.hh")
    )
  ))
  
  
  ###### DeCodeGenetics 1706_batch11 ###### 
  covars <- fread(here("data_raw", "cohorts", "demgene", "geno", "demgene.cov"))
  
  cases_controls <- fread(here("data_raw/cohorts/demgene/qc/1706_batch11", "1706_batch11_qc10.fam")) %>%
    left_join(covars, by = c("V1" = "FID", "V2" = "IID")) %>%
    filter(diagnosis %in% c("AD", "control")) %>%
    select(V1, V2)
  
  fwrite(
    x = cases_controls,
    file = here(
      "data_raw/cohorts/demgene/qc/1706_batch11",
      "1706_batch11_qc10_case_control_ids.txt"
    ),
    quote = F,
    sep = " ",
    col.names = F
  )
  
  ## keep only AD cases and controls
  system(
    paste0(
      plink1.9, 
      " --bfile ", here("data_raw", "cohorts", "demgene", "qc/1706_batch11", "1706_batch11_qc10"), 
      " --keep ", here("data_raw/cohorts/demgene/qc/1706_batch11","1706_batch11_qc10_case_control_ids.txt"), 
      " --make-bed", 
      " --out ", here("data_raw", "cohorts", "demgene", "qc/1706_batch11", "1706_batch11_qc10_noMCI")
    )
  )
  
  fam <- fread(here("data_raw/cohorts/demgene/qc/1706_batch11", "1706_batch11_qc10_noMCI.fam")) %>%
    left_join(covars, by = c("V1" = "FID", "V2" = "IID")) %>%
    mutate(V6 = case_when(
      diagnosis == "AD" ~ 2,
      diagnosis == "control" ~ 1,
      .default = -9
    ))
  
  ## save .fam file with phenotype
  fam %>%
    select(V1:V6) %>%
    fwrite(
      file = here(
        "data_raw/cohorts/demgene/geno",
        "demgene_1706_batch11.fam"
      ),
      quote = F,
      col.names = F,
      sep = " "
    )
  
  ## save .cov file
  fam %>%
    select(V1, V2, diagnosis, age, sex, batch, array) %>%
    rename(
      FID = V1,
      IID = V2
    ) %>%
    fwrite(
      file = here(
        "data_raw/cohorts/demgene/geno",
        "demgene_1706_batch11.cov"
      ),
      quote = F,
      col.names = T,
      sep = " "
    )
  
  ## copy .bed and .bim file
  system(paste0(
    "cp ",
    here(
      "data_raw/cohorts/demgene/qc/1706_batch11",
      "1706_batch11_qc10_noMCI.bim ",
      here("data_raw/cohorts/demgene/geno/demgene_1706_batch11.bim")
    )
  ))
  system(paste0(
    "cp ",
    here(
      "data_raw/cohorts/demgene/qc/1706_batch11",
      "1706_batch11_qc10_noMCI.bed ",
      here("data_raw/cohorts/demgene/geno/demgene_1706_batch11.bed")
    )
  ))
  system(paste0(
    "cp ",
    here(
      "data_raw/cohorts/demgene/qc/1706_batch11",
      "1706_batch11_qc10_noMCI.hh ",
      here("data_raw/cohorts/demgene/geno/demgene_1706_batch11.hh")
    )
  ))

}