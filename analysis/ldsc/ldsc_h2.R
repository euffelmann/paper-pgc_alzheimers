#### 0. setup ####
library(data.table)
library(dplyr)
library(here)
library(gsubfn)
library(openxlsx)
setwd(here("analysis/ldsc/ldsc_h2"))

#### 1. run ldsc for meta-analysis sumstats ####
K <- 0.05
for (analysis in c("main", "fema", "male", "noag", "apoe", "agex")) {
  for (anc in c("eur", "afr", "eas", "sas", "amr")) {
    for (proxy_casec in c("proxy", "case_control", "combined")) {
      
      ## ldsc ancestries
      if (anc == "sas") {
        anc_ldsc <- "CSA"
      } else {
        anc_ldsc <- toupper(anc)
      }
      
      pheno <- paste(analysis, proxy_casec, anc, "yeaapoe", sep = "_")
      sst_munged <- here("analysis/ldsc/ldsc_munge/meta", analysis, proxy_casec, paste0(pheno, "_munge.sumstats.gz"))
      
      if (file.exists(sst_munged) & !file.exists(here("analysis/ldsc/ldsc_h2/meta", paste0(pheno, "_h2.log")))) {
        
        system(
          paste0(
            "sbatch ", here("src/ldscH2.sh"), " ",
            "-f ", sst_munged, " ",
            "-e ", pheno, " ",
            "-o ", here("analysis/ldsc/ldsc_h2/meta"), " ",
            ## sample prevalence is set to  0.5 because I use neff as input to ldsc
            "-h '--pop-prev ", K, " --samp-prev 0.5' ",
            "-l ", here("data/reference_data/UKBB.ALL.ldscore/UKBB."), anc_ldsc, ".rsid "
          )
        )
      }
    }
  }
}


#### 2. collect results in excel table ####
ind <- 0
h2_df <- data.frame(proxy_casec = as.character())
K <- 0.05

params <- readWorkbook(here("analysis/cohort_summary/cohort_summary.xlsx"), sheet = "main") %>%
  mutate(
    subdir = ifelse(proxy_casecontrol == "proxy", "proxy_meta", proxy_casecontrol),
    anc_ldsc = ifelse(ancestry == "sas", "csa", ancestry)
  ) %>%
  mutate(
    sumstats_file = paste0(
      here("data/sumstats/"),
      subdir,
      "/clean_",
      cohort,
      "_",
      proxy_casecontrol,
      "_",
      analysis,
      "_",
      ancestry,
      ".txt.gz"
    )
  )

for (analysis in c("main", "fema", "male", "noag", "apoe", "agex")) {
  for (anc in c("eur", "afr", "eas", "sas", "amr")) {
    for (proxy_casec in c("proxy", "case_control", "combined")) { 
      for (apoe in c("yeaapoe")) {
        
        outfile <- paste(analysis, proxy_casec, anc, sep = "_")
        print(paste0("extracting results for: ", outfile))
        
        if (file.exists(paste0(here("analysis/ldsc/ldsc_h2/meta/"), outfile, "_", apoe, "_h2.log"))) {
          
          ## pull metrics out of ldsc log files
          h2 <- as.numeric(system(
            paste0( "grep 'Total Liability scale h2' ", here("analysis/ldsc/ldsc_h2/meta/"), outfile, "_", apoe, "_h2.log | awk '{print $5}'"),
            intern = T
          ))
          h2_se <- system(
            paste0( "grep 'Total Liability scale h2' ", here("analysis/ldsc/ldsc_h2/meta/"), outfile, "_", apoe, "_h2.log | awk '{print $6}'"),
            intern = T
          )
          h2_se <- as.numeric(gsubfn(".", list("(" = "", ")" = ""), h2_se))
          
          lambda_gc <- as.numeric(system(
            paste0("grep 'Lambda GC' ",here("analysis/ldsc/ldsc_h2/meta/"),outfile,"_",apoe,"_h2.log | awk '{print $3}'"),
            intern = T
          ))
          
          mean_chi2 <- as.numeric(system(
            paste0( "grep 'Mean Chi^2' ", here("analysis/ldsc/ldsc_h2/meta/"), outfile, "_", apoe, "_h2.log | awk '{print $3}'"),
            intern = T
          ))
          
          intercept <- as.numeric(system(
            paste0("grep 'Intercept' ",here("analysis/ldsc/ldsc_h2/meta/"),outfile,"_",apoe,"_h2.log | awk '{print $2}'"),
            intern = T
          ))
          
          intercept_se <- system(
            paste0("grep 'Intercept' ",here("analysis/ldsc/ldsc_h2/meta/"),outfile,"_",apoe,"_h2.log | awk '{print $3}'"),
            intern = T)
          intercept_se <- as.numeric(gsubfn(".", list("(" = "", ")" = ""), intercept_se))
          
          # ratio = (intercept-1)/(mean(chi^2)-1); this measures bias in mean^2
          ratio <- as.numeric(system(
            paste0("grep 'Ratio' ",here("analysis/ldsc/ldsc_h2/meta/"),outfile,"_",apoe,"_h2.log | awk '{print $2}'"),
            intern = T
          ))
          ratio_se <- system(
            paste0("grep 'Ratio' ",here("analysis/ldsc/ldsc_h2/meta/"),outfile,"_",apoe,"_h2.log | awk '{print $3}'"),
            intern = T)
          ratio_se <- as.numeric(gsubfn(".", list("(" = "", ")" = ""), ratio_se))   
          
          if (proxy_casec == "combined") {
            neff <- sum(params$neff[params$anc == anc])
          } else {
            neff <- sum(params$neff[params$anc == anc & params$proxy_casec == proxy_casec])
          }
          
          ## save metrics
          ind <- ind + 1
          h2_df[ind, "proxy_casec"] <- proxy_casec
          h2_df[ind, "analysis"] <- analysis
          h2_df[ind, "anc"] <- anc
          h2_df[ind, "apoe"] <- apoe
          h2_df[ind, "K"] <- as.numeric(K)
          h2_df[ind, "neff"] <- neff
          h2_df[ind, "h2"] <- h2
          h2_df[ind, "h2_se"] <- h2_se
          h2_df[ind, "lambda_gc"] <- lambda_gc
          h2_df[ind, "mean_chi2"] <- mean_chi2
          h2_df[ind, "intercept"] <- intercept
          h2_df[ind, "intercept_se"] <- intercept_se
          h2_df[ind, "ratio"] <- ratio
          h2_df[ind, "ratio_se"] <- ratio_se
          
        }
      }
    }
  }
}

h2_df %>%
  mutate(
    analysis = factor(analysis, levels = c("main", "noag", "agex", "fema", "male")),
    proxy_casec = factor(proxy_casec, levels = c("combined", "case_control", "proxy")),
    anc = factor(anc, levels = c("eur", "afr", "eas", "sas", "amr")),
    p = 2 * pnorm(q = abs(h2 / h2_se), lower.tail = FALSE),
  ) %>%
  relocate(
    p, .after = h2_se
  ) %>%
  arrange(analysis, proxy_casec, anc, apoe) %>%
  fwrite(file = here("analysis/ldsc/ldsc_h2.txt"), quote = F, sep = "\t", na = NA)


#### 3. run ldsc with --chisq-max 999999 ####


# , "fema", "male", "noag", "apoe", "agex"
# , "afr", "eas", "sas", "amr"
# "proxy", "case_control", 

K <- 0.05
for (analysis in c("main")) {
  for (anc in c("eur")) {
    for (proxy_casec in c("combined")) {
      
      ## ldsc ancestries
      if (anc == "sas") {
        anc_ldsc <- "CSA"
      } else {
        anc_ldsc <- toupper(anc)
      }
      
      pheno <- paste(analysis, proxy_casec, anc, "yeaapoe", sep = "_")
      sst_munged <- here("analysis/ldsc/ldsc_munge/meta", analysis, proxy_casec, paste0(pheno, "_munge.sumstats.gz"))
      
      if (file.exists(sst_munged) & !file.exists(here("analysis/ldsc/ldsc_h2/chi_sq_max", paste0(pheno, "_h2.log")))) {
        
        system(
          paste0(
            "sbatch ", here("src/ldscH2.sh"), " ",
            "-f ", sst_munged, " ",
            "-e ", pheno, " ",
            "-o ", here("analysis/ldsc/ldsc_h2/chi_sq_max"), " ",
            ## sample prevalence is set to  0.5 because I use neff as input to ldsc
            "-h '--pop-prev ", K, " --samp-prev 0.5 --chisq-max 999999' ",
            "-l ", here("data/reference_data/UKBB.ALL.ldscore/UKBB."), anc_ldsc, ".rsid "
          )
        )
      }
    }
  }
}