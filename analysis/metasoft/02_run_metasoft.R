#### 0. set-up ####
library(here)


#### 1. run metasoft ####

for (proxy_casec in c("case_control", "combined", "proxy")) {
  
  ## input
  metasoft <- here("src/Metasoft/Metasoft.jar")
  pvalue_table <- here("src/Metasoft/HanEskinPvalueTable.txt")
  input <- here(paste0("analysis/metasoft/input/", proxy_casec, "_neff0.6_nsumstats1.txt"))
  output <- here(paste0("analysis/metasoft/output/", proxy_casec, "_out.txt"))
  log <- here(paste0("analysis/metasoft/output/", proxy_casec, "_log.txt"))

  ## run
  system(paste0(
    "java -jar ", metasoft, 
    " -pvalue_table ", pvalue_table, 
    " -input ", input,
    " -output ", output,
    " -log ", log))  
}
