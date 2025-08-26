meta_analysis_f <- function(sumstats, out_dir, outfile, scheme = "STDERR", cols = "beta") {
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  ### run genome-wide association meta analysis with metal
  ### ref: Willer et al. 2010 (Bioinformatics)
  
  ## - sumstats:        vector of full sumstats file names
  ## - out_dir:         output sumstats
  ## - outfile:         file name of output
  ## - scheme:          standard error- or N-weighted meta analysis
  ##                    default is standard error weighted
  ##                    {SAMPLESIZE | STDERR}
  ## - cols:            column headers in the sumstats
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  
  #### 0. set-up (load libraries and install when necessary) ####
  if (!require(dplyr)) { install.packages("dplyr"); library(dplyr) }
  if (!require(data.table)) { install.packages("data.table"); library(data.table) }
  
  #### 1. run case_control and proxy meta-analysis ####
  for (proxy_cc in c("case_control", "proxy")) {
    
    temp_sumstats <- sumstats[grepl(pattern = proxy_cc, x = sumstats)]
    
    ## case_control and proxy meta-analysis input has beta columns
    if (scheme == "SAMPLESIZE") { cols <- "beta" }
    
    ## copy sumstats to current directory (scratch disk)
    system(paste0("rsync -av ", out_dir, "/", proxy_cc, " ."))
    
    ancestries <- unique(sub(".*_(eur|eas|afr|amr|sas).*", "\\1", temp_sumstats))
    for (ancestry in c(ancestries, "all")) {
      
      full_outfile <- paste0(outfile, "_", proxy_cc, "_", ancestry)
      
      ## define options for metal
      if (proxy_cc == "case_control") {
        opts <- paste(
          "CUSTOMVARIABLE neff", 
          "LABEL neff as neff",
          "CUSTOMVARIABLE n_case",
          "LABEL n_case as n_case",
          "CUSTOMVARIABLE n_control",
          "LABEL n_control as n_control", 
          sep = "\n")
      } else if (proxy_cc == "proxy" & ancestry != "all") {
        opts <- paste(
          "CUSTOMVARIABLE neff", 
          "LABEL neff as neff",
          "CUSTOMVARIABLE n_proxy_case",
          "LABEL n_proxy_case as n_case",
          "CUSTOMVARIABLE n_proxy_control",
          "LABEL n_proxy_control as n_control", 
          sep = "\n")
      } else if (proxy_cc == "proxy" & ancestry == "all") {
        opts <- paste(
          "CUSTOMVARIABLE neff", 
          "LABEL neff as neff",
          "CUSTOMVARIABLE n_proxy_case",
          "LABEL n_proxy_case as n_proxy_case",
          "CUSTOMVARIABLE n_proxy_control",
          "LABEL n_proxy_control as n_proxy_control", 
          sep = "\n")
      }
      
      if (ancestry == "all") {
        if (scheme == "SAMPLESIZE") { cols <- "z_score" }
        
        sumstats_subset <- NULL
        for (temp_ancestry in ancestries) {
          sumstats_subset <- c(sumstats_subset, paste0(proxy_cc, "/", outfile, "_", proxy_cc, "_", temp_ancestry, ".txt.gz"))
        }
        
      } else {
        ## subset to proxy_cc and ancestry in the loop
        sumstats_subset <- temp_sumstats[grepl(ancestry, temp_sumstats)]
      }
      
      ## add sample overlap correction for N-weighted meta-analysis 
      ## (note: leads to an extreme power loss, hence commented out)
      # if (scheme == "SAMPLESIZE") {
      #    opts <- paste(opts, "OVERLAP ON", sep = "\n")
      # }
      
      cat("attempting to run meta-analysis for the following files: \n.\n.\n.\n")
      print(sumstats_subset)
      cat(".\n.\n.\n")
      
      if (!file.exists(paste0(proxy_cc, "/", full_outfile, ".txt.gz"))) {
        
        ## run metal
        run_metal(
          sumstats = sumstats_subset,
          outfile = full_outfile,
          out_dir = paste0("./", proxy_cc),
          metal_dir = "METAL-2020-05-05/build/bin/metal",
          rsid_file = "g1000_afr.frq",
          opts = opts,
          cols = cols,
          scheme = scheme        
        )
        
        #### copy folder back to home directory ####
        system(paste0("rsync -av ", proxy_cc, " ", out_dir, "/"))
        
      } else {
        print(paste0("file already exists: ", proxy_cc, "/", full_outfile, ".txt.gz"))
      }
    }
    
    ##### add missing column #####
    if (proxy_cc == "case_control") {
      
      sumstats_cols = list.files(
        path = proxy_cc,
        full.names = T,
        pattern = "\\.gz$"
      )
      
      # Loop through each file
      for (file in sumstats_cols) {
        
        # Read the file with only column names (header only)
        current_cols <- names(fread(file, nrows = 0))
        
        # Initialize an empty list to store new columns
        new_cols <- list()
        
        # Add n_proxy_case column if it doesn't exist
        if (!"n_proxy_case" %in% current_cols) {
          new_cols$n_proxy_case <- 0
        }
        
        # Add n_proxy_control column if it doesn't exist
        if (!"n_proxy_control" %in% current_cols) {
          new_cols$n_proxy_control <- 0
        }
        
        # If any new columns were added, read the data and modify it
        if (length(new_cols) > 0) {
          data <- fread(file)
          data <- data %>% mutate(!!!new_cols)
          
          # Write the modified data back to the file
          fwrite(
            data,
            file = file,
            sep = "\t",
            na = "NA",
            quote = FALSE
          )
          
          print(paste0("Columns added to: ", file))
          
        } else {
          
          print(paste0("No columns needed to be added to: ", file))
          
        }
      }
    } else {
      
      sumstats_cols = list.files(
        path = proxy_cc,
        full.names = T,
        pattern = "\\.gz$"
      )
      
      # Loop through each file
      for (file in sumstats_cols) {
        
        # Read the file with only column names (header only)
        current_cols <- names(fread(file, nrows = 0))
        
        if (!"n_proxy_case" %in% current_cols | !"n_proxy_control" %in% current_cols) {
          fread(file) %>%
            rename(
              n_proxy_case = n_case,
              n_proxy_control = n_control) %>%
            mutate(
              n_case = 0,
              n_control = 0
            ) %>%
            fwrite(
              file = file,
              sep = "\t",
              na = "NA",
              quote = FALSE
            )
          print(paste0("Columns added to: ", file))
        } else if (!"n_case" %in% current_cols | !"n_control" %in% current_cols) {
          fread(file) %>%
            mutate(
              n_case = 0,
              n_control = 0
            ) %>%
            fwrite(
              file = file,
              sep = "\t",
              na = "NA",
              quote = FALSE
            )
          print(paste0("Columns added to: ", file))
        } else {
          print(paste0("No columns needed to be added to: ", file))
        }
      }
    }
    
    #### copy folder back to home directory ####
    system(paste0("rsync -av ", proxy_cc, " ", out_dir, "/"))
    
  }
  
  #### 2. run combined (case_control + proxy) meta-analysis #####
  ## make folder to store results
  system(paste0("rsync -av ", out_dir, "/proxy ."))
  system(paste0("rsync -av ", out_dir, "/case_control ."))
  system(paste0("rsync -av ", out_dir, "/combined ."))
  
  ## define options for metal
  opts <- paste(
    "CUSTOMVARIABLE neff", 
    "LABEL neff as neff",
    "CUSTOMVARIABLE n_proxy_case",
    "LABEL n_proxy_case as n_proxy_case",
    "CUSTOMVARIABLE n_proxy_control",
    "LABEL n_proxy_control as n_proxy_control", 
    "CUSTOMVARIABLE n_case",
    "LABEL n_case as n_case",
    "CUSTOMVARIABLE n_control",
    "LABEL n_control as n_control", 
    sep = "\n")
  
  ## add sample overlap correction for N-weighted meta-analysis
  if (scheme == "SAMPLESIZE") {
    # leads to an extreme power loss, hence commented out
    # opts <- paste(opts, "OVERLAP ON", sep = "\n")
    
    # combined meta-analysis input has z_score columns when N-weighted meta-
    # analysis was run
    cols <- "z_score"
  }
  
  ancestries = c(sub(".*_(eur|eas|afr|amr|sas|all).*","\\1", list.files("./case_control", pattern = "\\.gz$")),
                 sub(".*_(eur|eas|afr|amr|sas|all).*","\\1", list.files("./proxy", pattern = "\\.gz$")))
  ancestries <- unique(ancestries[ancestries != "all"])
  
  ## loop through ancestries, but make sure "all" is run last
  for (ancestry in c(ancestries, "all")) {
    
    full_outfile <- paste0(outfile, "_combined_", ancestry)
    
    if (!file.exists(paste0("combined/", full_outfile, ".txt.gz"))) {
      
      ## subset to proxy_cc and ancestry in the loop
      if (ancestry != "all") {
        sumstats_subset <- c(
          paste0("case_control/", outfile, "_case_control_", ancestry, ".txt.gz"),
          paste0("proxy/", outfile, "_proxy_", ancestry, ".txt.gz")
        )
      } else {
        sumstats_subset <- NULL
        for (temp_ancestry in ancestries) {
          sumstats_subset <- c(
            sumstats_subset,
            paste0("combined/", outfile, "_combined_", temp_ancestry, ".txt.gz")
          )
        }
      }
      
      ## run metal
      run_metal(
        sumstats = sumstats_subset,
        outfile = full_outfile,
        out_dir = "./combined",
        metal_dir = "METAL-2020-05-05/build/bin/metal",
        rsid_file = "g1000_afr.frq",
        opts = opts,
        cols = cols,
        scheme = scheme
      )
      
    }
  }
  
  ##### copy folder back to home directory #####
  system(paste0("rsync -av combined ", out_dir, "/"))
}
