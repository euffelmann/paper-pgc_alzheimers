neff_f <- function(n_case, n_control, proxy_casecontrol = "case_control") {
  ### function to calculate the effective sample size.
  ### proxy-case designs yield 1/4th the effective sample size of true-case 
  ### designs, therefore their neffs are divided by 4
  ### see (Liu et al., 2017; 10.1038/ng.3766)
  
  neff <- round(4 / ((1 / n_case) + (1 / n_control)))
  
  ## if proxy-case design, divide by 4
  if (proxy_casecontrol == "proxy") {
    neff <- neff / 4
  }
  return(neff)
}
