prs_r2obs_to_r2liab <- function(K,P,prs_r2obs){
  ## Lee et al. 2012 Genet Epidemiology
  t = -qnorm(K,mean=0,sd=1) # disease threshold
  z<-dnorm(t)               # height of the normal distribution at T
  i1<-z/K                   # mean liability of A1 (eg Falconer and Mackay) 
  k1=i1*(i1-t)              # reduction in variance in A1
  i0<--z/(1-K)              # mean liability of A0
  k0=i0*(i0-t)              # reduction in variance in A0
  
  theta = i1*(P-K)/(1-K)*(i1*(P-K)/(1-K)-t) # theta in equation (15) Lee et al. 2012 Genet Epidemiology
  cv = K*(1-K)/z^2*K*(1-K)/(P*(1-P))        # C in equation (15)
  R2 = prs_r2obs*cv/(1+prs_r2obs*theta*cv)
  
  return(R2)
}

prs_r2liab_to_r2obs <- function(K,P,prs_r2liab){
  ## Lee et al. 2012 Genet Epidemiology
  t = -qnorm(K,mean=0,sd=1) # disease threshold
  z<-dnorm(t)               # height of the normal distribution at T
  i1<-z/K                   # mean liability of A1 (eg Falconer and Mackay) 
  k1=i1*(i1-t)              # reduction in variance in A1
  i0<--z/(1-K)              # mean liability of A0
  k0=i0*(i0-t)              # reduction in variance in A0
  
  theta = i1*(P-K)/(1-K)*(i1*(P-K)/(1-K)-t) # theta in equation (15) Lee et al. 2012 Genet Epidemiology
  cv = K*(1-K)/z^2*K*(1-K)/(P*(1-P))        # C in equation (15)
  
  prs_r2obs = prs_r2liab / {cv-prs_r2liab*theta*cv}
  
  return(prs_r2obs)
}