NR <- function(para, map, family, data, ref, V, bet0, sample.info, outcome, type, silent,
               learning.rate=1, nRestarts=5, step.scale=0.5) {

  inv.V        <- solve(V)       ##### Pass inv.V into NR_main #####
  step         <- learning.rate
  flag         <- 0

  # First check that para defines a feasible point. If not, then redefine para
  para <- checkInitParms(para, map, data, ref, inv.V, bet0, outcome, sample.info, family) 

  for (i in 1:nRestarts) {
    ret <- try(NR_main(para, map, family, data, ref, inv.V, bet0, sample.info, outcome, type, silent,
                       learning.rate=step), silent=TRUE)
    if ((is.null(ret)) || ("try-error" %in% class(ret))) {
      step <- step*step.scale
    } else {
      flag <- 1
      break
    }
  }
 
  if (!flag) stop("NR did not converge")

  #print(c(ret$conv, ret$iter))

  ret

} # END NR

checkInitParms <- function(para, map, data, ref, inv.V, bet0, outcome, sample.info, family) {

  fam_gaussian <- family == 'gaussian'
  fam_binomial <- family == 'binomial'
  fam_cc       <- family == 'case-control'
  fam_cml      <- family == 'cml'

  if (fam_gaussian) {
    f0 <- obj.lm(para, map, data, ref, inv.V, bet0, outcome)
  } else if (fam_binomial) {
    f0 <- obj.lo(para, map, data, ref, inv.V, bet0, outcome)
  } else if (fam_cc) {
    f0 <- obj.cc(para, map, data, ref, inv.V, bet0, sample.info, outcome) 
  } else if (fam_cml) {
    f0 <- obj.cml(para, map, data, ref, inv.V, bet0, sample.info, outcome)
  }
  
  if (!is.finite(f0)) {
    # For gaussian, likely cause is log(pr), so set lam values to 0
    if (fam_gaussian) {
      para[map$lam] <- 0
    } else {
      stop("ERROR: initial estimates do not define a feasible point")
    }
  }

  para

} # END: checkInitParms

# Newton-Raphson algorithm
NR_main <- function(para, map, family, data, ref, inv.V, bet0, sample.info, outcome, type, silent,
                    learning.rate=1){
  
  #inv.V     <- solve(V)    # now passed in instead of V
  np        <- length(para)
  para.null <- rep(NA, np)
  
  i <- 0
  while(i<1000){
    #fn <- compute.obj(para, map, family, data, ref, inv.V, bet0, sample.info, outcome, type)
    
    s0 <- compute.score(para, map, family, data, ref, inv.V, bet0, sample.info, outcome, type)
    if (any(!is.finite(s0))) return(NULL)
    
    if(!silent) cat('iter = ', i+1, '\t', formatC(max(abs(s0)), digits = 2, format = 'e'), '           \r')
    
    if(all(abs(s0) < 1e-6)){
      break
    }
    
    h0 <- compute.hess(para, map, family, data, ref, inv.V, bet0, sample.info, outcome, type)
    
    t0 <- try(d0 <- solve(h0, s0), silent = TRUE)
    if('try-error' %in% class(t0)){
      stop('hess fails')
    }
    #d0 <- d0/sqrt(sum(p^2))
    
    para <- para - learning.rate*d0
    
    i <- i + 1
    
  }
  
  if(all(abs(s0) > 1e-6)){
    stop('NR does not converge')
  }
  
  list(coefficients = para, score = s0, conv = ifelse(all(abs(s0) < 1e-6), 1, 0))
  
}
