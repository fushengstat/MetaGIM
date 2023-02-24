# main procedure for gim 
metagim <- function(models,group,sample.info,form,family,data){
  ncase = sample.info$ncase
  nctrl = sample.info$nctrl

  nmodel = length(group)
  I0 <- rep(list(NULL), nmodel)
  J0 <- rep(list(NULL), nmodel)
  npara <- rep(list(NULL), nmodel)
  thetas <- rep(list(NULL), nmodel)
  V.opt <- rep(list(NULL), nmodel)
  
  n.para <- NULL
  score = list()
  nmap <- list()
  outcome <- list()
  
  V.the <- list()
  
  for (i in 1:nmodel){
    id = group[[i]]
    nid = length(id)
    md = list()
    for (j in 1:nid){
      idj = id[j]
      md[[j]] = models[[idj]]
    }
    nca = ncase[id,id,drop = FALSE]
    nct = nctrl[id,id,drop = FALSE]
    
    samplei.info = list(ncase = nca, nctrl = nct)
    
    
    fit = gimfit(form, family, data, md, samplei.info)
    
    I0[[i]] <- fit$Iv0
    J0[[i]] <- fit$Jv0
    thetas[[i]] <- fit$coefficients
    npara[[i]] <- fit$parms
    V.the[[i]] <- fit$vcov
    
    nmap[[i]] <- fit$map
    V.opt[[i]] <- fit$V.bet
    outcome[[i]] <- fit$outcome

    n.para[i] = length(fit$parms)
    score[[i]] = fit$score
  }
  
  ### compute new Sig0 for all summary data
  res0 <- organize(thetas,npara,nmap,group,sample.info)
  para = res0$para
  map = res0$map
  K = res0$K
  
  ## pr0,Delta based on the internal-based mle
  mf0 <- model.frame(form, data = data)
  outcome0 <- colnames(mf0)[1]
  n1 <- sum(data[, outcome0])
  n0 <- sum(1 - data[, outcome0])	
  fit0 <- glm(form, family = 'binomial', data = data)	  
  pr0 <- (1-fit0$fitted.values)/n0
  # sum(pr0) == 1
  Delta <- exp(fit0$linear.predictors) * n0/n1
  
  fp <- formula.parse(form, family, data, models, ref = NULL)
  ref <- fp$ref
  Sig0 = Sigma0(para, map, ref, K, sample.info, pr0, Delta)
  
  
  
  ### Update the Covariance matrix for all parameters
  lam = para[1]
  daty = data[, outcome0]
  Inew = Cov_Omega(npara,nmap,score,daty,lam,n.para,J0,Sig0)
 
  
  Vthe = Varcof_meta(nmap,J0,Inew)
  nthetas = theta_refine(thetas)
  ### One-step meta-analysis for theta
  # regular estimate of theta
  res = meta_analysis(nthetas,Vthe) 
  coef1 = res$coefficients
  vcov1  = res$vcov
  
  list(coef = coef1,vcov = vcov1)
}


metagim_rho <- function(models,group,sample.info,form,family,data,cut0){
  ncase = sample.info$ncase
  nctrl = sample.info$nctrl
  
  nmodel = length(group)
  I0 <- rep(list(NULL), nmodel)
  J0 <- rep(list(NULL), nmodel)
  npara <- rep(list(NULL), nmodel)
  thetas <- rep(list(NULL), nmodel)
  V.opt <- rep(list(NULL), nmodel)

  
  n.para <- NULL
  score = list()
  nmap <- list()
  outcome <- list()
  
  V.the <- list()
  for (i in 1:nmodel){
    id = group[[i]]
    nid = length(id)
    md = list()
    for (j in 1:nid){
      idj = id[j]
      md[[j]] = models[[idj]]
    }
    nca = ncase[id,id,drop = FALSE]
    nct = nctrl[id,id,drop = FALSE]
    
    samplei.info = list(ncase = nca, nctrl = nct)
    
    fit = gimfit(form, family, data, md, samplei.info)
    I0[[i]] <- fit$Iv0
    J0[[i]] <- fit$Jv0
    thetas[[i]] <- fit$coefficients
    npara[[i]] <- fit$parms
    V.the[[i]] <- fit$vcov
    
    nmap[[i]] <- fit$map
    V.opt[[i]] <- fit$V.bet
    outcome[[i]] <- fit$outcome
    
    n.para[i] = length(fit$parms)
    score[[i]] = fit$score
  }
  
  ### compute new Sig0 for all summary data
  res0 <- organize(thetas,npara,nmap,group,sample.info)
  para = res0$para
  map = res0$map
  K = res0$K
  
  ## pr0,Delta based on the internal-based mle
  mf0 <- model.frame(form, data = data)
  outcome0 <- colnames(mf0)[1]
  n1 <- sum(data[, outcome0])
  n0 <- sum(1 - data[, outcome0])	
  fit0 <- glm(form, family = 'binomial', data = data)	  
  pr0 <- (1-fit0$fitted.values)/n0
  # sum(pr0) == 1
  Delta <- exp(fit0$linear.predictors) * n0/n1
  
  fp <- formula.parse(form, family, data, models, ref = NULL)
  ref <- fp$ref
  Sig0 = Sigma0(para, map, ref, K, sample.info, pr0, Delta)
  
  
  
  ### Update the Covariance matrix for all parameters
  lam = para[1]
  daty = data[, outcome0]
  Inew = Cov_Omega(npara,nmap,score,daty,lam,n.para,J0,Sig0)
  
  
  Vthe = Varcof_meta(nmap,J0,Inew)
  nthetas = theta_refine(thetas)
  
  
  ### One-step meta-analysis for theta
  res = meta_analysis_rho(nthetas,Vthe,cut0)
  
  coef1 = res$coefficients
  vcov1  = res$vcov
  
  list(coef = coef1,vcov = vcov1)
}





