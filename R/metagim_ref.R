# main procedure for gim 
metagim_ref <- function(models,group,sample.info,form,family,data,ref){
  ncase = sample.info$ncase
  nctrl = sample.info$nctrl

  nmodel = length(group)
  I0 <- rep(list(NULL), nmodel)
  J0 <- rep(list(NULL), nmodel)
  npara <- rep(list(NULL), nmodel)
  thetas <- rep(list(NULL), nmodel)
  V.opt <- rep(list(NULL), nmodel)
  pr0 <- rep(list(NULL), nmodel)
  Delta <- rep(list(NULL), nmodel)
  score_ref <- rep(list(NULL), nmodel)
  score_int <- rep(list(NULL), nmodel)

  omes <- NULL
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
    
    fit = gimfit_ref(form, family, data, md, ref,  samplei.info)
    
    I0[[i]] <- fit$Iv0
    J0[[i]] <- fit$Jv0
    omes[i] <- fit$ome
    thetas[[i]] <- fit$coefficients
    npara[[i]] <- fit$parms
    V.the[[i]] <- fit$vcov

    nmap[[i]] <- fit$map
    V.opt[[i]] <- fit$V.bet
    outcome[[i]] <- fit$outcome

    n.para[i] = length(fit$parms)
    pr0[[i]] = fit$pr0
    Delta[[i]] = fit$Delta
    score_int[[i]] = fit$scorei
    score_ref[[i]] = fit$scorer
  }
  
  ### compute new Sig0 for all summary data
  res0 <- reorganize_ref(thetas,omes,npara,nmap,group,sample.info)
  para = res0$para
  map = res0$map
  K = res0$K
  
  ## pr0,Delta based on the internal-based mle
  pr0 <- ls_mean(pr0)
  # sum(pr0) == 1
  Delta <- ls_mean(Delta)
  
  ### compute Sig0 for all external summary data
  ref <- fit$ref
  Sig0 = Sigma0.ccr(para, map, ref, K, sample.info, pr0, Delta)
  
  
  ### Update the Covariance matrix for all parameters
  outcome0 = outcome[[1]]
  daty = data[, outcome0]
  Inew = Cov_Omega_ref(npara,nmap,score_ref,score_int,daty,n.para,J0,Sig0)

  
  Vthe = Varcof_meta(nmap,J0,Inew)
  nthetas = theta_refine(thetas)
  ### One-step meta-analysis for theta
  
  # regular estimate of theta
  res = meta_analysis(nthetas,Vthe) 
  
  coef1 = res$coefficients
  vcov1  = res$vcov
  
  list(coef = coef1,vcov = vcov1)  
}



metagim_ref_rho <- function(models,group,sample.info,form,family,data,ref,cut0){
  ncase = sample.info$ncase
  nctrl = sample.info$nctrl
  
  nmodel = length(group)
  I0 <- rep(list(NULL), nmodel)
  J0 <- rep(list(NULL), nmodel)
  npara <- rep(list(NULL), nmodel)
  thetas <- rep(list(NULL), nmodel)
  V.opt <- rep(list(NULL), nmodel)
  pr0 <- rep(list(NULL), nmodel)
  Delta <- rep(list(NULL), nmodel)
  score_ref <- rep(list(NULL), nmodel)
  score_int <- rep(list(NULL), nmodel)
  
  omes <- NULL
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
    
    fit = gimfit_ref(form, family, data, md, ref,  samplei.info)
    
    I0[[i]] <- fit$Iv0
    J0[[i]] <- fit$Jv0
    omes[i] <- fit$ome
    thetas[[i]] <- fit$coefficients
    npara[[i]] <- fit$parms
    V.the[[i]] <- fit$vcov
    
    nmap[[i]] <- fit$map
    V.opt[[i]] <- fit$V.bet
    outcome[[i]] <- fit$outcome
    
    n.para[i] = length(fit$parms)
    pr0[[i]] = fit$pr0
    Delta[[i]] = fit$Delta
    score_int[[i]] = fit$scorei
    score_ref[[i]] = fit$scorer
  }
  
  ### compute new Sig0 for all summary data
  res0 <- reorganize_ref(thetas,omes,npara,nmap,group,sample.info)
  para = res0$para
  map = res0$map
  K = res0$K
  
  ## pr0,Delta based on the internal-based mle
  pr0 <- ls_mean(pr0)
  # sum(pr0) == 1
  Delta <- ls_mean(Delta)
  
  ### compute Sig0 for all external summary data
  ref <- fit$ref
  Sig0 = Sigma0.ccr(para, map, ref, K, sample.info, pr0, Delta)
  
  
  ### Update the Covariance matrix for all parameters
  outcome0 = outcome[[1]]
  daty = data[, outcome0]
  Inew = Cov_Omega_ref(npara,nmap,score_ref,score_int,daty,n.para,J0,Sig0)
  
  
  Vthe = Varcof_meta(nmap,J0,Inew)
  nthetas = theta_refine(thetas)
  ### One-step meta-analysis for theta
  
  # regular estimate of theta
  res = meta_analysis_rho(nthetas,Vthe,cut0)
  
  coef1 = res$coefficients
  vcov1  = res$vcov
  
  list(coef = coef1,vcov = vcov1)
}

