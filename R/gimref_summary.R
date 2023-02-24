## summary of gim_ref fitting for each batch

gimfit_ref <- function(form, family, data, model, ref, samplei.info){
  n = nrow(ref)
  
  ncase = samplei.info$ncase
  nctrl = samplei.info$nctrl
  
  fit <- gim(form, family, data, model, ncase = ncase, nctrl = nctrl, ref = ref)
  outcomei = fit$outcome
  parai = fit$parms
  mapi = fit$map
  Vi = fit$V.bet
  refi = fit$ref
  fit$ome = parai[mapi$ome]
  
  
  lam <- parai[mapi$lam]
  # add a column of intercept
  data$'(Intercept)' <- 1
  
  tilt <- tilt.ccr(parai, mapi, data, refi)
  Delta <- tilt$Delta
  delta <- tilt$delta
  g <- gfunction.cc(parai, mapi, refi, Delta, delta, ncase, nctrl)
  pr <- as.vector(1/(1+g %*% lam))
  
  fit$pr0 = pr/n
  # fit$pr0 = 1/n
  fit$Delta = Delta
  
  ## sample-based score using the internal and reference dataset
  score = scorek_ref(parai, mapi, data, refi, model, Vi, samplei.info, outcomei)
  fit$scorei = score$scoreint
  fit$scorer = score$scoreref
  
  
  fit
}



