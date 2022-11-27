gimfit <- function(form, family, data, model, samplei.info){
  n = nrow(data)
  
 
  ncase = samplei.info$ncase
  nctrl = samplei.info$nctrl
  
  fit <- gim(form, family, data, model, ncase = ncase, nctrl = nctrl)
  outcomei = fit$outcome
  parai = fit$parms
  mapi = fit$map
  Vi = fit$V.bet
  refi = fit$ref

  d1 = length(parai)
  score = scorek(parai, mapi, refi, model, Vi, samplei.info, outcomei)
  
  fit$score = score
  
  fit
  
}