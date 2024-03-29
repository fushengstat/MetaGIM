# BW 2021-07-23 mat is now a list which also includes Iv0 and Jv0

mcov <- function(para, map, family, data, ref, model, sample.info, V, bet0, outcome, type){
  
  if(family == 'gaussian'){
    mat <- cov.lm(para, map, data, ref, model, sample.info, V, bet0, outcome)
  }
  
  if(family == 'binomial'){
    mat <- cov.lo(para, map, data, ref, model, sample.info, V, bet0, outcome)
  }
  
  if(family == 'case-control'){
    if(type == 'cc-ref'){
      mat <- cov.ccr(para, map, data, ref, sample.info, V, bet0, outcome)
    }else{
      mat <- cov.cc(para, map, data, ref, sample.info, V, bet0, outcome)
    }
  }
  
  if(family == 'cml'){
    mat <- cov.cml(para, map, data, ref, sample.info, V, bet0, outcome)
  }
  
  mat
  
}

