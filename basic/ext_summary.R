ext_summary <- function(data_ext){
  ext = data_ext
  y = ext$y
  
  K = ncol(ext)-1
  nms = colnames(ext)[-1]
  
  models = list()
  
  for(i in 1:K){
    formi = paste0("y", " ~ ", nms[i])
    mi <- glm(formi, data = ext, family = "binomial")
    
    models[[i]] = list(formula = formi, 
                       info = data.frame(var = names(coef(mi))[-1], 
                                         bet = coef(mi)[-1], stringsAsFactors = FALSE))
  }
  
  ncase <- matrix(sum(y), K, K)
  nctrl <- matrix(sum(!y), K, K)
  
  sample.info = list(ncase = ncase, nctrl = nctrl)
  
  list(models = models, sample.info = sample.info)
}