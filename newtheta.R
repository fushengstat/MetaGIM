### new theta without the intercepts

newtheta <- function(thetas){
  K = length(thetas)
  nthetas = rep(list(NULL),K)
  for(i in 1:K){
    nthetas[[i]] = thetas[[i]][-1]
  }
  
  nthetas
  
}