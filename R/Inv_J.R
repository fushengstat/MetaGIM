## return the inverse of J matrix
Jinv <- function(J){
  K = length(J)
  
  Jinv = list()
  
  for (i in 1:K){
    Jinv[[i]] = solve(J[[i]])
  }
  
  Jinv
}