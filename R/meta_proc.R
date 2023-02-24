## meta-analysis procedure
meta_analysis <- function(thetas,Vtheta){
  K = length(thetas)
  d = length(thetas[[1]])
  nms = names(thetas[[1]])

  theta = NULL  
  for (i in 1:K){
    theta = c(theta,thetas[[i]])
  }
  
  theta = matrix(theta,K*d,1)
  
  a = matrix(1,K,1)
  E = diag(d)
  S = kronecker(a,E)
  
  inv.W = solve(Vtheta)
  A1 = t(S) %*% inv.W %*% S
  A2 = t(S) %*% inv.W
  
  weight = solve(A1) %*% A2
  weight = matrix(weight,nrow = d)
  
  tem = weight %*% theta
  
  the_new = as.vector(tem)
  names(the_new) = nms

  V1 = solve(A1)
  colnames(V1) = nms
  rownames(V1) = nms
  list(coefficients = the_new, vcov = V1)
}


meta_analysis_rho <- function(thetas,Vtheta,cut0){
  K = length(thetas)
  d = length(thetas[[1]])
  nms = names(thetas[[1]])
  
  theta = NULL  
  for (i in 1:K){
    theta = c(theta,thetas[[i]])
  }
  
  theta = matrix(theta,K*d,1)
  
  a = matrix(1,K,1)
  E = diag(d)
  S = kronecker(a,E)
  
  Vsvd = svd(Vtheta)
  svalue = Vsvd$d
  smax = max(svalue)
  smin = min(svalue)
  tem0 = (smax-cut0*smin)/(cut0-1)
  rho = max(0,tem0)
  
  nv = nrow(Vtheta)
  Vtheta1 = Vtheta + rho*diag(nv)
  
  inv.W = solve(Vtheta1)
  A1 = t(S) %*% inv.W %*% S
  A2 = t(S) %*% inv.W
  
  weight = solve(A1) %*% A2
  weight = matrix(weight,nrow = d)
  
  tem = weight %*% theta
  
  the_new = as.vector(tem)
  names(the_new) = nms
  
  V1 = solve(A1) %*% A2 %*% Vtheta %*% t(A2) %*% solve(A1) 
  colnames(V1) = nms
  rownames(V1) = nms
  list(coefficients = the_new, vcov = V1)
}


