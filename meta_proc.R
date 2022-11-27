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


meta_analysis_new <- function(thetas,Vtheta){
  K = length(thetas)
  d = length(thetas[[1]])
  nms = names(thetas[[1]])
  
  nmcov = colnames(Vtheta)
  
  theta = NULL
  for (i in 1:K){
    theta = c(theta,thetas[[i]])
  }
  
  theta = matrix(theta,K*d,1)
  
  a = matrix(1,K,1)
  
  coef = NULL
  Vcoef = NULL
  
  cid = (0:(K-1))*d
  
  for (j in 1:d){
    idj = cid+j

    Theta = theta[idj]
    Theta_cov = Vtheta[idj,idj]
    
    inv.W = solve(Theta_cov)
    a1 = t(a) %*% inv.W %*% a
    A2 = t(a) %*% inv.W 
    
    coef[j] = solve(a1) %*% A2 %*% Theta
    Vcoef[j] = solve(a1)
  }
  
  coef = as.vector(coef)
  names(coef) = nms
  
  V1 = matrix(0,d,d)
  colnames(V1) = nms
  rownames(V1) = nms
  
  diag(V1) = Vcoef
  
  list(coefficients = coef, vcov = V1) 
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



meta_ind <- function(thetas,V.the){
  K = length(thetas)
  nthe = rep(0,K)
  Vars = rep(0,K)
  
  for(i in 1:K){
    nthe[i] = thetas[[i]]
    Vars[i] = V.the[[i]][-1,-1]
  }
  
  vcov = 1/sum(1/Vars)
  meta.w = 1/Vars*vcov
  theta = as.vector(meta.w%*%nthe)
  
  list(coefficients = theta, vcov = vcov)
}




