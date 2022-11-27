### data generation for simulations
Xgen <- function(n,p,Sigma){
  M = length(p)
  bk = nrow(Sigma)
  nbk = M/bk

  Sig = kronecker(diag(1, nbk), Sigma)
  Xtem = rmvnorm(n, mean = rep(0,M), sigma=Sig)
  X = matrix(1,n,M)
  
  for(i in 1:M){
    pi = p[i]
    ## 0,1,2
    # prob = c((1-pi)^2,2*pi*(1-pi),pi^2)
    pb0 = c((1-pi)^2,1-pi^2)
    cut = qnorm(pb0)
    id1 = which(Xtem[,i]<cut[1])
    id2 = which(Xtem[,i]>cut[2])
    X[id1,i] = 0
    X[id2,i] = 2
  }
  
  colnames(X) <- paste0('X', 1:M)
  
  X
}

Zgen <- function(X,weight){
  w = matrix(weight,ncol = 1)
  Z = X %*% w
  
  Z
}


score.update <- function(data,weight,group){
  data$y <- NULL
  data$score <- NULL
  
  n = nrow(data)
  nscore = length(group)
  w = matrix(weight,ncol = 1)
  score = matrix(0,n,nscore)
  X = as.matrix(data)
  for(i in 1:nscore){
    id = group[[i]]
    wei = w[id]
    Xi = X[,id]
    score[,i] = Zgen(Xi,wei)
  }
  score
}


ygen <- function(Z,a,b){
  n = length(Z)
  pr = plogis(a + b*Z)
  y = rbinom(n,1,pr)
  y
}


int_simu <- function(n,p,w,Sigma,a,b){
  nhalf = n/2
  X0 =  Xgen(n,p,Sigma)
  Z0 = Zgen(X0,w)
  y0 = ygen(Z0,a,b)
  df0 = cbind.data.frame(y = y0, score = Z0, X0)
  
  n_test = table(df0$y)
  indicator = sum(n_test<=nhalf)
  
  while(indicator<2){
    X1 =  Xgen(n,p,Sigma)
    Z1 = Zgen(X1,w)
    y1 = ygen(Z1,a,b)
    df1 = cbind.data.frame(y = y1, score = Z1, X1)
    
    df0 <-rbind(df0, df1)
    n_test = table(df0$y)
    indicator = sum(n_test>=nhalf)
    # cat(n_test,"\n")
  }
  
  id.int = NULL
  case_control = c(0,1)
  for(j in case_control){
    idj = which(df0$y==j)
    id.int = c(id.int,idj[1:nhalf])
  }
  int = df0[id.int,]
  
  int
}

ext_simu <- function(N,p,w,Sigma,a,b){
  nhalf = N/2
  X0 =  Xgen(N,p,Sigma)
  Z0 = Zgen(X0,w)
  y0 = ygen(Z0,a,b)
  df0 = cbind.data.frame(y = y0, score = Z0, X0)
  
  n_test = table(df0$y)
  indicator = sum(n_test<=nhalf)
  
  while(indicator<2){
    X1 = Xgen(N,p,Sigma)
    Z1 = Zgen(X1,w)
    y1 = ygen(Z1,a,b)
    df1 = cbind.data.frame(y = y1, score = Z1, X1)
    
    df0 <-rbind(df0, df1)
    n_test = table(df0$y)
    indicator = sum(n_test>=nhalf)
    # cat(n_test,"\n")
  }
  
  id.ext = NULL
  case_control = c(0,1)
  for(j in case_control){
    idj = which(df0$y==j)
    id.ext = c(id.ext,idj[1:nhalf])
  }
  ext = df0[id.ext,]
  
  ext
}


ref_simu <- function(N,p,w,Sigma,a,b){
  X0 =  Xgen(N,p,Sigma)
  Z0 = Zgen(X0,w)
  y0 = ygen(Z0,a,b)
  df0 = cbind.data.frame(y = y0, score = Z0, X0)
  

  indicator = 0
  
  while(indicator<1){
    X1 = Xgen(N,p,Sigma)
    Z1 = Zgen(X1,w)
    y1 = ygen(Z1,a,b)
    df1 = cbind.data.frame(y = y1, score = Z1, X1)
    
    df0 <-rbind(df0, df1)
    t1 = sum(df0$y==0)
    indicator = (t1>N)
    
    cat(t1,"\n")
  }
  
  
  idj = which(df0$y==0)
  id.ref = idj[1:N]
  ref = df0[id.ref,]
  
  ref
}


ref_simu1 <- function(N,p,w,Sigma,a,b,nper){
  ref0 = ref_simu(N,p,w,Sigma,a,b)
  ref_0 = ref0
  ref0$y = NULL
  ref0$score = NULL

  ref0 = as.matrix(ref0)
  M = length(w)
  
  Xr = matrix(0,N*nper,M)
  bk = nrow(Sigma)
  nbk = M/bk
  
  for(i in 1:nbk){
    id.b1 = (i-1)*bk+1
    id.b2 = i*bk
    set.seed(i)
    idi = sample(1:N,N*nper,replace=TRUE)
    Xr[,id.b1:id.b2] = ref0[idi,id.b1:id.b2]
    # cat(i,id.b1,id.b2,"\n")
    # cat(idi,"\n")
  }
  
  colnames(Xr) <- paste0('X', 1:M)
  
  yr = rep(0,N*nper)
  w = matrix(w,ncol = 1)
  Zr = Xr %*% w
  ref = cbind.data.frame(y = yr, score = Zr, Xr)

  list(ref = ref, ref0 = ref_0)
  
}


ref_simu_p <- function(N,p,w,Sigma,a,b,nper){
  ref0 = ref_simu(N,p,w,Sigma,a,b)
  ref0$y = NULL
  ref0$score = NULL
  
  ref0 = as.matrix(ref0)
  M = ncol(ref0)
  
  Xr0 = matrix(0,N*(nper-1),M)
  bk = nrow(Sigma)
  nbk = M/bk
  
  for(i in 1:nbk){
    id.b1 = (i-1)*bk+1
    id.b2 = i*bk
    
    va = rep(1:N,each = nper-1)
    set.seed(i)
    idx = sample(va)
    Xr0[,id.b1:id.b2] = ref0[idx,id.b1:id.b2]
  }
  
  Xr = rbind(ref0,Xr0)
  row.names(Xr) <- NULL
  
  yr = rep(0,N*nper)
  w = matrix(w,ncol = 1)
  Zr = Xr %*% w
  ref = cbind.data.frame(y = yr, score = Zr, Xr)
  
  ref
  
}




int_simu0<- function(n,p,w,Sigma,a,b){
  X0 =  Xgen(n,p,Sigma)
  Z0 = Zgen(X0,w)
  y0 = ygen(Z0,a,b)
  int = cbind.data.frame(y = y0, score = Z0, X0)
  
  int
}



