### return the mean of a series of list --vectors
ls_mean <- function(lst){
  K = length(lst)
  d = length(lst[[1]])
  nms = names(lst[[1]])
  
  t_lst = matrix(0,nrow = K, ncol = d)
  
  for(i in 1:K){
    t_lst[i,] = lst[[i]]
  }
  
  nls = colMeans(t_lst)
  names(nls) <- nms
  
  nls
}


### return the standard error of a series of list --vectors
ls_sd <- function(lst){
  K = length(lst)
  d = length(lst[[1]])
  nms = names(lst[[1]])
  
  t_lst = matrix(0,nrow = K, ncol = d)
  
  for(i in 1:K){
    t_lst[i,] = lst[[i]]
  }
  
 
  nsd = apply(t_lst,2,sd)
  names(nsd) <- nms
  
  
  nsd
}


### return the mean and standard error of a series of list --vectors
stat_mean <- function(lst){
  K = length(lst)
  d = length(lst[[1]])
  nms = names(lst[[1]])
  
  t_lst = matrix(0,nrow = K, ncol = d)
  
  for(i in 1:K){
    t_lst[i,] = lst[[i]]
  }
  
  nls = colMeans(t_lst)
  names(nls) <- nms
  
  nsd = apply(t_lst,2,sd)
  names(nsd) <- nms
  
  
  list(mean = nls, sd = nsd)
}



### return the mean of a series of list
mat_mean <- function(lst){
  K = length(lst)
  d = dim(lst[[1]])
  rnms = rownames(lst[[1]])
  cnms = colnames(lst[[1]])
  
  t_lst = matrix(0,nrow = d[1], ncol = d[2])
  
  for(i in 1:K){
    t_lst = t_lst + lst[[i]]
  }
  
  nls = t_lst/K
  
  rownames(nls) <- rnms
  colnames(nls) <- cnms
  
  nls
}

### return the mean of a series of list
mat_sd <- function(lst){
  K = length(lst)
  tem = lst[[1]]
  rnms = rownames(tem)
  cnms = colnames(tem)
  d1 = nrow(tem)
  d2 = ncol(tem) 
  
  t_lst = matrix(0,nrow = K, ncol = d1*d2)
  
  for(i in 1:K){
    t_lst[i,] = c(lst[[i]])
  }
  
  sd0 = apply(t_lst,2,sd)
  
  nsd = matrix(sd0,d1,d2)  
  rownames(nsd) <- rnms
  colnames(nsd) <- cnms
  
  nsd
}



### return the mean of a list of matrix
mattovec <- function(lst){
  K = length(lst)
  d = dim(lst[[1]])
  rnms = rownames(lst[[1]])
  
  t_lst = rep(0,d[1]*d[2])
  lst_vec = matrix(0,K,d[1]*d[2])
  
  for(i in 1:K){
    t_lst = t_lst + c(lst[[i]])
    lst_vec[i,] = c(lst[[i]])
  }
  
  nls = t_lst/K
  
  names(nls) <- rep(rnms,d[2])
  colnames(lst_vec) <- rep(rnms,d[2])
 
  list(mlst = nls, vlst = lst_vec)
}


### list to vector
lsttovec <- function(lst){
  K = length(lst)
  vec = NULL
  for(i in 1:K){
    vec = c(vec,lst[[i]])
  }
  
  vec
}


### list to matrix
lsttomat <- function(lst){
  K = length(lst)
  d = length(lst[[1]])
  
  nms = names(lst[[1]])
	
  mat = matrix(0,K,d)
  for(i in 1:K){
    mat[i,] = lst[[i]]
  }
  
  colnames(mat) <- nms
  
  mat
}



### take the upper triangle part of a matrix
uppertri <- function(A){
  d = ncol(A)
  
  a = NULL
  m = NULL
  
  for(i in 1:d){
    for (j in i:d){
      m = c(m,paste0("Z",i,j))
      a = c(a,A[j,i])
    }
  }
  
  names(a) = m
  
  a
}


### return the estimate of covariance matrix
cov_lst <- function(lst){
  K = length(lst)
  d = length(lst[[1]])
  nms = names(lst[[1]])
  
  mat = matrix(0,K,d)
  
  for(i in 1:K){
    mat[i,] = c(lst[[i]])
  }
  
  colnames(mat) <- nms
  cov = cov(mat)
  list(cov = cov,mat = mat)
}


getsd <- function(V.the){
  K = length(V.the)
  d = ncol(V.the[[1]])-1
  
  sd.matrix = matrix(0,K,d)
  
  for(i in 1:K){
    tem = diag(V.the[[i]])[-1]
    sd.matrix[i,] = sqrt(tem)
  }
  
  sd.matrix
}

