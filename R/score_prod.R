score_prod <- function(scorei,scorej,daty){
  n0 = sum(1-daty)
  n1 = sum(daty)
  n = n0 + n1
  
  idy1 = which(daty==1)
  idy0 = which(daty==0)
  size1 = length(idy1)
  size0 = length(idy0)
  

  scorei1 = matrix(colSums(scorei[idy1,]),nrow = 1)
  scorej1 = matrix(colSums(scorej[idy1,]),nrow = 1)

  scorei0 = matrix(colSums(scorei[idy0,]),nrow = 1)
  scorej0 = matrix(colSums(scorej[idy0,]),nrow = 1)
  
  
  tem = 1/size1* t(scorei1) %*% scorej1 + 1/size0* t(scorei0) %*% scorej0
  
  tem1 = t(scorei) %*% scorej
  
  list(tem = tem, tem1 = tem1)
}

