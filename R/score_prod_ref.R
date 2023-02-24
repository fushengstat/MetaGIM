score_prodref <- function(iscorei,iscorej,rscorei,rscorej,daty){
  n0 = sum(1-daty)
  n1 = sum(daty)
  n = n0 + n1
  
  idy1 = which(daty==1)
  idy0 = which(daty==0)
 
  size1 = length(idy1)
  size0 = length(idy0) 

  scorei1 = matrix(colSums(iscorei[idy1,]),nrow = 1)
  scorej1 = matrix(colSums(iscorej[idy1,]),nrow = 1)
  scorei0 = matrix(colSums(iscorei[idy0,]),nrow = 1)
  scorej0 = matrix(colSums(iscorej[idy0,]),nrow = 1)  
  
  temi = 1/size1* t(scorei1) %*% scorej1 + 1/size0* t(scorei0) %*% scorej0
  
  
  ## for the reference set
  nr = nrow(rscorei)
  rscore1 = matrix(colSums(rscorei),nrow = 1)
  rscore2 = matrix(colSums(rscorej),nrow = 1)
  
  temr = 1/nr * t(rscore1) %*% rscore2
  
  list(temi = temi, temr = temr)
}
