withinrange_vec <- function(x,mea_vec,se_vec){
  n = length(x)
  z = rep(0,n)
  for (i in 1:n){
    xi = x[i]
    mi = mea_vec[i]
    sei = se_vec[i]
    ai = mi - 1.96*sei
    bi = mi + 1.96*sei
    
    z[i] = ifelse((xi <= bi & xi >= ai),1,0)
  }
 
  z
}


withinrange_mat <- function(X,Mean,Se){
  n = nrow(X)
  d = ncol(X)  
  Z = matrix(0,n,d)
  for (i in 1:n){
    for (j in 1:d){
      xi = X[i,j]
      mi = Mean[i,j]
      sei = Se[i,j]
      ai = mi - 1.96*sei
      bi = mi + 1.96*sei
      
      Z[i,j] = ifelse((xi <= bi & xi >= ai),1,0)      
      
    }

  }
 
  Z
}