### get the index of all theta without the intercept
theid <- function(nmap){
  K = length(nmap)
  the_id <- NULL
  id = 0
  for (i in 1:K){
    n.para = max(nmap[[i]]$all.bet)
    tem = id + nmap[[i]]$the[-1]
    the_id = c(the_id,tem)
    id  = id + n.para
  }
  
  the_id
}



### new theta without the intercepts
theta_refine <- function(thetas){
  K = length(thetas)
  nthetas = rep(list(NULL),K)
  for(i in 1:K){
    nthetas[[i]] = thetas[[i]][-1]
  }
  
  nthetas  
}


## get the values of all bet0
get_bet0 <- function(model){
  nmodel <- length(model)
  
  bet0 <- NULL
  for(i in 1:nmodel){
    mdi = model[[i]]
    tem = mdi$info[2]
    tem = as.numeric(tem)
    names(tem) <- mdi$info[1]
    bet0 <- c(bet0, tem)
}
  bet0
}



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

