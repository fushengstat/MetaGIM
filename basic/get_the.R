### get the index of all theta
theid <- function(nmap){
  K = length(nmap)
  the_id <- NULL
  id = 0
  for (i in 1:K){
    n.para = max(nmap[[i]]$all.bet)
    tem = id + nmap[[i]]$the
    the_id = c(the_id,tem)
    id  = id + n.para
  }
  
  the_id
}

### get the index of all theta without the intercept

theid_new <- function(nmap){
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

