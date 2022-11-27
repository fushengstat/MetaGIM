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


