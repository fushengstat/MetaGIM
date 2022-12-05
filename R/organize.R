### return all paras/map in one single argument
reorganize <- function(thetas,npara,nmap,group,sample.info){
  nmodel <- length(npara)
  
  ncase <- sample.info$ncase
  nctrl <- sample.info$nctrl
  
  the = ls_mean(thetas) 
  n.the = length(the)
  
  n.lam <- NULL
  n.alp <- NULL
  n.bet <- NULL
  alp <- NULL
  bet <- NULL 

  lam0 <- NULL
  alam1 <- NULL
  blam1 <- NULL
  

  map <- list()
  map$alp <- list(0:0)
  map$bet <- list(0:0)
  
  cj = 0
  for (i in 1:nmodel){
    parai = npara[[i]]
    mapi = nmap[[i]]
    id = group[[i]]
    nca = ncase[id,id]
    nca = as.matrix(nca)
    nct = nctrl[id,id]
    nct = as.matrix(nct)
    
  
    id.lam = mapi$lam
    ## multipier to Delta-1
    lam0 = c(lam0,parai[1])
    
    ## select the part of g
    id.lam1 = id.lam[-1]
    lam1 = parai[id.lam1]
    id.alp = mapi$all.alp
    alp = c(alp,parai[id.alp])
    id.bet = mapi$all.bet
    bet = c(bet,parai[id.bet])
    
    offset = max(mapi$the)
    ci.alp = id.alp - offset
    alam1 <- c(alam1,lam1[ci.alp])
    ci.bet = id.bet - offset
    blam1 <- c(blam1,lam1[ci.bet])
    
    ni = length(mapi$alp)
    for (s in 1:ni){
      j = cj + s
      alp.id = mapi$alp[[s]]
      alp0 = parai[alp.id]
      bet.id = mapi$bet[[s]]
      bet0 = parai[bet.id]
      map$alp[[j+1]] = max(map$alp[[j]])+ 1:length(alp0)
      map$bet[[j+1]] = max(map$bet[[j]])+ 1:length(bet0)
    }
    cj = cj + ni
    
  }
  
  
  map$alp[[1]] <- NULL
  map$bet[[1]] <- NULL
  
  

  lam <- c(mean(lam0),alam1,blam1)
  nlam = length(lam)
  names(lam) <- paste0('lam', 1:nlam)
  
  para <- c(lam, the, alp, bet)
  
  nalp <- length(alp)

  map$lam <- 1:nlam
  map$the <- (nlam + 1) : (nlam + n.the)
  all.alp <- NULL
  all.bet <- NULL
  
  ## K = number of all external models
  K = length(map$alp)
  
  for(i in 1:K){
    map$alp[[i]] <- map$alp[[i]] + nlam + n.the
    all.alp <- c(all.alp, map$alp[[i]])
    map$bet[[i]] <- map$bet[[i]] + nlam + n.the + nalp
    all.bet <- c(all.bet, map$bet[[i]])
  }
  
  map$all.alp <- sort(unique(all.alp))
  map$all.bet <- sort(unique(all.bet))
  

  
  list(para = para, map = map, K = K)
  
}






