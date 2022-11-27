scorek_ref <- function(para, map, data, ref, model, V, sample.info, outcome){
  #return(grad(obj.cc, para, map = map, data = data, ref = ref, inv.V = inv.V, bet0 = bet0, outcome = outcome))
  inv.V <- solve(V)
  bet0 = get_bet0(model)
  
  ncase <- sample.info$ncase
  nctrl <- sample.info$nctrl
  
  nmodel <- length(map$bet)
  nlam <- max(map$lam)
  
  lam <- para[map$lam]
  xi <- lam[-1]

  ome <- para[map$ome]
  the <- para[map$the]
  fx <- as.matrix(data[, names(the), drop = FALSE])
  y <- data[, outcome]
  n1 <- sum(y)
  n0 <- sum(!y)
  rho <- n1/n0
  
  tilt <- tilt.ccr(para, map, data, ref)
  Delta <- tilt$Delta
  delta <- tilt$delta
  Delta0 <- tilt$Delta0
  
  g <- gfunction.cc(para, map, ref, Delta, delta, ncase, nctrl)
  
  pr <- as.vector(1/(1+g %*% lam))
  
  np <- length(para)
  
  ### reference set
  nr  = nrow(g)
  scr <- matrix(0, nrow = nr, ncol = np)
  colnames(scr) <- names(para)
  ### internal set
  ni  = nrow(fx)
  sci <- matrix(0, nrow = ni, ncol = np)
  colnames(sci) <- names(para)

  tmp0 = -g * pr
  scr[,map$lam] <- as.matrix(tmp0)
  
  pred <- 1 - 1/(1 + rho * Delta0)
  sci[,map$ome] <- y - pred
  
  

  
  dlogL <- as.matrix(fx * (y - pred))
  g.the <- gfunction.the.cc(para, map, ref, Delta, delta, ncase, nctrl)
  g.the.xi <- gfunction.the.xi.cc(g.the, xi)
  tmp <- g.the.xi + lam[1] * Delta * ref[, names(the)]
  tmp1 =  - tmp * pr
  sci[,map$the] <- as.matrix(dlogL)
  scr[,map$the] <- as.matrix(tmp1)

  if(!is.null(map$all.alp)){
    g.alp <- gfunction.alp.cc(para, map, ref, Delta, delta, ncase, nctrl)
    g.alp.xi <- gfunction.alp.xi.cc(g.alp, xi)
    tmp2 = -g.alp.xi * pr
    scr[,map$all.alp] <- as.matrix(tmp2)
  }
  

  g.bet <- gfunction.bet.cc(para, map, ref, Delta, delta, ncase, nctrl)
  g.bet.xi <- gfunction.alp.xi.cc(g.bet, xi)
  tmp3 = -g.bet.xi * pr
  scr[,map$all.bet] <- as.matrix(tmp3)
 
  
  
  list(scoreref = scr, scoreint = sci)
  
}

