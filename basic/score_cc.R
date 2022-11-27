scorek <- function(para, map, ref, model, V, sample.info, outcome){

  #return(grad(obj.cc, para, map = map, data = data, ref = ref, inv.V = inv.V, bet0 = bet0, outcome = outcome))
  inv.V <- solve(V)
  bet0 = get_bet0(model)
  
  ncase <- sample.info$ncase
  nctrl <- sample.info$nctrl
  
  nmodel <- length(map$bet)
  nlam <- max(map$lam)
  
  lam <- para[map$lam]
  xi <- lam[-1]
  
  the <- para[map$the]
  fx <- as.matrix(ref[, names(the), drop = FALSE])
  y <- ref[, outcome]
  
  tilt <- tilt.cc(para, map, ref)
  Delta <- tilt$Delta
  delta <- tilt$delta
  
  g <- gfunction.cc(para, map, ref, Delta, delta, ncase, nctrl)
  
  pr <- as.vector(1/(1+g %*% lam))

  np <- length(para)
  
  n  = nrow(g)
  sc <- matrix(0, nrow = n, ncol = np)
  colnames(sc) <- names(para)
  tmp0 = -g * pr
  sc[,map$lam] <- as.matrix(tmp0)

  dlogL <- as.matrix(fx * y)
  g.the <- gfunction.the.cc(para, map, ref, Delta, delta, ncase, nctrl)
  g.the.xi <- gfunction.the.xi.cc(g.the, xi)
  tmp <- g.the.xi + lam[1] * Delta * ref[, names(the)]
  tmp1 = dlogL - tmp * pr
  sc[,map$the] <- as.matrix(tmp1)

  
  if(!is.null(map$all.alp)){
    g.alp <- gfunction.alp.cc(para, map, ref, Delta, delta, ncase, nctrl)
    g.alp.xi <- gfunction.alp.xi.cc(g.alp, xi)
    tmp2 = -g.alp.xi * pr
    sc[,map$all.alp] <- as.matrix(tmp2)
  }
  
  # bet <- para[map$all.bet]
  # tmp = 1/n*t(bet - bet0) %*% inv.V
  # dqf <- matrix(rep(tmp,n),nrow = n,byrow = TRUE)
  g.bet <- gfunction.bet.cc(para, map, ref, Delta, delta, ncase, nctrl)
  g.bet.xi <- gfunction.alp.xi.cc(g.bet, xi)
  # tmp3 = -g.bet.xi * pr - dqf
  tmp3 = -g.bet.xi * pr
  sc[,map$all.bet] <- as.matrix(tmp3)
  
  sc
  
}

