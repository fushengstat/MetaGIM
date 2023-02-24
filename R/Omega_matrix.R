## return the covariance of the separate gim procedures
Cov_Omega <- function(npara,nmap,score,daty,lam,n.para,J0,Sig0){
  K = length(n.para)
  d.para = sum(n.para)
  Iid = cumsum(c(1,n.para))[-(K+1)]
  Jid = cumsum(n.para)
  n = nrow(score[[1]])
  
  I_new = matrix(0,d.para,d.para)
  I1 = matrix(0,d.para,d.para)
  Id <- list()
  nms <- NULL
  
  for (i in 1:K){
    Id[[i]] = c(Iid[i]:Jid[i])
    i1 = Iid[i]
    i2 = Jid[i]
    
    ilamid = nmap[[i]]$lam 
    itheid = nmap[[i]]$the 
    ibetid = nmap[[i]]$all.bet 
    
    slam = ilamid + i1-1
    sthe = itheid + i1-1
    sbet = ibetid + i1-1
    

    scorei = score[[i]]
    
    for (j in 1:K){
      j1 = Iid[j]
      j2 = Jid[j]

      jlamid = nmap[[j]]$lam 
      jtheid = nmap[[j]]$the 
      jbetid = nmap[[j]]$all.bet 
      
      tlam = jlamid + j1-1
      tthe = jtheid + j1-1
      tbet = jbetid + j1-1
      
      
      scorej = score[[j]]
      
      I_new[i1:i2,j1:j2]  = t(scorei) %*% scorej
      
      nres = score_prod(scorei,scorej,daty)
      tem = nres$tem
      
      ns1 = names(npara[[i]][ibetid])
      ns2 = names(npara[[j]][jbetid])
      nSig = Sig0[ns1,ns2,drop = FALSE]
      
      
      S1 = Sig0[ns1,ns1,drop = FALSE]
      S2 = Sig0[ns2,ns2,drop = FALSE]
      I_new[sbet,tbet] = I_new[sbet,tbet] + solve(S1) %*% nSig %*% solve(S2)
      
      

      I_new[i1:i2,j1:j2] = I_new[i1:i2,j1:j2] - tem
      

    }
    
    nms = c(nms,names(npara[[i]]))
  }
  colnames(I_new) = nms
  rownames(I_new) = nms
  
  I = I_new
  
  list(I = I_new, Id = Id)
}

