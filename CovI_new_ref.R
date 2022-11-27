## return the covariance of the separate gim procedures
nCov_ref <- function(npara,nmap,scoreref,scoreint,daty,n.para,J0,Sig0){
  K = length(n.para)
  d.para = sum(n.para)
  Iid = cumsum(c(1,n.para))[-(K+1)]
  Jid = cumsum(n.para)
  
  nr = nrow(scoreref[[1]])
  ni = nrow(scoreint[[1]])
  
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
    

    iscorei = scoreint[[i]]
    rscorei = scoreref[[i]]	
    
    for (j in 1:K){
      j1 = Iid[j]
      j2 = Jid[j]

      jlamid = nmap[[j]]$lam 
      jtheid = nmap[[j]]$the 
      jbetid = nmap[[j]]$all.bet 
      
      tlam = jlamid + j1-1
      tthe = jtheid + j1-1
      tbet = jbetid + j1-1
      
      
      iscorej = scoreint[[j]]
      rscorej = scoreref[[j]]	
	
      I_new[i1:i2,j1:j2]  = t(iscorei) %*% iscorej + t(rscorei) %*% rscorej 
      
      nres = score_prodref(iscorei,iscorej,rscorei,rscorej,daty)
      temr = nres$temr
      temi = nres$temi
      
      ns1 = names(npara[[i]][ibetid])
      ns2 = names(npara[[j]][jbetid])
      nSig = Sig0[ns1,ns2,drop = FALSE]
      
      # S1 = J0[[i]][ibetid,ibetid,drop = FALSE]
      # S2 = J0[[j]][jbetid,jbetid,drop = FALSE]
      # I_new[sbet,tbet] = S1 %*% nSig %*% S2
      
      # S1 = V.opt[[i]]
      # S2 = V.opt[[j]]
      
      S1 = Sig0[ns1,ns1,drop = FALSE]
      S2 = Sig0[ns2,ns2,drop = FALSE]
      I_new[sbet,tbet] = I_new[sbet,tbet] + solve(S1) %*% nSig %*% solve(S2)
      
      

      I_new[i1:i2,j1:j2] = I_new[i1:i2,j1:j2] - temr - temi
      

    }
    
    nms = c(nms,names(npara[[i]]))
  }
  colnames(I_new) = nms
  rownames(I_new) = nms
  
  I = I_new
  
  list(I = I_new, Id = Id)
}


