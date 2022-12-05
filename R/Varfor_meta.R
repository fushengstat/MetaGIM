Varcof_meta <- function(nmap,J0,Inew){
  K = length(nmap)
  In = Inew$I
  Id = Inew$Id
  Inv.J = Jinv(J0)
  
  Vall = In*0
  for (i in 1:K){
    i1 = Id[[i]]
    for(j in i:K){
      j1 = Id[[j]]
      Vall[i1,j1] = Inv.J[[i]] %*% In[i1,j1] %*% Inv.J[[j]]
      
      if(i == j){
        next
      }
      Vall[j1,i1] = t(Vall[i1,j1])
    }
    
  }
    
  id.the = theid(nmap)
  Vthe = Vall[id.the,id.the]
  
  Vthe = as.matrix(Vthe)
  
  nms = colnames(In[,id.the])
  
  rownames(Vthe) = nms
  colnames(Vthe) = nms
  
  Vthe
}


