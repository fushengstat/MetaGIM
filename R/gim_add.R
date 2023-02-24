# main procedure for gim-add procedure 
gim_add <- function(models,group,sample.info,form,family,data,scores){
  ncase = sample.info$ncase
  nctrl = sample.info$nctrl
  
  nmodel = length(group)
  thetas <- rep(list(NULL), nmodel)
  outcome <- list()
  V.the <- list()
  
  data$score = NULL
  
  for (i in 1:nmodel){
    id = group[[i]]
    nid = length(id)
    md = list()
    for (j in 1:nid){
      idj = id[j]
      md[[j]] = models[[idj]]
    }
    nca = ncase[id,id,drop = FALSE]
    nct = nctrl[id,id,drop = FALSE]
    
    samplei.info = list(ncase = nca, nctrl = nct)
    data$score = scores[,i]
    
    fit = gimfit(form, family, data, md, samplei.info)
    
    thetas[[i]] <- fit$coefficients
    V.the[[i]] <- fit$vcov
    outcome[[i]] <- fit$outcome
  }
  
  
  nthe = rep(0,nmodel)
  Vars = rep(0,nmodel)
  nms = names(thetas[[1]])
  
  for(i in 1:nmodel){
    nthe[i] = thetas[[i]][-1]
    Vars[i] = V.the[[i]][-1,-1]
  }
  
  vcov = 1/sum(1/Vars)
  meta.w = 1/Vars*vcov
  theta = as.vector(meta.w%*%nthe)
  
  vcov = as.matrix(vcov)
  names(theta) = nms[-1]
  colnames(vcov) = nms[-1] 
  rownames(vcov) = nms[-1] 
  
  list(coefficients = theta, vcov = vcov)
  
}


gimref_add <- function(models,group,sample.info,form,family,data,ref,scores,scorer){
  ncase = sample.info$ncase
  nctrl = sample.info$nctrl

  nmodel = length(group)
  thetas <- rep(list(NULL), nmodel)
  outcome <- list()
  V.the <- list()
  
  data$score = NULL
  ref$score = NULL
  
  for (i in 1:nmodel){
    id = group[[i]]
    nid = length(id)
    md = list()
    for (j in 1:nid){
      idj = id[j]
      md[[j]] = models[[idj]]
    }
    nca = ncase[id,id,drop = FALSE]
    nct = nctrl[id,id,drop = FALSE]
    
    samplei.info = list(ncase = nca, nctrl = nct)
    data$score = scores[,i]
    ref$score = scorer[,i]
    
    
    fit = gimfit_ref(form, family, data, md, ref, samplei.info)
    
  
    thetas[[i]] <- fit$coefficients
    V.the[[i]] <- fit$vcov
    outcome[[i]] <- fit$outcome
  }
  
 
  nthe = rep(0,nmodel)
  Vars = rep(0,nmodel)
  nms = names(thetas[[1]])
	
  for(i in 1:nmodel){
    nthe[i] = thetas[[i]][-1]
    Vars[i] = V.the[[i]][-1,-1]
  }

  vcov = 1/sum(1/Vars)
  meta.w = 1/Vars*vcov
  theta = as.vector(meta.w%*%nthe)
  
  vcov = as.matrix(vcov)
  names(theta) = nms[-1]
  colnames(vcov) = nms[-1] 
  rownames(vcov) = nms[-1] 
  
  list(coefficients = theta, vcov = vcov)
  
}





