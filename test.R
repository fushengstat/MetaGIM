rm(list = ls())
sf0 = "C:/Users/shengfu/Box/Metagim"
setwd(sf0)

SRCDIR <- "C:/Users/shengfu/Box/Metagim/basic"
for (SF in list.files(SRCDIR, full.names=TRUE)) source(SF) 


## we need the latest gim with version: gim_0.34.9  
library(gim)

## 4 varaibles, data--internal data; 
## score.add--score matrix based on the internal data and independent batch assignment
## models--the considered external model and summary statistics
## nsample-a list of two matrices specifying the number of cases/controls shared in 
## datasets that are used to fit the working models given in model.
load("test.Rdata")
 


form = "y ~ score"
family="case-control"

#batch assignment
## Random assignment
M = ncol(data)-2
bk = 10
set.seed(0)
M0 = sample(1:M)
group1 = split(M0, ceiling((1:M)/bk))

## true batch assignment as the data generating
group2 = split(1:M, ceiling((1:M)/bk))


## independent assignment 
M0 = (1:M) %% bk
batch0 = list()
ng = M/bk
for(i in 1:ng){
  batch0[[i]] = which(M0==(i%%bk))
}
batches = unlist(batch0)
group3 = split(batches,ceiling((1:M)/bk))



fit0 = metagim(models,group1,nsample,form,family,data)

## a is used for ill-conditioned variance-covariance matrix of theta
a = 100
fit1 = metagim_rho(models,group1,nsample,form,family,data,a)
fit2 = metagim_rho(models,group2,nsample,form,family,data,a)
fit3 = metagim_rho(models,group3,nsample,form,family,data,a)


# score.add is the independent scores based on the true batch assignment
fit4 = gim_add(models,group2,nsample,form,family,data,score.add)
