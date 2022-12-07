library(MetaGIM)

data(data, package="MetaGIM")

# Random assigned 10 batches for 100 summary data
set.seed(0)
M      <- ncol(data)-2
bk     <- 10
M0     <- sample(1:M)
group1 <- split(M0, ceiling((1:M)/bk))

fit = metagim(models, group1, nsample, y~score, "case-control", data)
fit$coefficients
fit$vcov
