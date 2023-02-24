\name{metagim_rho}
\alias{metagim_rho}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  metagim_rho
}
\description{Integrating High-Dimensional Summary Statistics}
\usage{
metagim_rho(models,group,sample.info,form,family,data,cut0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{models}{A list of external model and summary statistics.}
  \item{group}{A list of partitioned sample id numbers.}
  \item{sample.info}{A list of two matrices specifying the number of cases/controls shared in 
        datasets that are used to fit the working models given in \code{models}.}
  \item{form}{A formula.}
  \item{family}{Model family.}
  \item{data}{A data frame containing all variables that are specified in \code{formula} and \code{models}. }
  \item{cut0}{A value used for an ill-conditioned variance-covariance matrix of theta}
}

\examples{
  
  data(data, package="MetaGIM")

  # Random assignment into 10 batches
  set.seed(0)
  M      <- ncol(data)-2
  bk     <- 10
  M0     <- sample(1:M)
  group1 <- split(M0, ceiling((1:M)/bk))

  # True batch assignment 
  group2 <- split(1:M, ceiling((1:M)/bk))

  metagim_rho(models, group1, nsample, y~score, "case-control", data, 100)
  metagim_rho(models, group2, nsample, y~score, "case-control", data, 100)

}
