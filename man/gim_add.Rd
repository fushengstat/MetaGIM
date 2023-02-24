\name{gim_add}
\alias{gim_add}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  gim_add
}
\description{Integrating High-Dimensional Summary Statistics}
\usage{
  gim_add(models,group,sample.info,form,family,data,scores)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{models}{A list of external model and summary statistics.}
  \item{group}{A list of partitioned sample id numbers.}
  \item{sample.info}{A list of two matrices specifying the number of cases/controls shared in 
        datasets that are used to fit the working models given in \code{models}.}
  \item{form}{A formula.}
  \item{family}{Model family.}
  \item{data}{A data frame containing all variables that are specified in \code{formula} and \code{models}.}
  \item{scores}{A matrix of independent scores based on the true batch assignment.}
}

\examples{
  
  data(data, package="MetaGIM")

  # True batch assignment 
  M      <- ncol(data)-2
  bk     <- 10
  group2 <- split(1:M, ceiling((1:M)/bk))

  gim_add(models, group2, nsample, y~score, "case-control", data, score.add)
}

