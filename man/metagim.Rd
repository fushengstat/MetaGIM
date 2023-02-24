\name{metagim}
\alias{metagim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  MetaGIM
}
\description{Integrating High-Dimensional Summary Statistics}
\usage{
metagim(models,group,sample.info,form,family,data)
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
}

\examples{
  
  data(data, package="MetaGIM")

  # Random assignment into 10 batches
  set.seed(0)
  M      <- ncol(data)-2
  bk     <- 10
  M0     <- sample(1:M)
  group1 <- split(M0, ceiling((1:M)/bk))

  metagim(models, group1, nsample, y~score, "case-control", data)

}