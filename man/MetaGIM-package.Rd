\name{MetaGIM-package}
\alias{MetaGIM-package}
\docType{package}
\title{
  Integrative Analysis of Individual-Level Data and High-Dimensional Summary Statistics
}
\description{
  A general divide-and-conquer integration procedure, to combine high-dimensional summary data with individual-level for more efficient statistical inference.
}
\details{
  Researchers usually conduct statistical analyses based on models built on 
  raw data collected from individual participants (individual-level data). 
  There is a growing interest in enhancing inference efficiency by
  incorporating aggregated summary information from other sources, 
  such as summary statistics on genetic markers' marginal associations with 
  a given trait generated from genome-wide association studies. 
  However, combining high-dimensional summary data with individual-level data 
  using existing integrative procedures can be challenging. 
  This package overcomes this obstacle using a divide-and-conquer strategy
  by breaking the task into easier parallel jobs, 
  with each integrating a small proportion of summary data. 
}
\author{
Sheng Fu, Kai Yu
}
\references{
 Fu, S., Deng, L., Zhang, H.,  Wheeler, W., Qin, J., Yu, K. (2022)
 Integrative Analysis of Individual-Level Data and High-Dimensional Summary Statistics.
 Submitted

Zhang, H., Deng, L., Schiffman, M., Qin, J., Yu, K. (2020) 
  Generalized integration model for improved statistical inference by leveraging external summary data.\cr 
  Biometrika. asaa014, https://doi.org/10.1093/biomet/asaa014
}

\keyword{Generalized Integration Model, Case-control study, Divide-and-conquer scheme,
      Empirical likelihood, Estimating equations, High dimensionality, Summary data}

