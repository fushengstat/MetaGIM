r# MetaGIM: A procedure to incorporate individual-level data and high-dimensional summary statistics 

# Overview
***MetaGIM*** is a novel method that integrates individual-level data and high-dimensional summary statistics using a divide-and-conquer scheme and meta-like procedure.

# Installation
***MetaGIM*** PACKAGE can be installed via Github. To install the latest version of MetaGIM package via Github, run following commands in R:
```{r, include = FALSE}
library(devtools)
devtools::install_github("fushengstat/MetaGIM")
```

# Usage
This is a good start point for using ***MetaGIM*** package.
```{r}
  data(data, package="MetaGIM")

  # Random assignment into 10 batches
  set.seed(0)
  M      <- ncol(data)-2
  bk     <- 10
  M0     <- sample(1:M)
  group1 <- split(M0, ceiling((1:M)/bk))

  metagim(models, group1, nsample, y~score, "case-control", data)
```

Please see the [MetaGIM user manual](https://github.com/fushengstat/MetaGIM/blob/main/doc/MetaGIM-manual.pdf). 
 

<!---
# Information
Author: Han Zhang, Kai Yu, Sheng Fu \
Maintainer: Bill Wheeler <wheelerb@imsweb.com>
--->


# References
Fu S, Deng L, Zhang H, Wheeler W, Qin J, Yu, K (2023). Integrative analysis of individual-level data and high-dimensional summary statistics. Manuscript. \
Zhang, H., Deng, L., Schiffman, M., Qin, J., & Yu, K. (2020). Generalized integration model for improved statistical inference by leveraging external summary data. Biometrika, 107(3), 689-703. \
Zhang, H., Deng, L., Wheeler, W., Qin, J., & Yu, K. (2022). Integrative analysis of multiple case‐control studies. Biometrics, 78(3), 1080-1091.

