# MetaGIM: A procedure to incorporate individual-level data and high-dimensional summary statistics 

# Overview
**MetaGIM** is a novel method that integrates individual-level data and high-dimensional summary statistics using a divide-and-conquer scheme and meta-like procedure.

# Installation
**MetaGIM** PACKAGE can be installed via Github. To install the latest version of MetaGIM package via Github, run following commands in R:
```{r, include = FALSE}
library(devtools)
devtools::install_github("fushengstat/MetaGIM")
```

# Usage
This is a good start point for using **MetaGIM** package.
```{r,include = FALSE}
data(data, package="MetaGIM")

# Random assignment into 10 batches
set.seed(0)
M      <- ncol(data)-2                             # M is the number of markers
bk     <- 10                                       # bk is the bathc size of each block
group <- split(sample(1:M), ceiling((1:M)/bk))     # group is the batch assignment of M markers

# models is the collection of all smmmary data
# group is a batche assignment
# nsample is the sample sizes of case and controls
metagim(models, group, nsample, y~score, "case-control", data)
```

Please see the [MetaGIM user manual](https://github.com/fushengstat/MetaGIM/blob/main/doc/MetaGIM-manual.pdf) for more details. 
 

<!---
# Information
Author: Han Zhang, Kai Yu, Sheng Fu \
Maintainer: Bill Wheeler <wheelerb@imsweb.com>
--->


# References
Fu, S., Deng, L., Zhang, H., Qin, J., & Yu, K. (2023). <a href="https://doi.org/10.1093/bioinformatics/btad156" target="_blank">Integrative analysis of individual-level data and high-dimensional summary statistics</a>. Bioinformatics, 39(4), btad156. \
Zhang, H., Deng, L., Schiffman, M., Qin, J., & Yu, K. (2020). Generalized integration model for improved statistical inference by leveraging external summary data. Biometrika, 107(3), 689-703. \
Zhang, H., Deng, L., Wheeler, W., Qin, J., & Yu, K. (2022). Integrative analysis of multiple case‐control studies. Biometrics, 78(3), 1080-1091.

