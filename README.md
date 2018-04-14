# genSupPCA

This package implements two methods of supervised dimensionality reduction with data coming from exponential family distributions. 

1. The method from Rish, Irina, et al. "Closed-form supervised dimensionality reduction with generalized linear models." Proceedings of the 25th international conference on Machine learning. ACM, 2008.

1. Generalized supervised PCA, which extends generalized PCA from Landgraf, Andrew J., and Yoonkyung Lee. Generalized principal component analysis: Projection of saturated model parameters. Technical Report 892, Department of Statistics, The Ohio State University, 2015. 

This package is a work in progress.

## Installation

This package depends on another one of my packages (`FOptM`) for orthonormal optimization, which is on Github. It should be installed first before installing `genSupPCA`.

```r
# install.packages("devtools")
library(devtools)
install_github("andland/FOptM")
install_github("andland/genSupPCA")
```

