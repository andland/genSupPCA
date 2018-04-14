# genSupPCA

This package implements two methods of supervised dimensionality reduction with data coming from exponential family distributions. These methods can be thought of as alternatives to partial least squares, especially when the covariates or response takes discrete values. 

1. The method from Rish, Irina, et al. "Closed-form supervised dimensionality reduction with generalized linear models." Proceedings of the 25th international conference on machine learning. ACM, 2008. [link](https://www.cs.cmu.edu/~ggordon/rish-etal-SDR.pdf). This is implemented in the function `genSupMF`.

1. Generalized supervised PCA, which extends generalized PCA from Landgraf, Andrew J., and Yoonkyung Lee. Generalized principal component analysis: Projection of saturated model parameters. Technical Report 892, Department of Statistics, The Ohio State University, 2015. [link](https://www.asc.ohio-state.edu/lee.2272/mss/tr892.pdf). This is implemented in the function `genSupPCA`.

For both methods, you must specify the distribution of both the covariates (`x`) and the univariate response (`y`). The objective is a linear combination of the deviance for reconstructing the covariate matrix and the deviance for predicting the responce variable. `alpha` controls the weight between the two deviances, with `alpha = 0` equivalent to a standard GLM with no dimensionality reduction and `alpha = Inf` equivalent to generalized principal component regression with no supervision in the dimensionality reduction. 

The difference between the two methods is how the dimensionality reduction is done for the covariate matrix. `genSupMF` approximates the natural parameters of the covariate matrix by $A B^T$, where $A$ and $B$ are rank $k$ matrices. `genSupPCA` approximates the natural parameters of the covariate matrix by $\tilde{\Theta} U U^T$, where $\tilde{\Theta}$ is the matrix of saturated natural parameters and $U$ is a $k$-dimensional orthonormal loading matrix. In both methods, the natural parameters of the response are predicted by the scores multiplied by a regression coefficient vector ($\beta$): $A \beta$ for `genSupMF` and $\tilde{\Theta} U \beta$ for `genSupPCA`.

## Installation

This package depends on another one of my packages ([`FOptM`](https://github.com/andland/FOptM)) for orthonormal optimization, which is on Github. It should be installed first before installing `genSupPCA`.

```r
# install.packages("devtools")
library(devtools)
install_github("andland/FOptM")
install_github("andland/genSupPCA")
```

## Caution

This package is a work in progress.
