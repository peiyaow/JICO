# JICO
R package with implementations of the JICO algorithm: https://arxiv.org/abs/2209.12388

# How to Install
You can install the package through use devtools:
```r
devtools::install_github("peiyaow/JICO", upgrade_dependencies = FALSE)
```

# Overview
To see the full list of exported functions:

```{r}
library("JICO")
ls("package:JICO")
```

A quick overview of some of the key functions:

* `continuum.multigroup.iter`: This function iteratively solves the multi-group regression problem using the JICO algorithm

* `cv.continnum.iter`: This function performs K-fold cross validations to select the best tuning parameters for JICO.

# Toy Example

```{r}
library(JICO)

set.seed(76)
X1 = MASS::mvrnorm(50, rep(0, 200), diag(200)) # covariates of the first group
X2 = MASS::mvrnorm(50, rep(0, 200), diag(200)) # covariates of the second group
X.list = list(X1, X2)

Y1 = matrix(rnorm(50)) # responses for the first group
Y2 = matrix(rnorm(50)) # responses for the second group
Y.list = list(Y1, Y2)

ml.JICO = continuum.multigroup.iter(
  X.list, Y.list, gam=1e10, rankJ=1, rankA=c(1, 1), 
  maxiter = 300
)

cv.parameter.set = parameter.set.G_2(
  maxrankA = 1, maxrankJ = 1, gamma = 1e10
) # enumerate the set of tuning parameters

cv.ml.JICO = cv.continnum.iter(
  X.list, Y.list, parameter.set = cv.parameter.set, 
  criteria = "min", nfold = 5, maxiter = 300
) # fit the model and use CV to find the best parameters
```
