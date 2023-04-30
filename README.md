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
Generate the data using the PCR setting (section 4.1 in the paper):
```{r}
library(JICO)

n1 = 50
n2 = 50
n = n1 + n2
p = 200
r = 1
r1 = 1
r2 = 1

alpha = rep(1, r)
alpha1 = rep(1, r1)
alpha2 = rep(1, r2)

set.seed(78)
X1 = MASS::mvrnorm(n1, rep(0, p), diag(p)) # covariates of the first group
X2 = MASS::mvrnorm(n2, rep(0, p), diag(p)) # covariates of the second group
X.list = list(X1, X2)
X = rbind(X1, X2)

q = r
q1 = r1
q2 = r2
V = matrix(svd(X)$v[,1:q], ncol = q)%*%rep(1/sqrt(q), q)
V1 = matrix(svd(X1%*%(diag(p) - V%*%t(V)))$v[,1:q1], ncol = q1)%*%rep(1/sqrt(q1), q1)
V2 = matrix(svd(X2%*%(diag(p) - V%*%t(V)))$v[,1:q2], ncol = q2)%*%rep(1/sqrt(q2), q2)

e1 = rnorm(n1)*0.2
e2 = rnorm(n2)*0.2
Y1 = X1%*%V%*%alpha + X1%*%V1%*%alpha1 + e1 # responses for the first group
Y2 = X2%*%V%*%alpha + X2%*%V2%*%alpha2 + e2 # responses for the second group
Y.list = list(Y1, Y2)

ml.JICO = continuum.multigroup.iter(
  X.list, Y.list, gam=1e12, rankJ=1, rankA=c(1, 1),
  maxiter = 300
)

cv.parameter.set = parameter.set.G_2(
  maxrankA = 1, maxrankJ = 1, gamma = 1e12
) # enumerate the set of tuning parameters

cv.ml.JICO = cv.continnum.iter(
  X.list, Y.list, parameter.set = cv.parameter.set, 
  criteria = "min", nfold = 5, maxiter = 300,
  center.X = F, scale.X = F, center.Y = F, scale.Y = F
) # fit the model and use CV to find the best parameters

# > cv.ml.JICO$parameter
# $rankA
# [1] 1 1
# 
# $rankJ
# [1] 1
# 
# $gam
# [1] 1e+12
```

Section 
