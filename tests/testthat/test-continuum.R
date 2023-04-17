test_that("Main JICO Algorithm", {
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
  
  expect_type(ml.JICO, "list")
})

test_that("Cross Validation", {
  set.seed(76)
  X1 = MASS::mvrnorm(50, rep(0, 200), diag(200)) # covariates of the first group
  X2 = MASS::mvrnorm(50, rep(0, 200), diag(200)) # covariates of the second group
  X.list = list(X1, X2)
  
  Y1 = matrix(rnorm(50)) # responses for the first group
  Y2 = matrix(rnorm(50)) # responses for the second group
  Y.list = list(Y1, Y2)
  
  cv.parameter.set = parameter.set.G_2(
    maxrankA = 1, maxrankJ = 1, gamma = 1e10
  ) # enumerate the set of tuning parameters
  
  cv.ml.JICO = cv.continnum.iter(
    X.list, Y.list, parameter.set = cv.parameter.set, 
    criteria = "min", nfold = 5, maxiter = 300
  )
  
  expect_type(cv.ml.JICO, "list")
})
