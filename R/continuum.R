#' Compute the coefficients from the continuum regression (CR) algorithm
#'
#' This function converts the CR algorithm outputs to the regression coefficients 
#'
#' @param X The input feature matrix
#' @param Y The input response vector
#' @param C The weight matrix computed from CR algorithm
#' @param lambda Deprecated. Regularization parameter if L2 penalization is used for CR. 
#' JICO uses zero as default.
#' @return A list of regression coefficients to perform the prediction task.
C2beta = function(X, Y, C, lambda){
  n = nrow(X)
  X.mean = apply(X, 2, mean)
  Y.mean = mean(Y)
  S = t(X)%*%X
  s = t(X)%*%Y
  
  C = matrix(C, nrow = ncol(X))
  om = ncol(C) #omega
  
  alpha.C = SOLVE(t(C)%*%S%*%C + n*lambda*diag(om))%*%t(C)%*%s
  beta.C = C%*%alpha.C
  intercept = Y.mean - X.mean%*%beta.C
  
  return(list(intercept = intercept, beta = beta.C, alpha = alpha.C, coef = matrix(c(intercept, beta.C))))
}

#' Helper function to compute the SVD results 
#' 
#' This function computes the SVD results from a given matrix X. 
#' This is used as the initialization for the continuum regression.
#'
#' @param X The input feature matrix
#' @return A list of SVD results that are served as CR algorithm's inputs.
initialize.UDVZ = function(X){
  svd.X = svd(X)
  d = svd.X$d
  V = svd.X$v
  U = svd.X$u
  m = Matrix::rankMatrix(X)[1]
  
  if (m > 0){
    e = (d^2)[1:m]
    D = DIAG(d[1:m])
    E = DIAG(e)
    V = V[,1:m]
    U = U[,1:m]
  }else{
    U = matrix(,nrow=nrow(X), ncol=m)
    D = matrix(,nrow=m, ncol=m)
    E = matrix(,nrow=m, ncol=m)
    V = matrix(,nrow=ncol(X), ncol=m)
    e=0
  }
  Z = matrix(,nrow=m, ncol=0)
  return(list(U=as.matrix(U), D=D, V=as.matrix(V), Z=Z, E=E, e=e, m=m))
}

#' The continuum regression (CR) algorithm
#'
#' This function performs an iteration update of the JICO algorithm using the CR algorithm. 
#' Details can be found in Appendix B in the JICO paper: https://arxiv.org/pdf/2209.12388.pdf
#'
#' @param X The input feature matrix
#' @param Y The input response vector
#' @param lambda Deprecated. Regularization parameter if L2 penalization is used for CR. 
#' JICO uses zero as default.
#' @param gam The gamma parameter in the CR algorithm. Set gam=0 for OLS model, gam=0.5 for PLS model, 
#' gam >= 1e10 for PCR model
#' @param om The desired number of weight vectors to obtain in the CR algorithm, i.e. the predefined rank of 
#' joint or individual componenet
#' @param U_old The given inputs U from the previous JICO iteration step
#' @param D_old The given inputs D from the previous JICO iteration step
#' @param V_old The given inputs V from the previous JICO iteration step
#' @param Z_old The given inputs Z from the previous JICO iteration step
#' @param verbose Boolean. If it's desired to print out intermediate outputs
#' @return A list of CR outputs that serve as the input for the next JICO iteration
#' @export
continuum = function(X, Y, lambda, gam, om, 
                     U_old=matrix(,nrow=nrow(X), ncol=0), D_old=matrix(,nrow=0, ncol=0), V_old=matrix(,nrow=0, ncol=0), Z_old=matrix(,nrow=0, ncol=0), 
                     verbose = FALSE){
  # ----------- initialization ----------- #
  n = nrow(X)
  p = ncol(X)
  
  s = t(X)%*%Y
  UDVZ = initialize.UDVZ(X)
  D = UDVZ$D
  E = UDVZ$E
  V = UDVZ$V
  U = UDVZ$U
  e = UDVZ$e
  m = UDVZ$m
  d = t(V)%*%s
  E2 = D%*%t(U)%*%U_old%*%D_old
  
  # default to 1
  tau = 1
  
  # ----------- find a good initial value for rho ----------- #
  B = E2%*%Z_old
  fn = function(rho){
    A = diag(tau^2*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*tau^2*E
    M = solve(A)- solve(A)%*%B%*%SOLVE(t(B)%*%solve(A)%*%B)%*%t(B)%*%solve(A)
    q = tau*rho*d
    Mq = M%*%q
    z = Mq/norm(Mq, "2")
    return(t(z)%*%E%*%z + n*lambda - rho)
  }
  nleqslv.res = nleqslv::nleqslv(e[1]+n*lambda, fn, method = "Newton", global = "none", control = list(maxit = 150))
  rho = nleqslv.res$x
  if (nleqslv.res$termcd != 1 && verbose){
    print(paste0("Warning! The value is ", as.character(fn(rho))))
    print(nleqslv.res$termcd)
  }
  
  # ----------- iteration on the first step ----------- #
  A = diag(tau^2*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*tau^2*E
  M = solve(A)- solve(A)%*%B%*%SOLVE(t(B)%*%solve(A)%*%B)%*%t(B)%*%solve(A)
  q = tau*rho*d
  Mq = M%*%q
  Z = Mq/norm(Mq, "2")
  
  E_all = cbind(E2, E)
  Z_all = as.matrix(Matrix::bdiag(Z_old, Z))
  B = E_all%*%Z_all
  rho0 = rho
  
  # ----------- iteration on the following steps if om > 1 ----------- #
  while (ncol(Z) < om){
    nleqslv.res = nleqslv::nleqslv(rho0, fn, method = "Newton", global = "none", control = list(maxit = 150))
    rho = nleqslv.res$x
    if (nleqslv.res$termcd != 1 && verbose){
      print(paste0("Warning! The value is ", as.character(fn(rho))))
      print(nleqslv.res$termcd)
    }
    
    A = diag(tau^2*(gam*rho-(gam-1)*(n*lambda)), m) + (1-gam)*tau^2*E
    M = solve(A)- solve(A)%*%B%*%SOLVE(t(B)%*%solve(A)%*%B)%*%t(B)%*%solve(A)
    q = tau*rho*d
    Mq = M%*%q
    z = Mq/norm(Mq, "2")
    Z = cbind(Z, z)
    Z_all = as.matrix(Matrix::bdiag(Z_old, Z))
    B = E_all%*%Z_all
    rho0 = rho
  }
  
  # ----------- compute final results ----------- 
  C = V%*%Z
  C = C[, 0:min(ncol(C), om)]
  a = V%*%MASS::ginv(E)^(1/2)%*%Z
  
  return(list(C = as.matrix(C), a = a, V = as.matrix(V), Z = Z, E = E, D = D, U = as.matrix(U)))
}

#' The Joint and Individual Component Regression (JICO) algorithm
#'
#' This function iteratively solves the multi-group regression problem using the JICO algorithm. 
#' JICO paper: https://arxiv.org/pdf/2209.12388.pdf
#'
#' @param X.list The list of feature matrices from multiple groups.
#' @param Y.list The list of feature vectors from multiple groups.
#' @param lambda Deprecated. Regularization parameter if L2 penalization is used for CR. 
#' JICO uses zero as default.
#' @param gam The gamma parameter in the CR algorithm. Set gam=0 for OLS model, gam=0.5 for PLS model, 
#' gam >= 1e10 for PCR model.
#' @param rankJ The rank for the joint component.
#' @param rankA The ranks for individual components.
#' @param maxiter The maximum number of iterations to conduct before algorithm convergence.
#' @param conv The tolerance level for covergence.
#' @param center.X Boolean. If X should be preprocessed with centralization.
#' @param scale.X Boolean. If X should be preprocessed with scaling.
#' @param center.Y Boolean. If Y should be preprocessed with centralization.
#' @param scale.Y Boolean. If Y should be preprocessed with scaline.
#' @param orthIndiv Boolean. If we impose the orthogonality constraint on individual components.
#' @param I.initial The initial values for individual components.
#' @param sd The standard deviation used to generate random initial values for individual weight vectors.
#' @return The estimated parameters from JICO.
#' @examples 
#' set.seed(76)
#' X1 = MASS::mvrnorm(50, rep(0, 200), diag(200)) # covariates of the first group
#' X2 = MASS::mvrnorm(50, rep(0, 200), diag(200)) # covariates of the second group
#' X.list = list(X1, X2)
#'
#' Y1 = matrix(stats::rnorm(50)) # responses for the first group
#' Y2 = matrix(stats::rnorm(50)) # responses for the second group
#' Y.list = list(Y1, Y2)
#'
#' ml.JICO = continuum.multigroup.iter(
#'   X.list, Y.list, gam=1e10, rankJ=1, rankA=c(1, 1), 
#'   maxiter = 300
#' )
#' @export
continuum.multigroup.iter = function(X.list, Y.list, lambda = 0, gam, rankJ, rankA, maxiter = 1000, conv = 1e-7, 
                                     center.X = TRUE, scale.X = TRUE, center.Y = TRUE, scale.Y = TRUE, orthIndiv = FALSE,
                                     I.initial = NULL, sd = 0){
  G = length(X.list)
  centerValues.X <- list()
  scaleValues.X <- list()
  centerValues.Y <- list()
  scaleValues.Y <- list()
  n = c()
  for (g in 1:G){
    n[g] = nrow(X.list[[g]])
  }
  N = sum(n)
  p = ncol(X.list[[1]])
  
  n_cumsum = c(0,cumsum(n))
  index.list = lapply(1:G, function(i) (n_cumsum[i]+1):n_cumsum[i+1] )
  
  for (g in 1:G){
    # X
    if (center.X){
      centerValues.X[[g]] = apply(X.list[[g]], 2, mean)
    }else{
      centerValues.X[[g]] = rep(0, p)
    }
    if (scale.X){
      scaleValues.X[[g]] = norm(X.list[[g]], type = "f")
    }else{
      scaleValues.X[[g]] = 1
    }
    X.list[[g]] = sweep(X.list[[g]], 2, centerValues.X[[g]])
    X.list[[g]] = X.list[[g]]/scaleValues.X[[g]]
    
    # Y
    if (center.Y){
      centerValues.Y[[g]] = mean(Y.list[[g]])
    }else{
      centerValues.Y[[g]] = 0
    }
    if (scale.Y){
      scaleValues.Y[[g]] = norm(Y.list[[g]], type = "f")
    }else{
      scaleValues.Y[[g]] = 1
    }
    Y.list[[g]] = sweep(matrix(Y.list[[g]]), 2, centerValues.Y[[g]])
    Y.list[[g]] = sweep(matrix(Y.list[[g]]), 2, scaleValues.Y[[g]], FUN = "/")
  }
  X = do.call(rbind, X.list)
  Y = do.call(rbind, Y.list)
  
  # when gam == 0, it is equivalent to OLS. Either rankJ or rankA needs to be 0 and at most to be 1
  if (gam == 0){
    rankJ = min(1, rankJ)
    if (rankJ){
      rankA = rep(0, G)
    }else{
      rankA = sapply(rankA, function(r) min(1, r))
    }
  }
  
  nrun = 0
  converged = F
  
  Y.homo = Y
  Y.homo.list = Y.list
  
  C = matrix(stats::rnorm(p*rankJ, 0, sd), p, rankJ)
  Cind = lapply(1:G, function(g) matrix(stats::rnorm(p*rankA[g], 0, sd), p, rankA[g]))
  Cind_tot = do.call(cbind, Cind)
  P = Cind_tot%*%SOLVE(t(Cind_tot)%*%Cind_tot)%*%t(Cind_tot)
  
  if (is.null(I.initial)){
    X.heter = X%*%P
    X.heter.list = lapply(index.list, function(ixs) as.matrix(X.heter[ixs,]))
    X.homo.list = lapply(1:G, function(g) X.list[[g]] - X.heter.list[[g]])
  }else{
    X.heter.list = I.initial
    X.homo.list = lapply(1:G, function(g) X.list[[g]] - X.heter.list[[g]])
  }
  
  R = X
  R[,] = 0
  
  r = Y
  r[,] = 0
  
  U = list()
  W = list()
  
  ct.homo = matrix(0, nrow = rankJ, ncol = 1)
  ct.heter = lapply(1:G, function(g) matrix(0, nrow = rankA[g], ncol = 1))
  
  UDVZ.heter.list = lapply(1:G, function(g) initialize.UDVZ(X.heter.list[[g]]))
  U.heter.list = lapply(1:G, function(g) UDVZ.heter.list[[g]]$U)
  Z.heter.list = lapply(1:G, function(g) UDVZ.heter.list[[g]]$Z)
  V.heter.list = lapply(1:G, function(g) UDVZ.heter.list[[g]]$V)
  D.heter.list = lapply(1:G, function(g) UDVZ.heter.list[[g]]$D)
  
  U.heter = as.matrix(do.call(Matrix::bdiag, U.heter.list))
  Z.heter = as.matrix(do.call(Matrix::bdiag, Z.heter.list))
  V.heter = do.call(rbind, V.heter.list)
  D.heter = as.matrix(do.call(Matrix::bdiag, D.heter.list))
  
  X.homo = do.call(rbind, X.homo.list)
  UDVZ.homo = initialize.UDVZ(X.homo)
  U.homo = UDVZ.homo$U
  Z.homo = UDVZ.homo$Z
  V.homo = UDVZ.homo$V
  D.homo = UDVZ.homo$D
  U.homo.list = lapply(index.list, function(ixs) as.matrix(U.homo[ixs,]))
  
  while (nrun < maxiter & !converged){
    # initialization
    rlast = r
    Rlast = R
    
    ct.homo.last = ct.homo
    ct.heter.last = ct.heter
    
    X.homo.list = lapply(1:G, function(g) X.list[[g]] - X.heter.list[[g]])
    X.homo = do.call(rbind, X.homo.list)
    
    # joint
    if (rankJ){
      ml.homo = continuum(X.homo%*%(diag(p) - P), Y.homo, lambda = lambda, gam = gam, om = rankJ,
                          U_old=U.heter, D_old=D.heter, V_old=V.heter, Z_old=Z.heter)
      C = ml.homo$C
      
      U.homo = ml.homo$U
      Z.homo = ml.homo$Z
      V.homo = ml.homo$V
      D.homo = ml.homo$D
      U.homo.list = lapply(index.list, function(ixs) as.matrix(U.homo[ixs,]))
    }
    U = rlist::list.append(U, C)
    
    XC = X.homo%*%C
    ct.homo = (diag(t(XC)%*%XC)^(gam-1))*(t(XC)%*%Y.homo)^2
    
    beta.C = C2beta(X.homo, Y.homo, C, lambda = lambda)$beta
    Yhat.homo.list = lapply(1:G, function(g) X.homo.list[[g]]%*%beta.C)
    Y.heter.list = lapply(1:G, function(g) Y.list[[g]] - Yhat.homo.list[[g]])
    X.homo.list = lapply(1:G, function(g) X.homo.list[[g]]%*%C%*%SOLVE(t(C)%*%C)%*%t(C))
    X.heter.list = lapply(1:G, function(g) X.list[[g]] - X.homo.list[[g]])
    
    # individual
    temp = X.heter.list
    U.heter.list = list()
    Z.heter.list = list()
    V.heter.list = list()
    D.heter.list = list()
    for (g in 1:G){
      tempC = C
      # orthogonalization
      if (orthIndiv){
        if (nrun > 0){
          for (j in (1:G)[-g]){
            tempC = cbind(tempC, Cind[[j]])
          }
        }
      }
      temp[[g]] <- temp[[g]]%*%(diag(p) - tempC%*%SOLVE(t(tempC)%*%tempC)%*%t(tempC))
      
      ml.heter = continuum(temp[[g]], Y.heter.list[[g]], lambda = lambda, gam = gam, om = rankA[g], U_old=U.homo.list[[g]], D_old=D.homo, V_old=V.homo, Z_old=Z.homo)
      Cind[[g]] = ml.heter$C
      
      U.heter.list[[g]] = ml.heter$U
      Z.heter.list[[g]] = ml.heter$Z
      V.heter.list[[g]] = ml.heter$V
      D.heter.list[[g]] = ml.heter$D
    }
    
    U.heter = as.matrix(do.call(Matrix::bdiag, U.heter.list))
    Z.heter = as.matrix(do.call(Matrix::bdiag, Z.heter.list))
    V.heter = do.call(rbind, V.heter.list)
    D.heter = as.matrix(do.call(Matrix::bdiag, D.heter.list))
    
    W = rlist::list.append(W, Cind)
    Cind_tot = do.call(cbind, Cind)
    P = Cind_tot%*%SOLVE(t(Cind_tot)%*%Cind_tot)%*%t(Cind_tot)
    
    XC.list = lapply(1:G, function(g) X.heter.list[[g]]%*%Cind[[g]])
    
    ct.heter = lapply(1:G, function(g) 
      diag((t(XC.list[[g]])%*%XC.list[[g]])^(gam-1))*(t(XC.list[[g]])%*%Y.heter.list[[g]])^2)
    
    beta.Cind = lapply(1:G, function(g) C2beta(X.heter.list[[g]], Y.heter.list[[g]], Cind[[g]], lambda)$beta)
    Yhat.heter.list = lapply(1:G, function(g) X.heter.list[[g]]%*%beta.Cind[[g]])
    
    r.list = lapply(1:G, function(g) Y.heter.list[[g]] - Yhat.heter.list[[g]])
    r = do.call(rbind, r.list)
    
    Y.homo.list = lapply(1:G, function(g) Y.list[[g]] - Yhat.heter.list[[g]])
    Y.homo = do.call(rbind, Y.homo.list)
    X.heter.list = lapply(1:G, function(g) X.heter.list[[g]]%*%Cind[[g]]%*%SOLVE(t(Cind[[g]])%*%Cind[[g]])%*%t(Cind[[g]]))
    
    
    # compute residuals
    R.list = lapply(1:G, function(g) X.list[[g]] - X.homo.list[[g]] - X.heter.list[[g]])
    R = do.call(rbind, R.list)
    
    if (gam > 1e5){
      if (norm(Rlast - R, type = "f") <= conv){
        converged <- T
      }
    }else{
      CT1 = (ct.homo - ct.homo.last)/ct.homo.last
      CT2 = lapply(1:G, function(g) as.vector((ct.heter[[g]]-ct.heter.last[[g]])/ct.heter.last[[g]]))
      CT = c(CT1, do.call(c, CT2))
      if (max(abs(CT), na.rm = T) <= conv){
        converged <- T
      }
    }
    
    nrun = nrun + 1
  }
  if (converged) {
    cat(paste("Algorithm converged after ", nrun, 
              " iterations.\n"))
  }
  else {
    cat(paste("Algorithm did not converge after ", 
              nrun, " iterations.\n"))
  }
  
  beta.C0 = beta.C
  beta.Cind0 = beta.Cind
  beta.C = list()
  beta.Cind = list()
  intercept = list()
  for (g in 1:G){
    beta.C[[g]] = beta.C0/scaleValues.X[[g]]*scaleValues.Y[[g]]
    beta.Cind[[g]] = beta.Cind0[[g]]/scaleValues.X[[g]]*scaleValues.Y[[g]]
    intercept[[g]] = centerValues.Y[[g]] - t(beta.C[[g]])%*%centerValues.X[[g]] - t(beta.Cind[[g]])%*%centerValues.X[[g]]
  }
  return(list(C = C, Cind = Cind, 
              U = U, W = W, 
              centerValues.X = centerValues.X, scaleValues.X = scaleValues.X, 
              centerValues.Y = centerValues.Y, scaleValues.Y = scaleValues.Y,
              intercept = intercept, beta.C = beta.C, beta.Cind = beta.Cind, 
              beta = beta.C0, beta_i = beta.Cind0,
              J = X.homo.list, I = X.heter.list, 
              R = R, r = r,
              converged = converged, nrun = nrun))
}

#' Fit JICO with cross-validation to tune hyperparameters
#'
#' This function performs K-fold cross validations to select the best tuning parameters for JICO.
#'
#' @param X.list The list of feature matrices from multiple groups.
#' @param Y.list The list of feature vectors from multiple groups.
#' @param lambda Deprecated. Regularization parameter if L2 penalization is used for CR. 
#' JICO uses zero as default.
#' @param parameter.set The set of parameters to be tuned on. Containing choices of rankJ, rankA and gamma.
#' @param nfolds number of folds to perform CV
#' @param maxiter The maximum number of iterations to conduct before algorithm convergence.
#' @param center.X Boolean. If X should be preprocessed with centralization.
#' @param scale.X Boolean. If X should be preprocessed with scaling.
#' @param center.Y Boolean. If Y should be preprocessed with centralization.
#' @param scale.Y Boolean. If Y should be preprocessed with scaline.
#' @param orthIndiv Boolean. If we impose the orthogonality constraint on individual components.
#' @param plot Boolean. If we want to plot the rMSE vs different parameters
#' @param criteria criteria for selecting the best parameter. 
#' Use "min" to choose the parameter giving the best performance. 
#' Use "1se" to choose the simplest model that gives performance within 1se from the best one.
#' @param sd The standard deviation used to generate random initial values for individual weight vectors.
#' @return The parameter from the parameter.set that fit the training data the best.
#' @examples 
#' set.seed(76)
#' X1 = MASS::mvrnorm(50, rep(0, 200), diag(200)) # covariates of the first group
#' X2 = MASS::mvrnorm(50, rep(0, 200), diag(200)) # covariates of the second group
#' X.list = list(X1, X2)
#'
#' Y1 = matrix(stats::rnorm(50)) # responses for the first group
#' Y2 = matrix(stats::rnorm(50)) # responses for the second group
#' Y.list = list(Y1, Y2)
#' 
#' cv.parameter.set = parameter.set.G_2(
#'    maxrankA = 1, maxrankJ = 1, gamma = 1e10
#' ) # enumerate the set of tuning parameters
#' 
#' cv.ml.JICO = cv.continnum.iter(
#'   X.list, Y.list, parameter.set = cv.parameter.set, 
#'   criteria = "min", nfold = 5, maxiter = 300
#' ) # fit the model and use CV to find the best parameters
#' @export
cv.continnum.iter = function(X.list, Y.list, lambda = 0, parameter.set, nfolds = 10, maxiter = 100,
                             center.X = TRUE, scale.X = TRUE, center.Y = TRUE, scale.Y = TRUE, orthIndiv = FALSE, 
                             plot = F, criteria = c("min", "1se"), sd = 0){
  G = length(X.list)
  flds.list = lapply(1:G, function(g) createFolds(Y.list[[g]], k = nfolds))
  MSE.list = list()
  for (k in 1:nfolds){
    X.train.list = lapply(1:G, function(g) X.list[[g]][unlist(flds.list[[g]][-k]), ])
    X.val.list = lapply(1:G, function(g) X.list[[g]][unlist(flds.list[[g]][k]), ])
    Y.train.list = lapply(1:G, function(g) matrix(Y.list[[g]][unlist(flds.list[[g]][-k])]))
    Y.val.list = lapply(1:G, function(g) matrix(Y.list[[g]][unlist(flds.list[[g]][k])]))
    Y.val = do.call(rbind, Y.val.list)
    
    ml.list = lapply(parameter.set, function(parameter) 
      continuum.multigroup.iter(X.train.list, Y.train.list, maxiter = maxiter,     
                                gam = parameter$gam, rankJ = parameter$rankJ, rankA = parameter$rankA, 
                                center.X = center.X, scale.X = scale.X, center.Y = center.Y, scale.Y = scale.Y, orthIndiv = orthIndiv,
                                sd = sd))
    Yhat.list = lapply(ml.list, function(ml) 
      do.call(rbind, lapply(1:G, function(g) as.numeric(ml$intercept[[g]]) + X.val.list[[g]]%*%ml$beta.C[[g]] + X.val.list[[g]]%*%ml$beta.Cind[[g]])))
    MSE.list[[k]] = sapply(Yhat.list, function(Yhat) mean((Y.val - Yhat)^2))
  }
  MSE = do.call(rbind, MSE.list)
  rMSE = apply(MSE, 2, function(x) mean(sqrt(x)))
  if (plot){
    plot(rMSE)
  }
  if (criteria == "1se"){
    absBest = min(rMSE)
    MSEsd = apply(MSE, 2, function(x) sd(sqrt(x)))/sqrt(nfolds)
    
    gam_list = c()
    for (para in parameter.set){
      gam_list = c(gam_list, para$gam)
    }
    num_gam = length(unique(gam_list))
    
    rMSE.mtx = matrix(rMSE, ncol = num_gam) # num of cols = number of gams 
    absBest.ix = which.min(rMSE.mtx)
    col.ix = ceiling(absBest.ix/nrow(rMSE.mtx))
    row.ix = min(which((rMSE.mtx[, col.ix] - MSEsd[absBest.ix]) < absBest))
    ix = col.ix + (row.ix-1)*(num_gam)
    parameter = parameter.set[[ix]]
  }
  if (criteria == "min"){
    ix = which.min(rMSE)
    parameter = parameter.set[[ix]]
  }
  return(list(rMSE = rMSE, MSE = MSE, ix = ix, parameter = parameter))
}
