#' Generate parameter sets (G=2)
#'
#' This function generate set of hyperparameters when there 
#' are two groups. 
#'
#' @param maxrankA The maximum rank for individual component
#' @param maxrankJ The maximum rank for joint component
#' @param gamma The gamma parameter. Need to be fixed. 
#' @return A list of hyperparameter candidates
#' @export
parameter.set.G_2 = function(maxrankA, maxrankJ, gamma){
  parameter.set = list()
  for(rankA1 in 0:maxrankA)
    for(rankA2 in 0:maxrankA)
      for (rankJ in 0:maxrankJ){
        parameter.set = rlist::list.append(parameter.set, list(rankA = c(rankA1, rankA2), rankJ = rankJ, gam = gamma))
      }
  return(parameter.set)
}

#' Helper function to compute the inverse of input X matrix
#' 
#' This function computes the general inverse of X when it exists. 
#' If X contains a degenerated dimension, return the original X.
#'
#' @param x The input matrix X
#' @return Either the general inverse of X or the X itself
SOLVE = function(x){
  if (sum(dim(x))){
    return(MASS::ginv(x))
  }else{
    return(x)
  }
}

#' Generate parameter sets (G=3)
#'
#' This function generate set of hyperparameters when there 
#' are three groups. 
#'
#' @param maxrankA The maximum rank for individual component
#' @param maxrankJ The maximum rank for joint component
#' @param gamma The gamma parameter. Need to be fixed. 
#' @return A list of hyperparameter candidates
#' @export
parameter.set.G_3 = function(maxrankA, maxrankJ, gamma){
  parameter.set = list()
  for(rankA1 in 0:maxrankA)
    for(rankA2 in 0:maxrankA)
      for (rankA3 in 0:maxrankA){
        for (rankJ in 0:maxrankJ){
          parameter.set = rlist::list.append(parameter.set, list(rankA = c(rankA1, rankA2, rankA3), rankJ = rankJ, gam = gamma))
        }
      }
  return(parameter.set)
}

#' Generate parameter sets (equal individual ranks)
#'
#' This function generate set of hyperparameters when the individual ranks are the same
#'
#' @param G number of groups
#' @param maxrankA The maximum rank for individual component
#' @param maxrankJ The maximum rank for joint component
#' @param gamma.list The list of candidate gammas to be tuned
#' @return A list of hyperparameter candidates
#' @export
parameter.set.rankA_eq = function(G, maxrankA, maxrankJ, gamma.list){
  # generate set of hyperparameters when the individual ranks are the same
  
  # G: number of groups
  # maxrankA: the maximum rank for individual component
  # maxrankJ: the maximum rank for joint component
  # 
  parameter.set <- list() 
  for (gam in gamma.list)
    for(rankA in 0:maxrankA)
      for (rankJ in 0:maxrankJ){
        parameter.set = rlist::list.append(parameter.set, list(rankA = rep(rankA, G), rankJ = rankJ, gam = gam))
      }
  return(parameter.set)
}

#' Generate diagonal matrix
#'
#' This function returns a diagnoal matrix using the input vector or number as diagonal.
#'
#' @param e Diagonal element. Can be a vector or a number
#' @return A square diagonal matrix using the input as diagonal elements
DIAG = function(e){
  if (length(e) > 1){
    return(diag(e))
  }else{
    return(matrix(e))
  }
}

#' Utility function to create folds for stratified samples
#'
#' This function generate data folds for cross validation given stratified samples
#'
#' @param strat_id A vector of the stratified sample id. 
#' E.g. In total of 5 samples, first three from group 1, last two from group 2
#' -> c(1, 1, 1, 2, 2)
#' @param k Number of folds to create.
#' @return A list of sample indices in k folds.
createFolds <- function(strat_id, k) {
  # function to create folds for cross validation
  
  if(k > length(strat_id)) {
    k <- length(strat_id)
  }	
  perm <- sample(length(strat_id))
  strat_id <- factor(strat_id, levels=sample(unique(strat_id)))
  
  strat_order <- order(strat_id[perm])
  
  num_item_per_fold_ceil <- ceiling(length(strat_id) / k)
  
  fold_ids <- rep(seq_len(k), times= num_item_per_fold_ceil)
  fold_ids <- fold_ids[seq_along(strat_id)]
  
  folds <- split(perm[strat_order], fold_ids)
  names(folds) <- paste0("Fold", seq_len(k))	
  return(folds)
}

