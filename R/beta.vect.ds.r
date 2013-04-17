#' Function required by the function \code{datashield.gee}
#'
#' @title Generates combined final estimates and standard errors of the fitted model
#'
#' @param score
#' @param infomatrix
#' @param J.matrix
#' @param start.betas
#' @return a list that contains the beta estimates and standard errors
#'
beta.vect.ds <- function(score, infomatrix, J.matrix, start.betas){
  if(is.null(start.betas)) {
    start.betas <- rep(0, length(score[[1]]))
  }
  
  sum.scores <- as.list(rep(0, length(start.betas)))
  rows <- dim(infomatrix[[1]])[1]
  cols <- dim(infomatrix[[1]])[2]
  sum.I.matrices <- matrix(0, nrow=rows, ncol=cols)
  
  # loop through studies and sum up matrices and vectors separately
  for(i in 1:length(start.betas)){
    for(s in 1:length(score)){
      sum.scores[[i]] <- sum.scores[[i]] + score[[s]][i]
    }
  }
  
  for(i in 1: length(infomatrix)){
    sum.I.matrices <- sum.I.matrices + infomatrix[[i]]
  }
  
  # add up start values for beta, score vector and info matrix
  beta.vector <- start.betas + (solve(sum.I.matrices)%*%unlist(sum.scores))
  
  # compute the 'empirical' standard errors of the beta values (get the J.matrix first)
  rows <- dim(J.matrix[[1]])[1]
  cols <- dim(J.matrix[[1]])[2]
  sum.J.matrices <- matrix(0, nrow=rows, ncol=cols)
  for(i in 1: length(J.matrix)){
    sum.J.matrices <- sum.J.matrices + J.matrix[[i]]
  }
  v.beta.hat <- solve(sum.I.matrices)%*%sum.J.matrices%*%solve(sum.I.matrices)
  stderrs <- sqrt(diag(v.beta.hat))
  
  # output
  return(list(beta.vector, stderrs))
}