#' This function is required by the function \code{datashield.gee}
#'
#' @title Combines some estimates (alpha and phi) obtained from the different studies
#'
#' @param N 
#' @param npara
#' @param M_study
#' @param alphaM
#' @param sum_p
#' @return a list that contains alpha and phi
#'
alpha.phi.ds <- function(N, npara, M_study, alphaM, sum_p){
  N <- sum(N)
  npara <- npara[1]
  phi <- (N-npara)^(-1)*sum(sum_p)
  
  if(is.null(M_study)){
    alpha <- NULL  
  }else{
    M <- (phi*sum(M_study))-(phi*npara) 
    if(length(alphaM) > 1){
      dims <- c()
      for(i in 1:length(alphaM)){
        dims <- append(dims, length(alphaM[[i]]))
      }
      sum.alpha <- rep(0, max(dims))
      for(i in 1:max(dims)){
        for(s in 1:length(alphaM)){
          sum.alpha[i] <- sum(sum.alpha[i],alphaM[[s]][i], na.rm=T)
        }
      }
    }else{
      sum.alpha <- unlist(alphaM)
    }
    alpha <- M^(-1)*sum.alpha
  }  
  
  # ouput
  list(alpha=alpha, phi=phi)
}
