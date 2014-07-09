#'
#' @title Extracts the elements of an expression of type call
#' @description This is an internal function required by the client
#' function \code{ds.glm} to extract variables from a regression formula
#' @param input a call of the form  'x1+x2+...' or 'x1*x2*...', typically 
#' this is the rightern side of a regression formula (i.e. the linear predictor) 
#' with the intercept.
#' @keywords internal
#' @return a list which contains the individual elements of the input expression
#' @author Gaye, A.
#' 
glmhelper1 <- function(input){
  holder <- c()
  if(length(input) == 1 & class(input)=="name"){
    holder <- input
    return(holder)
  }else{
    while(length(input) == 3){
      if(as.character(input)[[1]] == "$"){
        covar <- input
        holder <- append(holder, covar)
        break 
      }else{
        covar <- input[[3]]
        holder <- append(holder, covar)
        input <- input[[2]] 
      }
    }
  }
  return(holder)
}