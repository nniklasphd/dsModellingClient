#'
#' @title Extracts object variables from a regression formula
#' @description This is an internal function required by the client
#' function \code{ds.glm} to extract variables from a regression formula
#' @param input a regression formula of the form 'y~x2+x2...' or 'y~x1*x2+x3...'
#' @keywords internal
#' @return a list which contains the individual objects of the linear predictor 
#' but not the 'intercept' if there is one.
#' @author Gaye, A.
#' 
glmhelper2 <- function(input){
  outvar <- terms(input)[[2]]
  explvars <- terms(input)[[3]]
  tempholder <- outvar
  if(as.character(explvars)[[1]] == "$"){
    all.terms <- explvars
  }else{
    all.terms <- unlist(as.list(explvars))
  }
  if(length(explvars) == 1){
    tempholder <- append(tempholder, explvars)
  }else{
    while(length(all.terms) > 1){
      if(as.character(all.terms)[[1]] == "$"){
        last.term <- all.terms
      }else{
        last.term <- all.terms[[3]]
      }
      output <- glmhelper1(last.term)
      if(class(output) == "name"){
        tempholder <- append(tempholder, output)
      }else{
        for(j in 1:length(output)){
          tempholder <- append(tempholder, output[[j]])
        }
      }
      ff <- all.terms
      all.terms <- all.terms[[2]]
    }
    # this is the last element in the formula
    # it can potentially be the 'intercept' value
    # if that it is the case do not add it to
    # the list of variables to check
    if(!(is.numeric(all.terms))){
      if(class(all.terms) == "name"){
        if(!(as.character(ff)[[1]] == "$")){
          tempholder <- append(tempholder, all.terms)
        }
      }
    }
  }
  qq <- as.character(ff[[3]])
  if(qq[1] == '*' ){
    tempholder <- append(tempholder, qq[2])
  }
  return(tempholder)
}
