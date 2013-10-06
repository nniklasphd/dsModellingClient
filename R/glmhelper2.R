#'
#' @title Extracts object variables from a regression formula
#' @details This is an internal function required by the client
#' function \code{ds.glm} to extract variables from a regression formula
#' @param input a regression formula of the form 'y~x2+x2...' or 'y~x1*x2+x3...'
#' @return a list which contains the individual objects of the linear predictor 
#' but not the 'intercept' if there is any.
#' @author Gaye, A.
#' 
glmhelper2 <- function(input){
  outvar <- terms(input)[[2]]
  explvars <- terms(input)[[3]]
  tempholder <- outvar
  all.terms <- unlist(as.list(explvars))
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
    all.terms <- all.terms[[2]]
  }
  # this is the last element in the formula
  # it can potentially be the 'intercept' value
  # if that it is the case do not add it to
  # the list of variables to check
  if(!(is.numeric(all.terms))){
    tempholder <- append(tempholder, all.terms)
  }
  # if the last term is the just the name of the assigned table do not add it
  if(as.character(all.terms)[[1]] == "$"){
    symbol <- strsplit(deparse(all.terms), "\\$", perl=TRUE)[[1]][1]
    if(!(as.character(all.terms)[[1]] == symbol)){
      tempholder <- append(tempholder, all.terms)
    }
  }
  # if the last term is the name of a variable created in the workspace, add it
  if(class(output) == "name"){
    tempholder <- append(tempholder, all.terms)
  }

  return(tempholder)
}