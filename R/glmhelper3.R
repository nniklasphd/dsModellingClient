#'
#' @title Extracts the elements of an expression of type call
#' @description This is an internal function required by the client
#' function \code{ds.glm} to carry checks.
#' @details the functions checks if all the variables in the lp formula exist 
#' (i.e. are defined) on the server site and if any is empty (contains NAs only).
#' If any of the checks fails the process stops and an error message is thrown.
#' @param dtname the name if any of the data frame that hold the variables in the regression 
#' formula.
#' @param form an object of type formula which describes the model to be fitted.
#' @param datasources a list of opal object(s) obtained after login to opal servers;
#' these objects also hold the data assigned to R, as a \code{dataframe}, from opal 
#' datasources.
#' @keywords internal
#' @return if no error message is thrown a list with the names of the variables as 
#' call objects is returned.
#' @author Gaye, A.
#' 
glmhelper3 <- function(dtname, form, datasources){
  stdnames <- names(datasources)
  temp <- glmhelper2(form)
  variables <- c()
  for(i in 1:length(temp)){
    variables <- append(variables, temp[[i]])
  }
  for(i in 1:length(variables)){
    for(j in 1: length(datasources)){
      inputterms <- unlist(strsplit(deparse(variables[[i]]), "\\$", perl=TRUE))
      if(length(inputterms) < 2){
        cally1 <- call('exists', as.character(variables[[i]]))
        d1 <- datashield.aggregate(datasources[j], cally1)[[1]]
        if(!d1){
          stop("The variable ", as.character(variables[[i]]),  " is not defined in ", stdnames[j], ".", call.=FALSE)
        }else{
          cally2 <- call('isNaDS', variables[[i]])
          d2 <- datashield.aggregate(datasources[j], cally2)[[1]]
          if(d2){ 
            stop("The variable ", as.character(variables[[i]]), " in ", stdnames[j], " is empty (all values are 'NA').", call.=FALSE)
          }
        }
        # turn the vector into numeric
        #cally <- paste0('as.numeric(as.character(', variables[[i]], '))' )
        #datashield.assign(datasources[j], as.character(variables[[i]]), as.symbol(cally))
      }else{
        if(is.null(dtname)){
          stop("dtname=NULL! You must set the argument dtname to ", inputterms[1], " because the variable ", inputterms[2], " is indicated as attached ", inputterms[1], " in ", as.character(variables[[i]]), "!", call.=FALSE)
        }
        cally2 <- call('isNaDS', variables[[i]])
        d2 <- datashield.aggregate(datasources[j], cally2)[[1]]
        if(d2){ 
          stop("The variable ", as.character(variables[[i]]), " in ", stdnames[j], " is empty (all values are 'NA').", call.=FALSE)
        }
        # turn the vector into numeric
        #cally <- call("as.numeric",call("as.character", variables[[i]]))
        #datashield.assign(datasources[j], inputterms[2], cally)
      }
    }
    
  
  }
  return (temp)
}