#' 
#' @title Deals with differing factor levels between studies
#' @description This is an internal function required by the client
#' function \code{ds.glm} to carry ensure all variables ahve the same number of levels
#' in all studies.
#' @details The function goes through the variables in a formula and looks for
#' those of factor type. The number of levels are then compared across all the studies.
#' If the number of levels differ, 'dummy' levels are added to the studies with a lower 
#' number of levels and the formula is amended to include the recoded variable(s).
#' @param form an object of type formula, the input formula
#' @param variables a list of variable names given as calls
#' @param datasources datasources a list of opal object(s) obtained after login to opal servers;
#' these objects also hold the data assigned to R, as a \code{dataframe}, from opal datasources
#' @keywords internal
#' @return an object of type formula, an amended formula if number of levels differ across studies
#' for any one or more variables and the initial formula if levels were the same across studies.
#'
glmhelper4 <- function(form, variables, datasources){
  temp <- variables
  formula <- form
  
  for(i in 1:length(temp)){
    # check if variable is categorical
    cally0 <- call("class", temp[[i]])
    aa <- datashield.aggregate(datasources, cally0)
    
    # get the name of the variable
    initvarname <- deparse(temp[[i]])
    qq <- unlist(strsplit(deparse(temp[[i]]), "\\$", perl=TRUE))
    if(length(qq) > 1){ varName <- qq[2] }else{ varName <- qq[1] }
    
    if(length(unique(unlist(aa))) == 1 & unique(unlist(aa)) == "factor"){
      # if variable is categorical check if number of levels are the same across studies
      records <- c()
      xx <- list()
      for(j in 1:length(datasources)){
        cally <- call("levels", temp[[i]])
        bb <- datashield.aggregate(datasources[j], cally)
        xx[[j]] <- bb[[1]]
        records <- append(records, length(bb[[1]]))
      }
      # if number of levels are different across studies identify studies 
      # that require additional empty levels and generate those levels
      if(length(unique(records)) > 1){
        indicator <- 1
        studies2level <- which(records != max(records, na.rm=TRUE))
        levels2generate <- which(records == max(records, na.rm=TRUE))[1]
        cally <- paste0("recodeLevelsDS(", temp[[i]], ", c(","'",paste(xx[[levels2generate]],collapse="','"),"')",")")
        datashield.assign(datasources[studies2level], varName, as.symbol(cally))
        
        # because in datasshield the formula pass on to the server site has to be same across all studies
        # we need to generate a variable with the same name as the one recoded above for all studies
        # we do this by just assigning the variable in the studies that alread have the 'right' number of levels
        datashield.assign(datasources[-studies2level], varName, temp[[i]])
        
        # amend the formula to use the recoded variable
        amendedform <- gsub(initvarname, varName, deparse(formula), fixed=TRUE)
        formula <- as.formula(amendedform)
        
        # turn the input variable into numeric (required for the glm process)
        cally1 <- paste0('as.numeric(as.character(', varName, '))' )
        datashield.assign(datasources[studies2level], varName, as.symbol(cally1))
        cally2 <- call('as.numeric',call('as.character', temp[[i]]))
        datashield.assign(datasources[-studies2level], varName, cally2)
      }else{
        # turn the input variable into numeric (required for the glm process)
        cally3 <- call('as.numeric',call('as.character', temp[[i]]))
        datashield.assign(datasources, varName, cally3)      
      }
    }else{
      # if the variable is not a factor in all the studies stop and throw an error message
      if(length(unique(unlist(aa))) > 1){
        stop(paste0("The variable ", varName, " is not a factor in all the studies"), call.=FALSE)
      }
      # turn the input variable into numeric (required for the glm process)
      cally4 <- call('as.numeric',call('as.character', temp[[i]]))
      datashield.assign(datasources, varName, cally4)
    }
    
  }

  return(formula)
}



