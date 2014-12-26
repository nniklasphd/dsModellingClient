#'
#' @title Checks if the elements in a regression formula are defined and not empty
#' @description This is an internal function required by the client
#' function \code{ds.glm} to ensure all the variable in the LP are defined and not empty,
#' i.e. are not missing at complete.
#' @param formula a character, a regression formula given as a string character
#' @param data a character, the name of an optional data frame containing the variables in 
#' in the \code{formula}.
#' @param datasources a list of opal object(s) obtained after login in to opal servers;
#' these objects hold also the data assign to R, as \code{dataframe}, from opal datasources.
#' @keywords internal
#' @return an integer 0 if check was passed and 1 if failed
#' @author Gaye, A.
#' 
glmhelper1 <- function(formula, data, datasources){
  
  # turn the formula into a character
  formula <- paste0(Reduce(paste, deparse(formula)))
  
  # replace the symbols '~', '+' and '*' by a separator
  formula <- gsub( " ", "", formula, fixed=TRUE)
  formula <- gsub( "~", "|", formula, fixed=TRUE)
  formula <- gsub( "+", "|", formula, fixed=TRUE)
  formula <- gsub( "*", "|", formula, fixed=TRUE)
  
  # split the input formula by "|" to obtain the names of the variables
  elts <- unlist(strsplit(formula, split="|", fixed=TRUE))
  
  # check that each variable is defined and not empty and each study. Stop the process if any check fails
  stdnames <- names(datasources)
  for(i in 1:length(elts)){
    if(is.na(as.numeric(elts[i], options(warn=-1)))){ # making sure an eventual intercept term is not included in the checks
      message(paste0("    ", elts[i], "..."))
      for(j in 1: length(datasources)){
        
        # check if the variable is defined on the server site
        myterms <- unlist(strsplit(elts[i], split='$', fixed=TRUE))
        if(length(myterms) > 1){
          cally <- call("exists", myterms[1])
          out <- datashield.aggregate(datasources[j], cally)
          if(!(out[[1]])){ 
            stop(paste0("'", myterms[1], "' is not defined in ", stdnames[j], "!"), call.=FALSE)
          }else{
            cally <- paste0("colnames(", myterms[1], ")")
            clnames <- unlist(datashield.aggregate(datasources[j], as.symbol(cally)))
            if(!(myterms[2] %in% clnames)){
              stop(paste0("'", myterms[2], "' is not defined in ", stdnames[j], "!"), call.=FALSE)
            }
          }         
        }else{
          if(!(is.null(data))){
            defined <- isDefined(datasources, paste0(data, "$", elts[i]))   
            call <- paste0("isNaDS(", paste0(data, "$", elts[i]), ")")
          }else{
            defined <- isDefined(datasources, elts[i]) 
            call <- paste0("isNaDS(", elts[i], ")")
          }
        }
        # check if variable is not missing at complete
        out1 <- datashield.aggregate(datasources[j], as.symbol(call))
        if(out1[[1]]){ 
          stop("The variable ", elts[i], " in ", stdnames[j], " is missing at complete (all values are 'NA').", call.=FALSE)
        }
      }
    }
  }
}