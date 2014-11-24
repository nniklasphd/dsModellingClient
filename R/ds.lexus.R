#' 
#' @title Generates an expanded version of a dataset that contains survival data
#' @description This function is meant to be used as part of a piecewise regression analysis.
#' @details It splits the survial interval time of subjects into sub-intervals and reports the failure 
#' status of the subjects at each sub-interval. Each of those sub-interval is given an id e.g. if the overall
#' interval of a subject is split into 4 sub-interval, those sub-intervals have ids 1, 2, 3 and 4; so this is 
#' basically the count of periods for each subject. The interval ids are held in a column named "TIMEID". 
#' The entry and exit times in the input table are used to compute the total survival time. 
#' By default all the covariates in the input table are included in the expanded output table but it is 
#' preferable to indicate the names of the covariates to be included via the argument 'variables'.
#' @param x a character, the name of the table that holds the original data, this is the data to be expanded.
#' @param intervalWidth, a numeric vector which gives the chosen width of the intervals ('pieces'). 
#' This can be one value (in which case all the intervals have same width) or several different values.
#' If no value(s) are provided a single default value is used. That default value is the set to be the 
#' 1/10th of the mean of the exit time values across all the studies.
#' @param idCol a character the name of the column that holds the individual IDs of the subjects.
#' @param entryCol a character, the name of the column that holds the entry times (i.e. start of follow up).
#' If no name is provided the default is to set all the entry times to 0 in a column named "STARTTIME".
#' A message is then printed to alert the user as this has serious consequences if the actual entry times are 
#' not 0 for all the subjects. 
#' @param exitCol a character, the name of the column that holds the exit times (i.e. end of follow up).
#' @param statusCol a character, the name of the column that holds the 'failure' status of each subject, 
#' tells whether or not a subject has been censored.
#' @param variables a character vector, the column names of the variables (covariates) to include in the 
#' final expanded table. The input table might have a large number of covariates and if only some of those
#' variables are relevant for the sought analysis it make sense to only include those. By default (i.e. if
#' no variables are indicated) all the covariates in the inout table are included and this will lengthen the
#' run time of the function. 
#' @param newobj the name of the output expanded table. By default the name is the name of the input table with 
#' the suffixe "_expanded".
#' @param datasources a list of opal object(s) obtained after login to opal servers;
#' these objects also hold the data assigned to R, as a \code{data frame}, from opal datasources
#' @return a dataframe, an expanded version of the input table.
#' @author Gaye, A.; Burton, P.
#' 
ds.lexus <- function(x=NULL, intervalWidth=NULL, idCol=NULL, entryCol=NULL, exitCol=NULL, statusCol=NULL, variables=NULL, newobj=NULL, datasources=NULL){
  
  # if no opal login details were provided look for 'opal' objects in the environment
  if(is.null(datasources)){
    findLogin <- getOpals()
    if(findLogin$flag == 1){
      datasources <- findLogin$opals
    }else{
      if(findLogin$flag == 0){
        stop(" Are yout logged in to any server? Please provide a valid opal login object! ", call.=FALSE)
      }else{
        message(paste0("More than one list of opal login object were found: '", paste(findLogin$opals,collapse="', '"), "'!"))
        userInput <- readline("Please enter the name of the login object you want to use: ")
        datasources <- eval(parse(text=userInput))
        if(class(datasources[[1]]) != 'opal'){
          stop("End of process: you failed to enter a valid login object", call.=FALSE)
        }
      }
    }
  }
  
  # check if user have provided a the name of the input dataset
  if(is.null(x)){
    stop("Please provide the name of the dataset to expand!", call.=FALSE)
  }
  
  # check if user have provided the name of the column that holds the subject ids
  if(is.null(idCol)){
    stop("Please provide the name of the column that holds the subject IDs!", call.=FALSE)
  }
  
  # check if user have provided the name of the column that holds failure information
  if(is.null(statusCol)){
    stop("Please provide the name of the column that holds 'failure' information!", call.=FALSE)
  }
  
  # check if user have provided the name of the column that holds exit times 
  if(is.null(exitCol)){
    stop("Please provide the name of the column that holds the exit times (i.e. end of follow up time)!", call.=FALSE)
  }
  
  # if no value provided for 'intervalWidth) generate one
  if(is.null(intervalWidth)){
    intervalWidth <- lexusHelper1(datasources, paste0(x,"$",exitCol))
  }
  
  if(is.null(newobj)){
    newobj <- paste0(x,"_expanded")
  }
  
  # call the server side function
  cally <- call("lexusDS", x, intervalWidth, idCol, entryCol, exitCol, statusCol, variables)
  datashield.assign(datasources, newobj, cally)
  
  # check that the new object has been created if and display a message accordingly
  finalcheck <- isAssigned(datasources, newobj)  
  
}