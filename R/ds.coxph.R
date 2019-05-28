#' 
#' @title Runs a cox proportional-hazards model.
#' @description A function that fits a Cox proportional hazards regression model.
#' @details It is a wrapper for the server side function.
#' @param survival_time a character, the name of a numerical vector
#' @param survival_event a character, the name of a numerical vector
#' @param terms a character, the name of a numerical vector
#' @param method method for tie handling (e.g. "breslow" or "efron")
#' @param data a character, the name of an optional data frame containing the variables in 
#' @param datasources a list of opal object(s) obtained after login in to opal servers;
#' these objects hold also the data assign to R, as \code{data frame}, from opal datasources.
#' @return coefficients a named vector of coefficients
#' @return residuals the 'working' residuals, that is the residuals in the final
#' iteration of the IWLS fit.
#' @return method the \code{method} object used.
#' @return linear.predictors the linear fit on link scale.
#' @author Inberg, G
#' @export
#' @examples {
#' 
#'   # load that contains the login details
#'   data(logindata)
#' 
#'   # login
#'   opals <- datashield.login(logins=logindata,assign=TRUE)
#' 
#'   # Example 1: run a Cox proportional-hazards model (e.g. survival vs race)
#'   mod <- ds.coxph("D$time", "D$censor", "D$race", "breslow")
#'   mod
#' 
#'   # clear the Datashield R sessions and logout
#'   datashield.logout(opals)
#' 
#' }
#' 
ds.coxph = function(survival_time = NULL, survival_event = NULL, terms = NULL, method = NULL, data = NULL, maxit = 600, datasources = NULL){
  result <- NULL
  
  # if no opal login details are provided look for 'opal' objects in the environment
  if (is.null(datasources)){
    datasources <- findLoginObjects()
  }
  
  if(!(is.null(data))){
    if (!isDefined(datasources, data)) {
      stop("Please provide a valid dataset", call.=FALSE)
    }
  }
  
  if(is.null(survival_time)){
    stop("Please provide the survival time column.", call.=FALSE)
  }
  
  if(is.null(survival_event)){
    stop("Please provide the survival event column.", call.=FALSE)
  }
  
  if(is.null(terms)){
    stop("Please provide the terms columns", call.=FALSE)
  }
  
  # iteration counter
  iteration.count <- 0
  
  # number of 'valid' studies (those that passed the checks) and vector of beta values
  numstudies <- length(datasources)
  
  # Initialization
  cally    <- call('coxphDS1', survival_time, survival_event, terms, method, data)
  data_sum <- datashield.aggregate(datasources, cally, async = TRUE)
  data_zzc <- Reduce(f="+", data_sum)
  inv_ZZc  <- solve(data_zzc)
   
  # Algorithm
  n_features      <- ncol(data_zzc)
  beta1           <- matrix(rep(0, len=n_features))
  converge.state  <- FALSE
  epsilon         <- 1E-6
  iteration.count <- 0
  while(!converge.state && iteration.count < maxit) {
    iteration.count <- iteration.count + 1
    print(paste("Iteration:", iteration.count))
    
    beta0 <- beta1;
    #IN LEGAL TRANSMISSION FORMAT ("0,0,0,0,0")
    beta0_str <- paste0(as.character(beta0), collapse=",")

    #NOW CALL SECOND COMPONENT OF coxphDS
    cally2 <- call('coxphDS2',  survival_time, survival_event, terms, method, beta0_str, data)
    study.summary <- datashield.aggregate(datasources, cally2, async = TRUE)

    G <- Reduce(f = "+", study.summary)
    beta1 = beta0 + inv_ZZc %*% as.vector(Conj(t.default(G)))
    converge.state <- (sum(abs(beta0 - beta1)) <= epsilon)
  }

  if (!converge.state) {
      warning(paste("Did not converge after", maxit, "iterations. Increase maxit parameter as necessary."))
      return(NULL)
  }
  result        <- as.vector(beta1)
  names(result) <- rownames(beta1)
  return(result)
}
