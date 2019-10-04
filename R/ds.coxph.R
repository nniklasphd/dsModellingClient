#' 
#' @title Runs the distributed Cox model learning algorithm. 
#' @description A function that fits a distributed Cox proportional hazards regression model.
#' @details It enables a distributed approach of the Cox proportional model. The advantage is that (sensitive) data
#' doesn't have to be uploaded to a central location for the calculations to be performed.
#' @param data a character, the name of the table that holds the data.
#' @param survival_time a character, the name of the survival time column.
#' @param survival_event a character, the name of a the survival event column.
#' @param terms a comma separated string containing the model features.
#' @param maxit numeric, max number of algorithm iterations.
#' @param datasources a list of opal object(s) obtained after login in to opal servers;
#' these objects hold also the data assign to R, as \code{data frame}, from opal datasources.
#' @return a named vector of feature weights
#' @author Niklas, N.
#' @export
#' @examples {
#' 
#'   # load that contains the login details
#'   data(logindata)
#' 
#'   # login
#'   opals <- datashield.login(logins=logindata,assign=TRUE)
#' 
#'   # Example 1: run a Cox proportional-hazards model
#'   beta <- ds.coxph(data = "D", survival_time = "TIME", survival_event = "CENSOR", 
#'                    terms = "AAE,BDS,HU,CU,IVDUPN,IVDURN,NPDT,RACE,TREAT,SITE", maxit = 500)
#' 
#'   # clear the Datashield R sessions and logout
#'   datashield.logout(opals)
#' 
#' }
#' 
ds.coxph = function(data = NULL, survival_time = NULL, survival_event = NULL, terms = NULL, maxit = 600, datasources = NULL){
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
  
  # Initialization step1
  cally         <- call('coxphDS1', data, survival_time)
  study.summary <- datashield.aggregate(datasources, cally, async = TRUE)
  data_times    <- sort(unique(unlist(opal:::.select(study.summary, 'time.values'))))
   
  # Initialization step2
  data_times_str <- paste0(as.character(data_times), collapse=",")
  cally2         <- call('coxphDS2', data, survival_time, survival_event, terms, data_times_str)
  study.summary  <- datashield.aggregate(datasources, cally2, async = TRUE)
  
  # Calculate study index, DI and sumZ
  DI    <- Reduce(f="+", opal:::.select(study.summary, 'DI'))
  sumZ  <- Reduce(f="+", opal:::.select(study.summary, 'sum.Z'))
  
  n_features <- length(unlist(strsplit(terms,split=",")))
  beta1           <- matrix(rep(0, len=n_features))
  converge.state  <- FALSE
  epsilon         <- 1E-6
  iteration.count <- 0
	
  # Calculate weights until convergence or until max iterations
  while(!converge.state && iteration.count < maxit) {
    iteration.count <- iteration.count + 1
    beta0           <- beta1;
    beta0_str       <- paste0(as.character(beta0), collapse=",")
    cally3          <- call('coxphDS3', data, survival_time, survival_event, terms, beta0_str, data_times_str)
    study.summary   <- datashield.aggregate(datasources, cally3, async = TRUE)
    
    ebz  <- Reduce(f="+", opal:::.select(study.summary, 'ebz'))
    zebz  <- Reduce(f="+", opal:::.select(study.summary, 'zebz'))
    zzebz  <- Reduce(f="+", opal:::.select(study.summary, 'zzebz'))
	
	
    gradient <- matrix(colSums(sumZ)-colSums(zebz/ebz*DI),n_features,1)
	
	for(i in 1:n_features)
	{
		for(j in 1:n_features)
		{
			zzebz[,i,j] <- zzebz[,i,j]/ebz - zebz[,i]*zebz[,j]/ebz^2
		}
	}
	
	neghessian <- apply(zzebz*DI,c(2,3),sum)
	beta1 <- beta0 + solve(neghessian + diag(10^(-6),n_features)) %*% gradient
    
    
    converge.state <- (sum(abs(beta0 - beta1)) <= epsilon)
  }
  if (!converge.state) {
      warning(paste("Did not converge after", maxit, "iterations. Increase maxit parameter as necessary."))
  }
  se <- sqrt(diag(solve(neghessian)))
  zvalue <- beta1/se
  pvalue <- rep(0,length(zvalue))
  for(i in 1:length(pvalue)){pvalue[i] <- 2*pnorm(zvalue[i],lower.tail=(zvalue[i]<0))}
  #pvalue <- 2*pnorm(zvalue)
  return(list(beta1,se,zvalue,pvalue))
}
