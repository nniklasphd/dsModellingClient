#' 
#' @title Runs the distributed Cox model learning algorithm. 
#' @description A function that fits a distributed Cox proportional hazards regression model.
#' @details It enables a distributed approach of the Cox proportional model. The advantage is that (sensitive) data
#' doesn't have to be uploaded to a central location for the calculations to be performed.
#' @param survival_time a character, the name of a numerical vector
#' @param survival_event a character, the name of a numerical vector
#' @param terms a character, the name of a numerical vector
#' @param data a character, the name of an optional data frame containing the variables in 
#' @param datasources a list of opal object(s) obtained after login in to opal servers;
#' these objects hold also the data assign to R, as \code{data frame}, from opal datasources.
#' @return a named vector of feature weights
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
#'   # Example 1: run a Cox proportional-hazards model
#'   beta <- ds.coxph(data = "D", survival_time = "D$TIME", survival_event = "D$CENSOR", 
#'                    terms = "D$AAE,D$BDS,D$HU,D$CU,D$IVDUPN,D$IVDURN,D$NPDT,D$RACE,D$TREAT,D$SITE", maxit = 600)
#' 
#'   # clear the Datashield R sessions and logout
#'   datashield.logout(opals)
#' 
#' }
#' 
ds.coxph = function(survival_time = NULL, survival_event = NULL, terms = NULL, data = NULL, maxit = 600, datasources = NULL){
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
  
  # Initialization step1
  cally         <- call('coxphDS1', survival_time, survival_event, terms, data)
  study.summary <- datashield.aggregate(datasources, cally, async = TRUE)
  data_zzc      <- Reduce(f="+", as.vector(opal:::.select(study.summary, 'ZZvc')))
  study_length  <- lapply(opal:::.select(study.summary, 'time.values'), length)
  data_times    <- unique(unlist(opal:::.select(study.summary, 'time.values')))
  inv_ZZc       <- solve(data_zzc)
  
  # Initialization step2
  data_times_str <- paste0(as.character(data_times), collapse=",")
  cally2         <- call('coxphDS2', survival_time, survival_event, terms, data_times_str, data)
  study.summary  <- datashield.aggregate(datasources, cally2, async = TRUE)
  
  # Calculate study index, DI and sumZ
  index <- Reduce(f="+", opal:::.select(study.summary, 'index'))
  DI    <- Reduce(f="+", opal:::.select(study.summary, 'DI'))
  sumZ  <- Reduce(f="+", opal:::.select(study.summary, 'sum.Z'))
  study_index <- study_DI <- study_sumZ <- list()
  index <- cumsum(c(0, index[1:(length(index)-1)])) + 1
  for (s in 1:numstudies) {
    if (s == 1) {
      idx = (index <= study_length[[s]]) 
      study_index[[s]] <- index[idx]
    } else {
      idx = (index > study_length[[s-1]])
      study_index[[s]] = index[idx] - study_length[[s-1]] 
    }
    study_DI[[s]]   <- DI[idx]
    study_sumZ[[s]] <- sumZ[idx,]
  }

  n_features      <- ncol(data_zzc)
  beta1           <- matrix(rep(0, len=n_features))
  converge.state  <- FALSE
  epsilon         <- 1E-6
  iteration.count <- 0
  # Calculate weights until convergence or until max iterations
  while(!converge.state && iteration.count < maxit) {
    iteration.count <- iteration.count + 1
    beta0           <- beta1;
    beta0_str       <- paste0(as.character(beta0), collapse=",")
    cally3          <- call('coxphDS3', survival_time, survival_event, terms, beta0_str, data)
    study.summary   <- datashield.aggregate(datasources, cally3, async = TRUE)
    
    thetac_addition     <- 0
    thetaZtmpc_addition <- 0
    ZBc <- thetaZtmpc <- thetac <- thetaZc <- Gvc <- list()
    ## calculate Gvc
    for (s in numstudies:1) {
      ZBc[[s]]        <- study.summary[[s]]$exp.Zc.beta
      thetaZtmpc[[s]] <- study.summary[[s]]$theta.Ztmpc
      thetaZtmpc[[s]][nrow(thetaZtmpc[[s]]),] <- thetaZtmpc_addition + thetaZtmpc[[s]][nrow(thetaZtmpc[[s]]),]
      thetaZtmpc[[s]]     <- apply(apply(apply(thetaZtmpc[[s]], 2, rev), 2, cumsum), 2, rev)
      thetaZtmpc_addition <- thetaZtmpc[[s]][1,]
      
      thetac[[s]]     <- rev(t(thetac_addition + apply(apply(ZBc[[s]], 2, rev), 2, cumsum)))
      thetac_addition <- thetac[[s]][[1]]
      thetaZtmpc[[s]] <- thetaZtmpc[[s]] / do.call("cbind", rep(list(thetac[[s]]), n_features))
      thetaZc[[s]]    <- thetaZtmpc[[s]][study_index[[s]],]
      thetaZc[[s]]    <- thetaZc[[s]] * do.call("cbind", rep(list(study_DI[[s]]), n_features))
      Gvc[[s]]        <- study_sumZ[[s]] - thetaZc[[s]]
    }
    G              <- Reduce(f = "+", lapply(Gvc, colSums))
    beta1          <- beta0 + inv_ZZc %*% as.vector(Conj(t.default(G)))
    converge.state <- (sum(abs(beta0 - beta1)) <= epsilon)
  }
  if (!converge.state) {
      warning(paste("Did not converge after", maxit, "iterations. Increase maxit parameter as necessary."))
      return(NULL)
  }
  return(beta1)
}
