#' 
#' @title Runs a combined GLM analysis of non-pooled data
#' @param datasources a list of opal object(s) obtained after login to opal servers;
#' these objects also hold the data assigned to R, as a \code{dataframe}, from opal datasources.
#' @param formula an object of class \code{formula} which describes the model to be fitted
#' @param family a description of the error distribution function to use in the model
#' @param maxit the number of iterations of IWLS used
#' @details It enables a parallelized analysis of individual-level data sitting 
#' on distinct servers by sending 
#' instructions to each computer requesting non-disclosing summary statistics.
#' The summaries are then combined to estimate the parameters of the model; these
#' parameters are the same as those obtained if the data were 'physically' pooled.
#' @return coefficients a named vector of coefficients
#' @return residuals the 'working' residuals, that is the residuals in the final
#' iteration of the IWLS fit.
#' @return fitted.values the fitted mean values, obtained by transforming the
#' linear predictors by the inverse of the link function.
#' @return rank the numeric rank of the fitted linear model.
#' @return family the \code{family} object used.
#' @return linear.predictors the linear fit on link scale.
#' @return aic A version of Akaike's An Information Criterion, which tells how 
#' well the model fits
#' @authors Burton, P. (p.burton@bristol.ac.uk);  Gaye, A. (amadou.gaye@bristol.ac.uk); Laflamme, P.
#' @export
#' @examples {
#' 
#' # load the file that contains the login details
#' data(logindata)
#' 
#' # login and assign some variables to R
#' myvar <- list("DIS_DIAB","PM_BMI_CONTINUOUS","LAB_HDL", "GENDER")
#' opals <- datashield.login(logins=logindata,assign=TRUE,variables=myvar)
#' 
#' # Example 1: run a GLM without interaction (e.g. diabetes prediction using BMI and HDL levels and GENDER)
#' mod <- ds.glm(datasources=opals,formula=D$DIS_DIAB~D$PM_BMI_CONTINUOUS+D$LAB_HDL+D$GENDER,family=quote(binomial))
#'  
#' # Example 2: run the above GLM model with an intercept (eg. intercept = 1)
#'  mod <- ds.glm(datasources=opals,formula=D$DIS_DIAB~1+D$PM_BMI_CONTINUOUS+D$LAB_HDL+D$GENDER,family=quote(binomial))
#'  
#' # Example 3: run the above GLM with interaction HDL and GENDER
#'  mod <- ds.glm(datasources=opals,formula=D$DIS_DIAB~D$PM_BMI_CONTINUOUS+D$LAB_HDL*D$GENDER,family=quote(binomial))
#'  
#' # Example 4: now run the same GLM but with interaction between BMI and HDL 
#'  mod <- ds.glm(datasources=opals,formula=D$DIS_DIAB~D$PM_BMI_CONTINUOUS*D$LAB_HDL+D$GENDER,family=quote(binomial))
#' }
#'
ds.glm <- function(datasources=NULL, formula=NULL, family=NULL, maxit=15) {
  
  if(is.null(datasources)){
    cat("\n\n ALERT!\n")
    cat(" No valid opal object(s) provided.\n")
    cat(" Make sure you are logged in to valid opal server(s).\n")
    stop(" End of process!\n\n", call.=FALSE)
  }
  
  if(is.null(formula)){
    cat("\n\n ALERT!\n")
    cat(" Please provide a valid regression formula\n")
    stop(" End of process!\n\n", call.=FALSE)
  }
  
  if(is.null(family)){
    cat("\n\n ALERT!\n")
    cat(" Please provide a valid 'family' argument\n")
    stop(" End of process!\n\n", call.=FALSE)
  }
  
  # call the helper function that extracts the outcome and covariates 
  # from the regression formula as objects; these objects are required 
  # by the function that carries out the preliminary checks
  vars2check <-  dsmodellingclient:::glmhelper2(formula)
  
  # call the function that checks the variables are available and not empty
  datasources <- ds.checkvar(datasources, vars2check)
  
  # number of 'valid' studies (those that passed the checks) and vector of beta values
  numstudies <- length(datasources)
  beta.vect.next <- NULL
  
  #Iterations need to be counted. Start off with the count at 0
  #and increment by 1 at each new iteration
  iteration.count<-0
  
  #Provide arbitrary starting value for deviance to enable subsequent calculation of the
  #change in deviance between iterations
  dev.old<-9.99e+99
  
  #Convergence state needs to be monitored.
  converge.state<-FALSE
  
  #Define a convergence criterion. This value of epsilon corresponds to that used
  #by default for GLMs in R (see section S3 for details)
  epsilon<-1.0e-08
  
  f<-NULL
  
  while(!converge.state && iteration.count < maxit) {
    
    iteration.count<-iteration.count+1
    
    cat("--------------------------------------------\n")
    cat("Iteration", iteration.count, "\n")
    if(is.null(beta.vect.next)){
       beta.vect.temp <- NULL
     }else{
       beta.vect.temp <- paste0(beta.vect.next, collapse=",")
     }
    cally <- as.call(list(quote(glm.ds), formula, family, beta.vect.temp))
    
    study.summary <- datashield.aggregate(datasources, cally);
    
    .select <- function(l, field) {
     lapply(l, function(obj) {obj[[field]]})
    }
    
    info.matrix.total<-Reduce(f="+", .select(study.summary, 'info.matrix'))
    score.vect.total<-Reduce(f="+", .select(study.summary, 'score.vect'))
    dev.total<-Reduce(f="+", .select(study.summary, 'dev'))
    
    if(iteration.count==1) {
      # Sum participants only during first iteration.
      nsubs.total<-Reduce(f="+", .select(study.summary, 'numsubs'))
      # Save family
      f <- study.summary[[1]]$family
    }
    
    #Create variance covariance matrix as inverse of information matrix
    variance.covariance.matrix.total<-solve(info.matrix.total)
    
    # Create beta vector update terms
    beta.update.vect<-variance.covariance.matrix.total %*% score.vect.total
  
    #Add update terms to current beta vector to obtain new beta vector for next iteration
    if(is.null(beta.vect.next)) {
      beta.vect.next<-beta.update.vect
    } else {
      beta.vect.next<-beta.vect.next+beta.update.vect
    }
    
    #Calculate value of convergence statistic and test whether meets convergence criterion
    converge.value<-abs(dev.total-dev.old)/(abs(dev.total)+0.1)
    if(converge.value<=epsilon)converge.state<-TRUE
    if(converge.value>epsilon)dev.old<-dev.total
    
    #For ALL iterations summarise model state after current iteration
    cat("\nSUMMARY OF MODEL STATE after iteration No",iteration.count,
        "\n\nCurrent deviance",dev.total,"on",
        (nsubs.total-length(beta.vect.next)), "degrees of freedom",
        "\nConvergence criterion    ",converge.state," (", converge.value,")\n\n")
    
    cat("beta\n")
    print(as.vector(beta.vect.next))
    
    cat("Information matrix overall\n")
    print(info.matrix.total)
    
    cat("Score vector overall\n")
    print(score.vect.total)
    
    cat("Current Deviance\n")
    print(dev.total)
    cat("--------------------------------------------\n")    
  }
  
  #If convergence has been obtained, declare final (maximum likelihood) beta vector,
  #and calculate the corresponding standard errors, z scores and p values
  #(the latter two to be consistent with the output of a standard GLM analysis)
  #Then print out final model summary
  if(converge.state)
  {
    beta.vect.final<-beta.vect.next
    
    scale.par <- 1
    if(f$family== 'gaussian') {
      scale.par <- dev.total / (nsubs.total-length(beta.vect.next))
    }
    
    se.vect.final <- sqrt(diag(variance.covariance.matrix.total)) * sqrt(scale.par)
    z.vect.final<-beta.vect.final/se.vect.final
    pval.vect.final<-2*pnorm(-abs(z.vect.final))
    parameter.names<-names(score.vect.total[,1])
    model.parameters<-cbind(beta.vect.final,se.vect.final,z.vect.final,pval.vect.final)
    dimnames(model.parameters)<-list(parameter.names,c("Estimate","Std. Error","z-value","p-value"))
    
    glmds <- list(
      formula=formula,
      coefficients=model.parameters,
      dev=dev.total,
      nsubs=nsubs.total,
      df=(nsubs.total-length(beta.vect.next)),
      iter=iteration.count
    )
    
    class(glmds) <- 'glmds'
    
    glmds   
  } else {
    warning(paste("Did not converge after", maxit, "iterations. Increase maxit parameter as necessary."))
    NULL
  }
}
