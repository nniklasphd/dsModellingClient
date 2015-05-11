#' 
#' @title Runs a combined GLM analysis of non-pooled data
#' @description A function fit generalized linear models
#' @details It enables a parallelized analysis of individual-level data sitting 
#' on distinct servers by sending 
#' @param formula a character, a formula which describes the model to be fitted
#' @param family a description of the error distribution function to use in the model
#' @param startBetas starting values for the parameters in the linear predictor
#' @param offset  a character, null or a numeric vector that can be used to specify an a priori known component 
#' to be included in the linear predictor during fitting.
#' @param weights  a character, the name of an optional vector of 'prior weights' to be used in the fitting 
#' process. Should be NULL or a numeric vector.
#' @param data a character, the name of an optional data frame containing the variables in 
#' in the \code{formula}. The process stops if a non existing data frame is indicated. 
#' @param checks a boolean, if TRUE (default) checks that takes 1-3min are carried out to verify that the 
#' variables in the model are defined (exist) on the server site and that they have the correct characteristics
#' required to fit a GLM. The default value is FALSE because checks lengthen the runtime and are mainly meant to be 
#' # used as help to look for causes of eventual errors.
#' @param maxit the number of iterations of IWLS used
#' instructions to each computer requesting non-disclosing summary statistics.
#' The summaries are then combined to estimate the parameters of the model; these
#' parameters are the same as those obtained if the data were 'physically' pooled.
#' @param CI a numeric, the confidence interval.
#' @param viewIter a boolean, tells whether the results of the intermediate iterations
#' should be printed on screen or not. Default is FALSE (i.e. only final results are shown).
#' @param datasources a list of opal object(s) obtained after login to opal servers;
#' these objects also hold the data assigned to R, as a \code{dataframe}, from opal datasources.
#' @return coefficients a named vector of coefficients
#' @return residuals the 'working' residuals, that is the residuals in the final
#' iteration of the IWLS fit.
#' @return fitted.values the fitted mean values, obtained by transforming the
#' linear predictors by the inverse of the link function.
#' @return rank the numeric rank of the fitted linear model.
#' @return family the \code{family} object used.
#' @return linear.predictors the linear fit on link scale.
#' 
#' @author Burton,P;Gaye,A;Laflamme,P
#' @seealso \link{ds.lexis} for survival analysis using piecewise exponential regression
#' @seealso \link{ds.gee} for generalized estimating equation models
#' @export
#' @examples {
#' 
#'  # load the file that contains the login details
#'  data(glmLoginData)
#' 
#'  # login and assign all the variables to R
#'  opals <- datashield.login(logins=glmLoginData, assign=TRUE)
#' 
#'  # Example 1: run a GLM without interaction (e.g. diabetes prediction using BMI and HDL levels and GENDER)
#'  mod <- ds.glm(formula='D$DIS_DIAB~D$GENDER+D$PM_BMI_CONTINUOUS+D$LAB_HDL', family='binomial')
#'  mod

#'  # Example 2: run the above GLM model without an intercept
#'  # (produces separate baseline estimates for Male and Female)
#'  mod <- ds.glm(formula='D$DIS_DIAB~0+D$GENDER+D$PM_BMI_CONTINUOUS+D$LAB_HDL', family='binomial')
#'  mod

#'  # Example 3: run the above GLM with interaction between GENDER and PM_BMI_CONTINUOUS
#'  mod <- ds.glm(formula='D$DIS_DIAB~D$GENDER*D$PM_BMI_CONTINUOUS+D$LAB_HDL', family='binomial')
#'  mod

#'  # Example 4: Fit a standard Gaussian linear model with an interaction
#'  mod <- ds.glm(formula='D$PM_BMI_CONTINUOUS~D$DIS_DIAB*D$GENDER+D$LAB_HDL', family='gaussian')
#'  mod

#'  # Example 5: now run a GLM where the error follows a poisson distribution
#'  # P.S: A poisson model requires a numeric vector as outcome so in this example we first convert
#'  # the categorical BMI, which is of type 'factor', into a numeric vector
#'  ds.asNumeric('D$PM_BMI_CATEGORICAL','BMI.123')
#'  mod <- ds.glm(formula='BMI.123~D$PM_BMI_CONTINUOUS+D$LAB_HDL+D$GENDER', family='poisson')
#'  mod
#'  
#'  # clear the Datashield R sessions and logout
#'  datashield.logout(opals) 
#' }
#'
ds.glm <- function(formula=NULL, data=NULL, family=NULL, offset=NULL, weights=NULL, checks=FALSE, maxit=15, CI=0.95, viewIter=FALSE, datasources=NULL) {
  
 # if no opal login details are provided look for 'opal' objects in the environment
  if(is.null(datasources)){
    datasources <- findLoginObjects()
  }
  
  # verify that 'formula' was set
  if(is.null(formula)){
    stop(" Please provide a valid regression formula!", call.=FALSE)
  }

# check if user gave offset or weights directly in formula, if so the argument 'offset' or 'weights'
# to provide name of offset or weights variable
    if(sum(as.numeric(grepl('offset', formula, ignore.case=TRUE)))>0 ||
       sum(as.numeric(grepl('weights', formula, ignore.case=TRUE)))>0)
	{
       cat("\n\n WARNING: you may have specified an offset or regression weights")
       cat("\n as part of the model formula. In ds.glm (unlike the usual glm in R)")
       cat("\n you must specify an offset or weights separately from the formula")
       cat("\n using the offset or weights argument.\n\n")
	}

  formula <- as.formula(formula)

  
  # check that 'family' was set
  if(is.null(family)){
    stop(" Please provide a valid 'family' argument!", call.=FALSE)
  }
  
  # if the argument 'data' is set, check that the data frame is defined (i.e. exists) on the server site
  if(!(is.null(data))){
    defined <- isDefined(datasources, data)
  }
  
  # beginning of optional checks - the process stops if any of these checks fails #
  if(checks){
    message(" -- Verifying the variables in the model")
    # call the function that checks the variables in the formula are defined (exist) on the server site and are not missing at complete
 
   glmChecks(formula, data, offset, weights, datasources)
  }else{
    #message("WARNING:'checks' is set to FALSE; variables in the model are not checked and error messages may not be intelligible!")
  }

#MOVE ITERATION COUNT BEFORE ASSIGNMENT OF beta.vect.next  
#Iterations need to be counted. Start off with the count at 0
  #and increment by 1 at each new iteration
  iteration.count<-0

  # number of 'valid' studies (those that passed the checks) and vector of beta values
  numstudies <- length(datasources)
 

#ARBITRARY LENGTH FOR START BETAs AT THIS STAGE BUT IN LEGAL TRANSMISSION FORMAT ("0,0,0,0,0")
 beta.vect.next <- rep(0,5)
 beta.vect.temp <- paste0(as.character(beta.vect.next), collapse=",")
  

#IDENTIFY THE CORRECT DIMENSION FOR START BETAs VIA CALLING FIRST COMPONENT OF glmDS
 
   cally1 <- call('glmDS1', formula, family, data)
   
   study.summary <- datashield.aggregate(datasources, cally1)
#  num.par.glm<-study.summary$study1$dimX[2]
   num.par.glm<-study.summary[[1]][[1]][[2]]

   beta.vect.next <- rep(0,num.par.glm)
   beta.vect.temp <- paste0(as.character(beta.vect.next), collapse=",")
 

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
    
    message("Iteration ", iteration.count, "...")
   
#NOW CALL SECOND COMPONENT OF glmDS TO GENERATE SCORE VECTORS AND INFORMATION MATRICES
    cally2 <- call('glmDS2', formula, family, beta.vect=beta.vect.temp, offset, weights, data)

      study.summary <- datashield.aggregate(datasources, cally2)

      
      .select <- function(l, field) {
      lapply(l, function(obj) {obj[[field]]})
    }
    
    info.matrix.total<-Reduce(f="+", .select(study.summary, 'info.matrix'))
    score.vect.total<-Reduce(f="+", .select(study.summary, 'score.vect'))
    dev.total<-Reduce(f="+", .select(study.summary, 'dev'))

    message("CURRENT DEVIANCE:      ", dev.total)
    
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
	if(iteration.count==1)
	{
	 beta.vect.next<-rep(0,length(beta.update.vect))
	}
	
	beta.vect.next<-beta.vect.next+beta.update.vect

 	beta.vect.temp <- paste0(as.character(beta.vect.next), collapse=",")


    
    #Calculate value of convergence statistic and test whether meets convergence criterion
    converge.value<-abs(dev.total-dev.old)/(abs(dev.total)+0.1)
    if(converge.value<=epsilon)converge.state<-TRUE
    if(converge.value>epsilon)dev.old<-dev.total
    
    if(viewIter){
      #For ALL iterations summarise model state after current iteration
      message("SUMMARY OF MODEL STATE after iteration ", iteration.count)
      message("Current deviance ", dev.total," on ",(nsubs.total-length(beta.vect.next)), " degrees of freedom")
      message("Convergence criterion ",converge.state," (", converge.value,")")
      
      message("\nbeta: ", paste(as.vector(beta.vect.next), collapse=" "))
      
      message("\nInformation matrix overall:")
      message(paste(capture.output(info.matrix.total), collapse="\n"))
      
      message("\nScore vector overall:")
      message(paste(capture.output(score.vect.total), collapse="\n"))
      
      message("\nCurrent deviance: ", dev.total, "\n")
    }
  }
  if(!viewIter){
    #For ALL iterations summarise model state after current iteration
    message("SUMMARY OF MODEL STATE after iteration ", iteration.count)
    message("Current deviance ", dev.total," on ",(nsubs.total-length(beta.vect.next)), " degrees of freedom")
    message("Convergence criterion ",converge.state," (", converge.value,")")
    
    message("\nbeta: ", paste(as.vector(beta.vect.next), collapse=" "))
    
    message("\nInformation matrix overall:")
    message(paste(capture.output(info.matrix.total), collapse="\n"))
    
    message("\nScore vector overall:")
    message(paste(capture.output(score.vect.total), collapse="\n"))
    
    message("\nCurrent deviance: ", dev.total, "\n")
  }
  
  #If convergence has been obtained, declare final (maximum likelihood) beta vector,
  #and calculate the corresponding standard errors, z scores and p values
  #(the latter two to be consistent with the output of a standard GLM analysis)
  #Then print out final model summary
  if(converge.state)
  {
    family.identified<-0
    
    beta.vect.final<-beta.vect.next
    
    scale.par <- 1
    if(f$family== 'gaussian') {
      scale.par <- dev.total / (nsubs.total-length(beta.vect.next))
    }
    
    family.identified<-1
    se.vect.final <- sqrt(diag(variance.covariance.matrix.total)) * sqrt(scale.par)
    z.vect.final<-beta.vect.final/se.vect.final
    pval.vect.final<-2*pnorm(-abs(z.vect.final))
    parameter.names<-names(score.vect.total[,1])
    model.parameters<-cbind(beta.vect.final,se.vect.final,z.vect.final,pval.vect.final)
    dimnames(model.parameters)<-list(parameter.names,c("Estimate","Std. Error","z-value","p-value"))
    
    if(CI > 0)
    {
      ci.mult <- qnorm(1-(1-CI)/2)
      low.ci.lp <- model.parameters[,1]-ci.mult*model.parameters[,2]
      hi.ci.lp <- model.parameters[,1]+ci.mult*model.parameters[,2]
      estimate.lp <- model.parameters[,1]
      
      
      
      if(family=="gaussian"){
        estimate.natural <- estimate.lp
        low.ci.natural <- low.ci.lp
        hi.ci.natural <- hi.ci.lp
        name1 <- paste0("low",CI,"CI")
        name2 <- paste0("high",CI,"CI")
        ci.mat <- cbind(low.ci.lp,hi.ci.lp)
        dimnames(ci.mat) <- list(NULL,c(name1,name2))   
      }
      
      if(family=="binomial"){
        family.identified  <-  1
        num.parms <- length(low.ci.lp)
        estimate.natural <- exp(estimate.lp)/(1+exp(estimate.lp))
        low.ci.natural <- exp(low.ci.lp)/(1+exp(low.ci.lp))
        hi.ci.natural <- exp(hi.ci.lp)/(1+exp(hi.ci.lp))
        if(num.parms > 1){
          estimate.natural[2:num.parms] <- exp(estimate.lp[2:num.parms])
          low.ci.natural[2:num.parms] <- exp(low.ci.lp[2:num.parms])
          hi.ci.natural[2:num.parms] <- exp(hi.ci.lp[2:num.parms])
          name1 <- paste0("low",CI,"CI.LP")
          name2 <- paste0("high",CI,"CI.LP")
          name3 <- paste0("P_OR")
          name4 <- paste0("low",CI,"CI.P_OR")
          name5 <- paste0("high",CI,"CI.P_OR")
        }       
        ci.mat <- cbind(low.ci.lp,hi.ci.lp,estimate.natural,low.ci.natural,hi.ci.natural)
        dimnames(ci.mat) <- list(NULL,c(name1,name2,name3,name4,name5))
        
      }
      
      if(family=="poisson"){
        family.identified <- 1
        num.parms <- length(low.ci.lp)
        estimate.natural <- exp(estimate.lp)
        low.ci.natural <- exp(low.ci.lp)
        hi.ci.natural <- exp(hi.ci.lp)
        name1 <- paste0("low",CI,"CI.LP")
        name2 <- paste0("high",CI,"CI.LP")
        name3 <- paste0("EXPONENTIATED RR")
        name4 <- paste0("low",CI,"CI.EXP")
        name5 <- paste0("high",CI,"CI.EXP")
        ci.mat <- cbind(low.ci.lp,hi.ci.lp,estimate.natural,low.ci.natural,hi.ci.natural)
        dimnames(ci.mat) <- list(NULL,c(name1,name2,name3,name4,name5))        
      }
      
      if(family.identified==0)
      {
        estimate.natural <- estimate.lp
        low.ci.natural <- low.ci.lp
        hi.ci.natural <- hi.ci.lp
        name1 <- paste0("low",CI,"CI")
        name2 <- paste0("high",CI,"CI")
        ci.mat <- cbind(low.ci.lp,hi.ci.lp)
        dimnames(ci.mat) <- list(NULL,c(name1,name2))   
      }
      
    }
    
    model.parameters<-cbind(model.parameters,ci.mat)
    
    if(!is.null(offset)&&!is.null(weights)){
       formulatext <- paste0(Reduce(paste, deparse(formula)), paste0(" + offset(", offset, ")"), paste0(" + weights(", weights, ")"))
       }

    if(!is.null(offset)&&is.null(weights)){
       formulatext <- paste0(Reduce(paste, deparse(formula)), paste0(" + offset(", offset, ")"))
       }
    if(is.null(offset)&&!is.null(weights)){
       formulatext <- paste0(Reduce(paste, deparse(formula)), paste0(" + weights(", weights, ")"))
       }

    if(is.null(offset)&&is.null(weights)){
       formulatext <- Reduce(paste, deparse(formula))
       }
    
    
      glmds <- list(
      formula=formulatext,
	family=f,
      coefficients=model.parameters,
      dev=dev.total,
      df=(nsubs.total-length(beta.vect.next)),
 	nsubs=nsubs.total,
      iter=iteration.count
    )
    
    class(glmds) <- 'glmds'
    
    return(glmds)
  } else {
    warning(paste("Did not converge after", maxit, "iterations. Increase maxit parameter as necessary."))
    return(NULL)
  }
  
}

