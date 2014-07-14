#' 
#' @title Fits a Generalized Linear Model (GLM) model
#' @description A function that fits generalized linear models
#' @details It enables a parallelized analysis of individual-level data sitting 
#' on distinct servers by sending commands to each data computer to fit a regression 
#' model. The estimates returned are then combined and updated coefficients estimate sent
#' back for a new fit. This iterative process goes on until convergence is achieved.
#' @param x the name if any of the data frame that hold the variables in the regression formula
#' @param formula an object of class \code{formula} which describes the model to be fitted
#' @param family a character, the description of the error distribution function to use in the model
#' @param startCoeff a numeric vector, the starting values for the beta coefficients. If starting
#' values are not provided they are set to 0 for each beta at the start of the iterations.
#' @param maxit the number of iterations of IWLS used
#' instructions to each computer requesting non-disclosing summary statistics.
#' The summaries are then combined to estimate the parameters of the model; these
#' parameters are the same as those obtained if the data were 'physically' pooled
#' @param CI a numeric, the confidence interval.
#' @param viewIter a boolean, tells whether the results of the intermediate iterations
#' should be printed on screen or not. Default is FALSE (i.e. only final results are shown)
#' @param datasources a list of opal object(s) obtained after login to opal servers;
#' these objects also hold the data assigned to R, as a \code{dataframe}, from opal datasources
#' @return coefficients a named vector of coefficients
#' @return residuals the 'working' residuals, that is the residuals in the final
#' iteration of the IWLS fit.
#' @return fitted.values the fitted mean values, obtained by transforming the
#' linear predictors by the inverse of the link function
#' @return rank the numeric rank of the fitted linear model
#' @return family the \code{family} object used
#' @return linear.predictors the linear fit on link scale
#' @return aic A version of Akaike's An Information Criterion, which tells how 
#' well the model fits
#' @author Burton,P;Gaye,A;Laflamme,P
#' @export
#' @examples {
#' 
#' # load the file that contains the login details
#' data(logindata)
#' 
#' # login and assign some variables to R
#' myvar <- list('DIS_DIAB', 'PM_BMI_CONTINUOUS', 'LAB_HDL', 'PM_BMI_CATEGORICAL', 'GENDER')
#' opals <- datashield.login(logins=logindata,assign=TRUE,variables=myvar)
#' 
#' # Example 1: run a GLM without interaction (e.g. diabetes prediction using BMI and HDL levels and GENDER)
#' mod <- ds.glm(x='D', formula=D$DIS_DIAB~D$PM_BMI_CONTINUOUS+D$LAB_HDL+D$GENDER,family='binomial')
#'  
#' # Example 2: run the above GLM model with an intercept (eg. intercept = 1)
#' mod <- ds.glm(x='D', formula=D$DIS_DIAB~1+D$PM_BMI_CONTINUOUS+D$LAB_HDL+D$GENDER,family='binomial')
#'  
#' # Example 3: run the above GLM with interaction HDL and GENDER
#' mod <- ds.glm(x='D', formula=D$DIS_DIAB~D$PM_BMI_CONTINUOUS+D$LAB_HDL*D$GENDER,family='binomial')
#'  
#' # Example 4: now run a GLM but with interaction where the error follows a poisson distribution
#' mod <- ds.glm(x='D', formula=D$PM_BMI_CATEGORICAL~D$PM_BMI_CONTINUOUS+D$LAB_HDL+D$GENDER,family='poisson')
#'  
#' # clear the Datashield R sessions and logout
#' datashield.logout(opals) 
#'  
#' }
#' 
#' @references Jones EM, 'DataSHIELD-shared individual-level analysis without sharing the data: a biostatistical
#' perspective', Norsk epidemiologi - Norwegian journal of epidemiology 2012;21(2): 231-9.
#' 
ds.glm <- function(x=NULL, formula=NULL, family=NULL, startCoeff=NULL, maxit=15, CI=0.95, viewIter=FALSE, datasources=NULL) {
  
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
  
  # check if user have provided a formula
  if(is.null(formula) | class(formula) != 'formula'){
    stop("Please provide a valid formula!", call.=FALSE)
  }
  
  # family 
  families <- c("binomial", "gaussian", "Gamma", "poisson")
  if(is.null(family) | !(family %in% families)){
    stop("Please provide a valid 'family' parameter: 'binomial', 'gaussian', 'Gamma' or 'poisson'.", call.=FALSE)
  }
  
  # if the user provides starting values check that vector is of the same length than the variables in the formula
  if(!(is.null(startCoeff))){
    l1 <- length(startCoeff)
    l2 <- length(all.vars(formula))
    if(l1 != l2){
      stop("The number starting beta values must be the same as the number of variables in the formula!", call.=FALSE)
    }
  }
  
  # check if all the variables in the lp formula are exist on the server site and if any is empty
  message("Checking the input variables are defined and in the right format...")
  stdnames <- names(datasources)
  temp <- glmhelper2(formula)
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
        datashield.assign(datasources[j], as.character(variables[[i]]), quote(as.numeric(as.character(variables[[i]]))))
      }else{
        if(is.null(x)){
          stop("x=NULL! You must set the argument x to ", inputterms[1], " because the variable ", inputterms[2], " is indicated as attached ", inputterms[1], " in ", as.character(variables[[i]]), "!", call.=FALSE)
        }
        cally2 <- call('isNaDS', variables[[i]])
        d2 <- datashield.aggregate(datasources[j], cally2)[[1]]
        if(d2){ 
          stop("The variable ", as.character(variables[[i]]), " in ", stdnames[j], " is empty (all values are 'NA').", call.=FALSE)
        }
        # turn the vector into numeric
        cally <- call("as.numeric",call("as.character", variables[[i]]))
        datashield.assign(datasources[j], inputterms[2], cally)
      }
    }
  }
  
  # the variables have been turned into numeric to avoid issues with factors if for example
  # a poisson distribution is used for the fitted model. The numeric variables were saved 
  # with the same names; now we re-construct the linear predictor formula using the saved 
  # loose variables (i.e. without the 'dataFrameName$' bit if that was how the formula was given)
  formulaChr <- as.character(formula)
  # a formula has always 3 parts: '~', the outcome and the covariates part
  bits <- c()
  for(i in 2:3){
    bit <- gsub(paste0(x,"([$])"), '', formulaChr[i])
    bits <- append(bits, bit)
  }
  formula <- as.formula(paste0(bits[1], "~", bits[2]))
    
  # number of 'valid' studies (those that passed the checks) and vector of beta values
  numstudies <- length(datasources)
  beta.vect.next <- startCoeff
  
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
    
    if(iteration.count == 0){
      if(is.null(beta.vect.next)){
        beta.vect.temp <- NULL
      }else{
        beta.vect.temp <- paste0(startCoeff, collapse=",")
      }
    }else{
      beta.vect.temp <- paste0(beta.vect.next, collapse=",")
    }

    # call for parallel glm and retrieve results when available
    cally <- as.call(list(quote(glmDS), formula, family, beta.vect.temp))
    study.summary <- datashield.aggregate(datasources, cally)
    
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
    
    # update terms to current beta vector to obtain new beta vector for next iteration
    if(iteration.count == 0) {
      beta.vect.next<-beta.update.vect
    } else {
      beta.vect.next <- beta.vect.next+beta.update.vect
    }
    iteration.count<-iteration.count+1
    message("Iteration ", iteration.count, "...")
    
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
    family.identified <- 0
    beta.vect.final <- beta.vect.next
    
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
        low.ci.natural <- exp(low.ci.lp)/(1+low.ci.lp)
        hi.ci.natural <- exp(hi.ci.lp)/(1+hi.ci.lp)
        if(num.parms > 1){
          estimate.natural[2:num.parms] <- exp(estimate.lp[2:num.parms])
          low.ci.natural[2:num.parms] <- exp(low.ci.lp[2:num.parms])
          hi.ci.natural[2:num.parms] <- exp(hi.ci.lp[2:num.parms])
          name1 <- paste0("low",CI,"CI.LP")
          name2 <- paste0("high",CI,"CI.LP")
          name3 <- paste0("P_or_OR")
          name4 <- paste0("low",CI,"CI.P_OR")
          name5 <- paste0("high",CI,"CI.P_OR")
        }       
        ci.mat <- cbind(low.ci.lp,hi.ci.lp,estimate.natural,low.ci.natural,hi.ci.natural)
        dimnames(ci.mat) <- list(NULL,c(name1,name2,name3,name4,name5))
        
      }
      
      if(family=="poisson"){
        family.identified <- 1
        num.parms <- length(low.ci.lp)
        estimate.natural <- exp(estimate.lp)/(1+exp(estimate.lp))
        low.ci.natural <- exp(low.ci.lp)/(1+low.ci.lp)
        hi.ci.natural <- exp(hi.ci.lp)/(1+hi.ci.lp)
        estimate.natural[2:num.parms] <- exp(estimate.lp[2:num.parms])
        low.ci.natural[2:num.parms] <- exp(low.ci.lp[2:num.parms])
        hi.ci.natural[2:num.parms] <- exp(hi.ci.lp[2:num.parms])
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
    
    model.parameters <- cbind(model.parameters,ci.mat)
    
    glmds <- list(
      formula=formula,
      coefficients=model.parameters,
      dev=dev.total,
      nsubs=nsubs.total,
      df=(nsubs.total-length(beta.vect.next)),
      iter=iteration.count
    )
    
    class(glmds) <- 'glmds'
    
    return(glmds)
  } else {
    warning(paste("Did not converge after", maxit, "iterations. Increase maxit parameter as necessary."))
    return(NULL)
  }
  
}
