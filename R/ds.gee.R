#' 
#' @title Fits a Generalized Estimating Equation (GEE) model
#' @description A function that fits generalized estimated equations to deal with correlation 
#' structures arising from repeated measures on individuals, or from clustering as in family data. 
#' @details It enables a parallelized analysis of individual-level data sitting on distinct servers 
#' by sending commands to each data computer to fit a GEE model model. The estimates returned are 
#' then combined and updated coefficients estimate sent back for a new fit. This iterative process 
#' goes on until convergence is achieved. The input data should not contain missing values.
#' The data must be in a data.frame obejct and the variables must be refer to through the data.frame.
#' @param data the name of the data frame that hold the variables in the regression formula.
#' @param formula a string character, the formula which describes the model to be fitted.
#' @param family a character, the description of the error distribution:  'binomial', 'gaussian', 
#' 'Gamma' or 'poisson'.
#' @param corStructure a character, the correlation structure: 'ar1', 'exchangeable', 'independence', 
#' 'fixed' or 'unstructure'.
#' @param clusterID a character, the name of the column that hold the cluster IDs
#' @param startCoeff a numeric vector, the starting values for the beta coefficients. 
#' @param userMatrix a list of user defined matrix (one for each study). These matrices are 
#' required if the correlation structure is set to 'fixed'.
#' @param maxit an integer, the maximum number of iteration to use for convergence.
#' @param display a boolean to display or not the intermediate results. Default is FALSE.
#' @param datasources a list of opal object(s) obtained after login to opal servers;
#' these objects also hold the data assigned to R, as a \code{dataframe}, from opal datasources.
#' @return a list which contains the final coefficient estimates (beta values), the pooled alpha
#' value and the pooled phi value.
#' @author Gaye, A.; Jones EM.
#' @export
#' @examples {
#' 
#'   # load the login data file for the correlated data
#'   data(geeLoginData)
#'   
#'   # login and assign all the stored variables
#'   opals <- datashield.login(logins=geeLoginData,assign=TRUE)
#'   
#'   # set some parameters for the function 9the rest are set to default values)
#'   myformula <- 'response~1+sex+age.60'
#'   myfamily <- 'binomial'
#'   startbetas <- c(-1,1,0)
#'   clusters <- 'id'
#'   mycorr <- 'ar1'
#'   
#'   # run a GEE analysis with the above specifed parameters
#'   xx <- ds.gee(data='D',formula=myformula,family=myfamily,corStructure=mycorr,clusterID=clusters,startCoeff=startbetas)
#'   
#' # clear the Datashield R sessions and logout
#' datashield.logout(opals) 
#' 
#' }
#' 
#' @references Jones EM, Sheehan NA, Gaye A, Laflamme P, Burton P. Combined analysis of correlated data
#' when data cannot be pooled. Stat 2013; 2: 72-85.
#'
ds.gee <- function(data=NULL, formula=NULL, family=NULL, corStructure='ar1', clusterID=NULL, startCoeff=NULL, userMatrix=NULL, 
                   maxit=20, display=FALSE, datasources=NULL){
  
  # turn the input formula into an object of type 'formula', it is given as a string character
  formula <- as.formula(formula)
  
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
  
  # check if user have provided the name of the dataset and if the dataset is defined
  if(is.null(data)){
    stop("data=NULL; please provide the name of the dataset that holds the variables!", call.=FALSE)
  }else{
    defined <- isDefined(datasources, data)
  }
  
  # check if user have provided a formula
  if(is.null(formula) | class(formula) != 'formula'){
    stop("Please provide a valid formula!", call.=FALSE)
  }
  
  # check if user have provided a formula
  if(is.null(clusterID)){
    stop("Please provide the name of the column that holds the cluster IDs!", call.=FALSE)
  }
  
  # if no start values have been provided by the user throw an alert and stop process.
  # it is possible to set all betas to 0 here but sometimes that can cause the program
  # to crash so it is safer to let the use choose sensible start values
  l1 <- length(startCoeff)
  l2 <- length(all.vars(formula))
  if(is.null(startCoeff)) {
    message("No starting values provided for the beta coefficients. The starting values will be set to 0 each.", call.=FALSE)
    startCoeff <- rep(0, l2)
  }else{
    if(l1 != l2){
      stop("The number starting beta values is incorrect!", call.=FALSE)
    }
  }
  
  # correlation structure
  corStructures <- c('ar1', 'exchangeable', 'independence', 'fixed', 'unstructure')
  if(is.null(corStructure) | !(corStructure %in% corStructures)){
    stop("Please provide a valid 'corStructure' parameter: 'ar1', 'exchangeable', 'independence', 'fixed' or 'unstructure'.", call.=FALSE)
  }
  # if the correlation structure is set to "fixed" and the user has not provided 'user defined matrices"
  # or the user has not supplied one matrix for each study throw an alert and stop process
  if(corStructure == "fixed"){
    if(is.null(userMatrix) | length(userMatrix) < length(datasources)) {
      stop("Please provide one correlation matrix for each of the studies!", call.=FALSE)
    }
  }
  
  # family 
  families <- c("binomial", "gaussian", "Gamma", "poisson")
  if(is.null(family) | !(family %in% families)){
    stop("Please provide a valid 'family' parameter: 'binomial', 'gaussian', 'Gamma' or 'poisson'.", call.=FALSE)
  }
  
  # check if the input dataframe is complete and if the variables in the lp formula are in the input dataframe
  cally <- call('complete.cases', as.symbol(data))
  datashield.assign(datasources, 'Dcomplete', cally)
  
  stdnames <- names(datasources)
  temp <- glmhelper2(formula)
  variables <- c()
  for(i in 1:length(temp)){
    variables <- append(variables, temp[[i]])
  }
  colsD <- datashield.aggregate(datasources, paste0("colnames(", x, ")"))[[1]]
  for(i in 1:length(variables)){
    for(j in 1: length(datasources)){
      lengthDcomplete <- datashield.aggregate(datasources[j],paste0("length(Dcomplete)"))[[1]]
      nrowD <- datashield.aggregate(datasources[j], paste0("dim(", data, ")"))[[1]][1]
      if(lengthDcomplete != nrowD){
        stop("The input dataset ", data,  " in ", stdnames[j] , " contains one or more missing values. Only complete datasets are allowed in GEE analysis.", call.=FALSE)
      }else{
        inputterms <- unlist(strsplit(deparse(variables[[i]]), "\\$", perl=TRUE))
        if(length(inputterms) > 1){
          if(!(inputerms[2] %in% colsD)){
            stop("The variable ", as.character(variables[[i]]),  " is not in the dataset ", data, " in ", stdnames[j], call.=FALSE)
          }
        }else{
          if(!(as.character(variables[[i]]) %in% colsD)){
            stop("The variable ", as.character(variables[[i]]),  " is not in the dataset ", data, " in ", stdnames[j], call.=FALSE)
          }
        }
      }
    }
  }
  
  # loop until convergence is achieved or the maximum number of iterations is reached  
  for(r in 1:maxit){
    
    if(r == 1){ 
      startCoeff <- startCoeff
      betas <- paste(startCoeff,collapse=',')
      beta1 <- startCoeff[1]
    }else{
      startCoeff <- beta.vector[[1]]
      betas <- paste(startCoeff,collapse=',')
      beta1 <- beta.vector[[1]][1]
    }
    
    # vectors to hold relevant output
    alphaMs <- vector("list", length(datasources))
    Ns <- c()
    Ms <- c()
    sum_ps <- c()
    
    # call the server side funtion 'alphaPhiDS'
    for(i in 1:length(datasources)){
      # run function and store relevant output
      cally <- as.call(list(quote(alphaPhiDS), as.symbol(data), formula, family, clusterID,corStructure,betas))
      temp <- datashield.aggregate(datasources[i], cally)
      output1 <- temp[[1]]
      Ns <- append(Ns, output1$N)
      num.para <- output1$npara
      Ms <- append(Ms, output1$M_study)      
      alphaMs[[i]] <- output1$alphaM
      sum_ps <- append(sum_ps, output1$sum_p) 
    }    

    # run 'geehelper1' using the output values from 'alphaphiDS' as arguments
    output2 <- geehelper1(N=Ns, npara=num.para, M_study=Ms, alphaM=alphaMs, sum_p=sum_ps, corStructure=corStructure)
    alpha.comb <- output2$alpha
    phi.comb <- output2$phi    
    
    # list to hold relevant output
    score.vects <- vector("list", length(datasources)) 
    info.mats <- vector("list", length(datasources)) 
    J.mats <- vector("list", length(datasources))
    
    # # call the server side funtion 'scoreVectDS'
    for(i in 1:length(datasources)){
      userMat <- userMatrix[[i]]
      cally <- as.call(list(quote(scoreVectDS), as.symbol(data), formula, family, clusterID, corStructure, alpha.comb,phi.comb,betas,userMat))
      temp <- datashield.aggregate(datasources[i], cally)
      output3 <- temp[[1]]
      score.vects[[i]] <- output3$score.vector
      info.mats[[i]] <- output3$info.matrix
      J.mats[[i]] <- output3$J.matrix
    }
    
    # call the internal function 'geehelper2'
    beta.vector <- geehelper2(score=score.vects, infoMatrix=info.mats, Jmatrix=J.mats, startCoeff=startCoeff)
    
    # calculate convergence term
    beta2 <- beta.vector[[1]][1]
    xx <- abs(beta2 - beta1)
    
    # print summaries
    message(paste0(" Iteration ", r))
    if(display){
     for(tt in 1:length(startCoeff)){
        if(tt==1){message(" ")}
        if(xx < 0.000001 | r == maxit){
          message(paste0("beta",tt-1 , ": ", round(startCoeff[tt],5), ", Std.Error =", round(beta.vector[[2]][tt],5)))
        }else{
          message(paste0("beta",tt-1 , ": ", round(startCoeff[tt],5)))
        }
      }
    }
    if(xx < 0.000001 | r == maxit){
      message(paste0(" Convergence value = ", xx)) 
      message(paste0(" Converged = ", xx < 0.000001))   
      break
    }
  }

  # finalize output
  betaValues <- round(beta.vector$beta.values,5)
  stdErrors <- round(beta.vector$standard.errors,5)
  mainOutput <- data.frame(cbind(betaValues,stdErrors))
  colnames(mainOutput) <- c("coefficients", "standard.errors")
  rownames(mainOutput) <- rownames(betaValues)
  
  return(list(estimates=mainOutput, alpha=round(alpha.comb,3), phi=round(phi.comb,3)))
}
