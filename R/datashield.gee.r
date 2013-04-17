#' This function allows for the fitting of generalized estimating equations
#' to correlated and non-pooled data. It enables a parallelized analysis of
#' individual-level data sitting on distinct computers/servers by sending 
#' instructions to each computer requesting non-disclosing summary statistics.
#' The sumaries are then combined to estimate the parameters of the model; these 
#' parameters are the same as those obtained if the data were 'physically' pooled.
#'
#' @title Runs a combined GEE analysis of non-pooled data
#'
#' @param opals character strings that represent the URL of the servers where 
#' the study datasets are stored
#' @param formula an object of class \code{formula} which describes the model to be fitted
#' @param cluster a vector that holds the clusters (i.e.family IDs)
#' @param family a description of the error distribution function to use in the model
#' @param link a description of the link function to be used in the model
#' @param data a list of character strings that represent the names of the datasets
#' @param corstr a character string that specifies the correlation structure to use
#' @param start.betas starting values for the parameters estimates
#' @param user.mats user-defined correlation matrices (one per dataset) used 
#' if the paramter \code{corstr} is set to 'fixed'
#' @param control a list of parameters for controlling the fitting process
#' @details This function requires to log into the opal servers where the data 
#' are stored (one dataset per study). Set the parameter \code{corstr} to set 
#' to 1, 2, 3 , 4 or 5 for respectively 'ar1', 'exchangeable', 'independence', 
#' 'fixed', or 'unstructured' correlation structures. Set the parameter \code{link}
#' to 1, 2, 3, 4, or 5 for respectively the 'logit', 'identity', 'inverse', 'log' or 
#' 'probit' link functions. The \code{binomial family} accepts the links 'logit', 
#' 'probit' and 'log'; the \code{gaussian family} accepts the links 'identity', 'log' 
#' and 'inverse' and the \code{Gamma family} acepts the links 'inverse', 'identity' and 'log'.
#' @return coefficients estimates of the parameters of the fitted model
#' @return stderrs standard errors of the estimates
#' @return alpha correlation parameter, a correlation between any two correlated observations
#' @return phi dispersion parameter
#' @return working.corr.matrix the working correlation matrix
#' @references Zeger, S. L., Liang, K. and Albert, P. S. (1988). 
#' Models for longitudinal data: a generalized estimating equation approach. 
#' Biometrics, 44, p1049-1060.
#' @examples
#' \dontrun{
#' # put example here
#'}
#' @export
#'
datashield.gee <- function(opals, formula, cluster, family, link=1, data=data, corstr=1, start.betas=NULL, user.mats=NULL, control=list(maxit=10, seed=17632, display=FALSE)){
    
    ######## FEW CHECKS TO ENSURE ALL REQUIRED ARGUMENTS HAVE BEEN SUPPLIED ########
    # check if user have provided the column names of the response, cluster IDs, and covariates
    # if not throw an alert and stop process
    if(is.null(formula)){
      stop(call.=FALSE, "\nPLEASE PROVIDE THE FORMULA\n")
    }
    # if no start values have been provided by the user throw an alert and stop process.
    # it is possible to set all betas to 0 here but sometimes that can cause the program
    # to crash so it is safer to let the use choose sensible start values
    if(is.null(start.betas)) {
      stop(call.=FALSE, "\nPLEASE PROVIDE INITIAL VALUES FOR THE BETA COEFFICIENTS\n")
    }
    
    # if the correlation structure is set to "fixed" and the user has not provided 'user defined matrices"
    # or the user has not supplied one matrix for each study throw an alert and stop process
    if(corstr == 4) {
      if(is.null(user.mats) | length(user.mats) < length(data)) {
        stop(call.=FALSE, "\n", paste("PLEASE PROVIDE ONE CORRELATION MATRIX FOR EACH OF THE ", length(opals), " STUDIES", sep=""), "\n")
      }
    }
    
    cor.call <- corstr
     
    beta.vector <- start.betas
    for(r in 1:control$maxit){
      set.seed(control$seed)
      if(r == 1){ 
        start.betas <- start.betas
        initial.betas <- start.betas
        beta1 <- start.betas[1]
      }else{
        start.betas <- beta.vector[[1]]
        beta1 <- beta.vector[[1]][1]
      }
      ############ LOOP THROUGH THE STUDIES AND RUN gee.alphaphi.ds ##############
      # vectors to hold relevant output
      alphaMs <- vector("list", length(opals))
      Ns <- c()
      Ms <- c()
      sum_ps <- c()
      
      call <- as.call(list(quote(gee.alphaphi.ds), formula, family, as.vector(link), data, cluster, as.vector(cor.call), as.vector(start.betas)));
      # run 'gee.alphaphi.ds'
      study.alphaphi<-datashield.aggregate(opals, call);
      for(i in 1:length(opals)){
        # run function and store relevant output
        output1 <- study.alphaphi[[i]]
        alphaMs[[i]] <- output1$alphaM
        Ns <- append(Ns, output1$N)
        Ms <- append(Ms, output1$M_study)
        sum_ps <- append(sum_ps, output1$sum_p) 
        num.para <- output1$npara
      }
      
      ############ COMBINE OUTPUT VALUES OF THE DIFFERENT STUDIES ##############
      # run 'alpha.phi.ds' using the output values from 'gee.alphaphi.ds' as arguments
      output2 <- alpha.phi.ds(N=Ns, npara=num.para, M_study=Ms, alphaM=alphaMs, sum_p=sum_ps)
      
      alpha.comb <- output2$alpha
      phi.comb <- output2$phi
      
      
      ############ LOOP THROUGH THE STUDIES AND RUN gee.beta.ds ##############
      # list to hold relevant output
      score.vects <- vector("list", length(opals)) 
      info.mats <- vector("list", length(opals)) 
      J.mats <- vector("list", length(data))
      work.cor.mats <- vector("list", length(data))
      
      for(i in 1:length(opals)){
        user.mat<-user.mats[[i]]
        makeMat<-as.call(list(quote(matrix), as.vector(user.mat), attributes(user.mat)$dim))
        call<-as.call(list(quote(gee.beta.ds), formula, family, as.vector(link), data, cluster, alpha.comb, phi.comb, as.vector(cor.call), as.vector(start.betas), makeMat));
        # run 'gee.beta.ds'
        output3 <- datashield.aggregate(opals[[i]], call);
        score.vects[[i]] <- output3$score.vector
        info.mats[[i]] <- output3$info.matrix
        J.mats[[i]] <- output3$J.matrix
        work.cor.mats[[i]] <- output3$working.cor.matrix
        #cat("\n")
        #print(work.cor.mats[[i]])
        #cat("\n")
      }
      
      # get the working correlation matrix to output
      diag.size <- rep(0, length(work.cor.mats))
      for(i in 1:length(work.cor.mats)){
        diagonal <- diag(work.cor.mats[[i]])
        diag.size[i] <- length(diagonal)
      }
      diag.index <- which(diag.size == max(diag.size, na.rm=TRUE))
      work.corr.mat <- work.cor.mats[[diag.index[1]]]
      #cat("\n", diag.size, "\n")
      
      ############ GENERATE THE FINAL OUTPUT (COMBINE OUTPUT ALL THE ANALYSES - ONE ANALYSIS PER STUDY)##############
      # run 'beta.vects.ds'
      beta.vector <- beta.vect.ds(score=score.vects, infomatrix=info.mats, J.matrix=J.mats, start.betas=start.betas)
      
      
      ########### CALCULATE CONVERGENCE TERM ############
      beta2 <- beta.vector[[1]][1]
      xx <- abs(beta2 - beta1)
      
      ############ PRINT SUMMARIES ##############
      if(control$display){
        #if(r==1){cat("\n Start beta values=", initial.betas, "\n")}
        cat("\n Iteration", r)
        for(tt in 1:length(start.betas)){
          if(tt==1){cat("\n")}
          if(xx < 0.000001 | r == control$maxit){
            cat("beta",tt-1 , "= ", round(start.betas[tt],5), ", Std.Error =", round(beta.vector[[2]][tt],5), "\n")
          }else{
            cat("beta",tt-1 , "= ", round(start.betas[tt],5), "\n")
          }
        }
      }
      
      
      if(xx < 0.000001 | r == control$maxit){
        cat("\n Convergence value =", xx, "\n") 
        cat(" Converged = ", xx < 0.000001, "\n")   
        break
      }
    }
    
    if(corstr == 4 ){ # user defined
      alpha.val <- work.corr.mat[col(work.corr.mat) < row(work.corr.mat)]
    }else{
      if(corstr == 3){
        alpha.val <- 0
      }else{
        alpha.val <- alpha.comb
      }
    }
    #cat("\n\n alpha=", alpha.val, "\n\n")
    final.betas <- beta.vector[[1]]
    stderrors <- beta.vector[[2]]
    return(list(coefficients=round(final.betas,5), stderrs=round(stderrors,5), alpha=round(alpha.val,5), phi=round(phi.comb,5), working.corr.matrix=work.corr.mat))
  }
  