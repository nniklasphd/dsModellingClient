#' 
#' @title Generates an expanded version of a dataset that contains survival data
#' @description This function is meant to be used as part of a piecewise regression analysis.
#' @details It splits the survival interval time of subjects into pre-specified sub-intervals that are each assumed
#' to encompass a constant base-line hazard which means a constant instantaneous risk of death). In the expanded dataset
#' a row is included for every interval in which a given individual is followed regardless how short that follow-up may be.
#' Each row includes: (1) a variable indicating failure status for a particular interval in that interval also known
#' as censoring status (1=failed, died, relapsed, developed a disease etc, 0= e.g. lost-to-follow-up or
#' passed right through the interval without failing); (2) an exposure-time variable indicating the duration of
#' exposure-to-risk-of-failure the corresponding individual experienced in that interval before he/she failed or
#' was censored). To illustrate, an individual who survives through 5 such intervals and then dies/fails in the 6th interval
#' will be allocated a 0 value for the failure status/censoring variable in the first five intervals and a 1 value
#' in the 6th, while the exposure-time variable will be equal to the total length of the relevant interval
#' in each of the first five intervals, and the additional length of time they survived in the sixth interval
#' before they failed or were censored. If they survive through the first interval and they are censored in the
#' second interval, the failure-status variable will take the value 0 in both intervals. The expanded data set
#' also includes a unique ID in a form such as 77.13 which identifies that row of the dataset as relating to the
#' 77th individual in the input data set (in whatever order they have been placed by Opal [which is often different to the
#' original numeric order of the IDs that were actually specified to Opal]) and his/her experience (exposure-time and
#' failure status)in the 14th interval. Note that .N indicates the N+1th interval because interval 1 has no
#' suffix. The expanded data set also includes the original ID as specified to Opal,
#' the total survival time (equivalent to the sum of the exposure times across all intervals) and the ultimate
#' failure-status in the final interval to which they were exposed. The value of each of these variables is repeated in
#' every row to which the corresponding individual contributes data (i.e. to every row corresponding to an interval
#' in which that individual was followed). Finally, the expanded dataset contains any other variables pertaining
#' to each individual that the user would like to  carry forwarded to a survival analysis based on the expanded
#' data. The values of each of these variables are also repeated in every row corresponding to an interval in
#' which that individual was followed.
#' @param data a character, the name of a dataframe containing the data to be expanded. Often, the dataframe
#' will also hold the original total-survival-time and final-censoring variables but the lexis function is deliberately
#' set up so those can also be specified as coming either from a different data frame or from the root area of your
#' analysis (i.e. lying outside any dataframe)table that holds the original data, this is the data to be expanded.
#' @param intervalWidth, a numeric vector specifying the length of each interval. If the total sum of the
#' duration across all intervals is less than the maximum follow-up of any individual in any contributing
#' study, a final interval will be added by ds.lexis extending from the end of the last interval specified
#' to the maximum follow-up time. If a single numeric value is specified rather than a vector, ds.lexis
#' will keep adding intervals of the length specified until the maximum follow-up time in any single
#' study is exceeded. 
#' @param idCol a character the name of the column that holds the individual IDs of the subjects. This may
#' be numeric or character. Note that when an particular variable is identified as being the ID to Opal
#' when the data are first transferred to Opal, the ID often ends up being of class character and will typically
#' be sorted in 'alphabetic' order (treating each digit as a character) rather than numeric. So, in a dataset
#' containing sequential IDs 1-1000, the order allocated by Opal will be mean that the first thirteen rows
#' in the original dataset (before expansion) will correspond to the original IDs:
#' 1,10,100,101,102,103,104,105,106,107,108,109,11 (analogous to b, ba, baa, bab, bac, bad, bae, baf, bag,
#' bah, bai, baj, bb ...) in an alphabetic listing: NOT to the expected order 1,2,3,4,5,6,7,8,9,10,11,12,13 ...
#' This alphabetic order or the ID listing will then carry forward to the expanded dataset. But the
#' nature and order of the original ID variable held in idCol doesn't matter to ds.lexis. Provided every
#' individual appears only once in the original data set (before expansion) the order does not matter
#' because ds.lexis works on its own unique numeric vector that is allocated from 1:M (where there are
#' M individuals) in whatever order they appear in the original dataset. It is this ds.lexis sequentially
#' allocated numeric ID vector that is ultimately combined with the interval period number to
#' produce the expanded unique ID variable in the expanded data set (e.g. see 77.13 under 'details')   
#' @param entryCol a character, the name of the column that holds the entry times (i.e. start of follow up).
#' Rather than using a total survival time variable to identify the intervals to which any given individual
#' is exposed, ds.lexis requires an initial entry time and a final exit time. If the data you wish to expand
#' contain only a total survival time variable and (as is most common) every individual starts follow-up
#' at time 0, the entry times should all be specified as zero, and the exit times as the total survival time.
#' So, entryCol should either be the name of the column holding the entry time of each individual, or else
#' if no entryCol is specified it will be defaulted to zero anyway and put into a variable called STARTTIME
#' in the expanded data set.
#' @param exitCol a character, the name of the column that holds the exit times (i.e. end of follow up).
#' If the entry times are set, or defaulted, to zero, the exitCol variable should contain the total
#' survival times.
#' @param statusCol a character, the name of the column that holds the failure/censoring status of each subject
#' (see under details).
#' @param variables a character vector, the column names of the variables to include in the 
#' final expanded table. If the 'variables' argument is not set (is null) but the 'data' argument is set
#' the expanded data set will contain all variables in the dataframe identified by the 'data' argument.
#' If neither the 'data' or 'variables' arguments are set, the expanded data set will only include the
#' ID, exposure time and failure/censoring status variables which may still be useful for plotting survival data
#' once these become available.
#' @param newobj the name of the new dataframe containing the expanded data set. If you specify a name,
#' that name will be used, but if no name is specified it will be defaulted to 'name of dataframe specified by
#' data argument' with '_expanded' as a suffix. If you use the client side function ds.ls() after running
#' ds.lexis the new dataframe you have created should be listed.
#' @param datasources As in all client-side functions, a list of opal object(s) obtained after login to opal servers;
#' these objects also hold the data assigned to R, as a \code{data frame}, from opal datasources
#' @return a dataframe, an expanded version of the input table.
#' @author Gaye, A, Burton PR
#' @seealso \code{ds.glm} for genralized linear models
#' @seealso \code{ds.gee} for generalized estimating equation models
#' @export
#' @examples {
#' EXAMPLE 1
#' In this example, the data to be expanded are held in a dataframe called 'CM'. The survival time intervals are to
#' be 0<t<=2.5; 2.5<t<=5.0, 5.0<t<=7.5, up to the final interval of duration 2.5 that includes the maximum survival time.
#' The original ID, entry-time, exit-time and censoring variables are all included in the original dataframe in variables
#' CM$ID, CM$STARTTIME, CM$ENDTIME and CM$CENS. The expanded dataframe will be created with the name "EM$new" and the
#' data will be held in a dataframe called 'expanded.table' inside EM$new.
#' 
#' ds.lexis(data = "CM", intervalWidth = 2.5, idCol = "CM$ID",
#'   entryCol = "CM$STARTTIME", exitCol = "CM$ENDTIME", statusCol = "CM$CENS",
#'   newobj = "EM.new")
#' ds.ls() #to confirm expanded dataframe created
#' 
#' Please note, if the censoring variable had instead been held in a variable
#' 'CENSOR' (outside any dataframe) or in DF$died (inside a different dataframe called DF) then it would have been perfectly
#' acceptable to specify statusCol="CENSOR" and/or statusCol="DF$died"
#' 
#' For illustration, the following is a schema of a typical set of variables you get
#' in an expanded dataframe. Depending how the arguments of ds.lexis are specified
#' some variables may be repeated:
#' 
#' A			B		C		D	E		F		G		H	I		J	 K		  L		M 
#' 657         	657     2.054   1   657		983  	2  		1   2.054	0     13      1      2.3
#' 658         	658     0.348   0   658     984  	2    	0   0.348	0      8      1     -2.7
#' 659  	  	659		2.500   0   659     985  	2     	1   9.300	0    -21      1     -0.6
#' 659.1		659     2.500   0   659		985  	2      	1   9.300	0    -21      1     -0.6
#' 659.2		659     2.500   0   659		985  	2       1   9.300	0    -21      1     -0.6
#' 659.3		659     1.800 	1   659		985     2       1   9.300	0    -21      1     -0.6
#' 
#' A = Expanded Unique ID - sequential ID allocated by ds.lexis combined with indicator of interval
#' to which that row corresponds. No indicator = interval 1, .N = N+1th interval (e.g. see subject 659)
#' B = Sequential ID allocated by ds.lexis
#' C = Exposure time in corresponding interval (e.g. see subject 659)
#' D = Failure/censoring status in corresponding interval (e.g. see subject 659)
#' E = Repeat of Sequential ID allocated by ds.lexis
#' F = Original ID allocated by Opal (note not in same order as Sequential ID)
#' G = Study indicator (in this example these individuals are in study 2)
#' H = Final failure/censoring status (note gets repeated across all intervals for any individual)
#' I = Total survival time (note gets repeated across all intervals for any individual)
#' J = Additional variable carried into expanded dataframe. In this case the variable is sex (0=male)
#'       note repeated across all intervals for any individual.
#' K = Additional variable carried into expanded dataframe.
#'       note repeated across all intervals for any individual.
#' L = Additional variable carried into expanded dataframe.
#'       note repeated across all intervals for any individual.
#' M = Additional variable carried into expanded dataframe.
#'       note repeated across all intervals for any individual.
#' 
#' 
#' EXAMPLE 2
#' In this example, the survival time intervals are to be 0<t<=1; 1<t<=2.0, 2.0<t<=5.0, 5.0<t<=11.0,
#' then 11.0<t<=maximum survival time in any study (if no survival time exceeds 11, the fifth interval
#' will not appear in any row). No newobj is specified so the output dataframe will be named 'CM_expanded'.
#' ds.lexis(data = "CM", intervalWidth = c(1,1,3,6), idCol = "CM$ID",
#'   entryCol = "CM$STARTTIME", exitCol = "CM$ENDTIME", statusCol = "CM$CENS")
#' ds.ls() #to confirm expanded dataframe created
#' }
#'
ds.lexis <- function(data=NULL, intervalWidth=NULL, idCol=NULL, entryCol=NULL, exitCol=NULL, statusCol=NULL, variables=NULL, newobj=NULL, datasources=NULL){
  
  # if no opal login details are provided look for 'opal' objects in the environment
  if(is.null(datasources)){
    datasources <- findLoginObjects()
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
  
  # if no value provided for 'intervalWidth' instruct user to specify one
  if(is.null(intervalWidth)){
    stop("Please provide a single numeric value or vector to identify the survival time intervals", call.=FALSE)
  }
  
  # if no value spcified for output (expanded) data set, then specify a default
  if(is.null(newobj)){
    newobj <- paste0(data,"_expanded")
  }

#FIRST CALL TO SERVER SIDE TO IDENTIFY THE MAXIMUM FOLLOW UP TIME IN ANY
#SOURCE. THE MAXIMUM I EACH SOURCE IS MASKED BY A RANDOM POSITIVE INCREMENT
  calltext1 <- call("lexisDS1", exitCol)

  maxtime<-datashield.aggregate(datasources, calltext1)

  nummax<-length(maxtime)

  temp1<-rep(NA,nummax)

  for(j in 1:nummax){
	temp1[j]<-unlist(maxtime[[j]][1])
  }

#IDENTIFY MAXIMUM OF THE MAXIMUM FOLLOW-UP TIMES
  maxmaxtime<-max(temp1)
 

#SECOND CALL TO SERVER SIDE USES maxmaxtime AND intervalWidth TO SET
#FOLLOW-UP TIME BREAKS IN EACH STUDY (ALL THE SAME)
  # call the main server side function
  calltext2 <- call("lexisDS2", data, intervalWidth, maxmaxtime, idCol, entryCol, exitCol, statusCol, variables)
  	datashield.assign(datasources, newobj, calltext2)

#RETURN COMPLETION INFORMATION TO CLIENT SIDE
  completion.flag1<-"ASSIGN FUNCTION COMPLETED - USE ds.ls() TO CONFIRM"
  completion.flag2<-"END OF LAST FOLLOW-UP PERIOD SET (RANDOMLY) AT maxmaxtime"
  out.obj<-list(completion.flag1=completion.flag1,completion.flag2=completion.flag2,maxmaxtime=maxmaxtime)
  return(out.obj)
}
