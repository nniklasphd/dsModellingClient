\name{ds.glm}
\alias{ds.glm}
\title{Fits a Generalized Linear Model (GLM) model}
\usage{
ds.glm(formula = NULL, family = NULL, maxit = 15, CI = 0.95,
  viewIter = FALSE, datasources = NULL)
}
\arguments{
  \item{formula}{an object of class \code{formula} which
  describes the model to be fitted}

  \item{family}{a description of the error distribution
  function to use in the model}

  \item{maxit}{the number of iterations of IWLS used
  instructions to each computer requesting non-disclosing
  summary statistics. The summaries are then combined to
  estimate the parameters of the model; these parameters
  are the same as those obtained if the data were
  'physically' pooled.}

  \item{CI}{a numeric, the confidence interval.}

  \item{viewIter}{a boolean, tells whether the results of
  the intermediate iterations should be printed on screen
  or not. Default is FALSE (i.e. only final results are
  shown).}

  \item{datasources}{a list of opal object(s) obtained
  after login to opal servers; these objects also hold the
  data assigned to R, as a \code{dataframe}, from opal
  datasources.}
}
\value{
coefficients a named vector of coefficients

residuals the 'working' residuals, that is the residuals in
the final iteration of the IWLS fit.

fitted.values the fitted mean values, obtained by
transforming the linear predictors by the inverse of the
link function.

rank the numeric rank of the fitted linear model.

family the \code{family} object used.

linear.predictors the linear fit on link scale.

aic A version of Akaike's An Information Criterion, which
tells how well the model fits
}
\description{
A function that fits generalized linear models
}
\details{


It enables a parallelized analysis of individual-level data
sitting on distinct servers by sending commands to each
data computer to fit a regression model. The estimates
returned are then combined and updated coefficients
estimate sent back for a new fit. This iterative process
goes on until convergence is achieved.
}
\examples{
{

# load the file that contains the login details
data(logindata)

# login and assign some variables to R
myvar <- list('DIS_DIAB','PM_BMI_CONTINUOUS','LAB_HDL', 'GENDER')
opals <- datashield.login(logins=logindata,assign=TRUE,variables=myvar)

# Example 1: run a GLM without interaction (e.g. diabetes prediction using BMI and HDL levels and GENDER)
mod <- ds.glm(formula=D$DIS_DIAB~D$PM_BMI_CONTINUOUS+D$LAB_HDL+D$GENDER,family='binomial')

# Example 2: run the above GLM model with an intercept (eg. intercept = 1)
 mod <- ds.glm(formula=D$DIS_DIAB~1+D$PM_BMI_CONTINUOUS+D$LAB_HDL+D$GENDER,family='binomial')

# Example 3: run the above GLM with interaction HDL and GENDER
 mod <- ds.glm(formula=D$DIS_DIAB~D$PM_BMI_CONTINUOUS+D$LAB_HDL*D$GENDER,family='binomial')

# Example 4: now run the same GLM but with interaction between BMI and HDL
 mod <- ds.glm(formula=D$DIS_DIAB~D$PM_BMI_CONTINUOUS*D$LAB_HDL+D$GENDER,family='binomial')
}
}
\author{
Burton,P;Gaye,A;Laflamme,P
}
\references{
Jones EM, 'DataSHIELD-shared individual-level analysis
without sharing the data: a biostatistical perspective',
Norsk epidemiologi - Norwegian journal of epidemiology
2012;21(2): 231-9.
}

