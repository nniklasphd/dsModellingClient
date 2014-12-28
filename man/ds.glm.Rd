\name{ds.glm}
\alias{ds.glm}
\title{Runs a combined GLM analysis of non-pooled data}
\usage{
ds.glm(formula = NULL, data = NULL, family = NULL, offset = NULL,
  weights = NULL, checks = TRUE, maxit = 15, CI = 0.95,
  viewIter = FALSE, datasources = NULL)
}
\arguments{
  \item{formula}{a character, a formula which describes the
  model to be fitted}

  \item{family}{a description of the error distribution
  function to use in the model}

  \item{offset}{a character, null or a numeric vector that
  can be used to specify an a priori known component to be
  included in the linear predictor during fitting.}

  \item{weights}{a character, the name of an optional
  vector of 'prior weights' to be used in the fitting
  process. Should be NULL or a numeric vector.}

  \item{data}{a character, the name of an optional data
  frame containing the variables in in the \code{formula}.
  The process stops if a non existing data frame is
  indicated.}

  \item{checks}{a boolean, if TRUE (default) checks that
  takes 1-3min are carried out to verify that the variables
  in the model are defined (exist) on the server site and
  that they have the correct characteristics required to
  fit a GLM. If FALSE (not recommended if you are not an
  experienced user) no checks are carried except some very
  basic ones and eventual error messages might not give
  clear indications about the cause(s) of the error.}

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
A function fit generalized linear models
}
\details{


It enables a parallelized analysis of individual-level data
sitting on distinct servers by sending
}
\examples{
{

# load the file that contains the login details
data(glmLoginData)

# login and assign some variables to R
myvar <- list("DIS_DIAB","PM_BMI_CONTINUOUS","LAB_HDL", "GENDER")
opals <- datashield.login(logins=glmLoginData,assign=TRUE,variables=myvar)

# Example 1: fit a GLM without interaction (e.g. diabetes prediction using BMI and HDL levels and GENDER)
mod <- ds.glm(formula='DIS_DIAB~PM_BMI_CONTINUOUS+LAB_HDL+GENDER', data='D', family='binomial')

# Example 2: fit the above GLM model with an intercept (eg. intercept = 1)
 mod <- ds.glm(formula='DIS_DIAB~1+PM_BMI_CONTINUOUS+LAB_HDL+GENDER', data='D', family='binomial')

# Example 3: fit the above GLM with interaction HDL and GENDER
 mod <- ds.glm(formula='DIS_DIAB~PM_BMI_CONTINUOUS+LAB_HDL*GENDER', data='D', family='binomial')

# Example 4: fit a GLM with an offset parameter: here use continuous BMI as offset parameter
 mod <- ds.glm(formula='LAB_HDL~GENDER', data="D", family='gaussian', offset="PM_BMI_CONTINUOUS")

# Example 5: fit a GLM with an offset parameter: here use continuous BMI values as weigths
 mod <- ds.glm(formula='LAB_HDL~GENDER', data="D", family='gaussian', weights="PM_BMI_CONTINUOUS")

# clear the Datashield R sessions and logout
datashield.logout(opals)
}
}
\author{
Burton,P;Gaye,A;Laflamme,P
}
\seealso{
\code{ds.lexus} for 'ds.glm' in survival analysis using
piecewise exponential regression
}

