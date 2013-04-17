\name{datashield.gee}
\alias{datashield.gee}
\title{Runs a combined GEE analysis of non-pooled data}
\usage{
  datashield.gee(opals, formula, cluster, family, link = 1,
    data = data, corstr = 1, start.betas = NULL,
    user.mats = NULL,
    control = list(maxit = 10, seed = 17632, display = FALSE))
}
\arguments{
  \item{opals}{character strings that represent the URL of
  the servers where the study datasets are stored}

  \item{formula}{an object of class \code{formula} which
  describes the model to be fitted}

  \item{cluster}{a vector that holds the clusters
  (i.e.family IDs)}

  \item{family}{a description of the error distribution
  function to use in the model}

  \item{link}{a description of the link function to be used
  in the model}

  \item{data}{a list of character strings that represent
  the names of the datasets}

  \item{corstr}{a character string that specifies the
  correlation structure to use}

  \item{start.betas}{starting values for the parameters
  estimates}

  \item{user.mats}{user-defined correlation matrices (one
  per dataset) used if the paramter \code{corstr} is set to
  'fixed'}

  \item{control}{a list of parameters for controlling the
  fitting process}
}
\value{
  coefficients estimates of the parameters of the fitted
  model

  stderrs standard errors of the estimates

  alpha correlation parameter, a correlation between any
  two correlated observations

  phi dispersion parameter

  working.corr.matrix the working correlation matrix
}
\description{
  This function allows for the fitting of generalized
  estimating equations to correlated and non-pooled data.
  It enables a parallelized analysis of individual-level
  data sitting on distinct computers/servers by sending
  instructions to each computer requesting non-disclosing
  summary statistics. The sumaries are then combined to
  estimate the parameters of the model; these parameters
  are the same as those obtained if the data were
  'physically' pooled.
}
\details{
  This function requires to log into the opal servers where
  the data are stored (one dataset per study). Set the
  parameter \code{corstr} to set to 1, 2, 3 , 4 or 5 for
  respectively 'ar1', 'exchangeable', 'independence',
  'fixed', or 'unstructured' correlation structures. Set
  the parameter \code{link} to 1, 2, 3, 4, or 5 for
  respectively the 'logit', 'identity', 'inverse', 'log' or
  'probit' link functions. The \code{binomial family}
  accepts the links 'logit', 'probit' and 'log'; the
  \code{gaussian family} accepts the links 'identity',
  'log' and 'inverse' and the \code{Gamma family} acepts
  the links 'inverse', 'identity' and 'log'.
}
\examples{
\dontrun{
# put example here
}
}
\references{
  Zeger, S. L., Liang, K. and Albert, P. S. (1988). Models
  for longitudinal data: a generalized estimating equation
  approach. Biometrics, 44, p1049-1060.
}

