\name{glmhelper3}
\alias{glmhelper3}
\title{Extracts the elements of an expression of type call}
\usage{
glmhelper3(dtname, form, datasources)
}
\arguments{
  \item{dtname}{the name if any of the data frame that hold
  the variables in the regression formula.}

  \item{form}{an object of type formula which describes the
  model to be fitted.}

  \item{datasources}{a list of opal object(s) obtained
  after login to opal servers; these objects also hold the
  data assigned to R, as a \code{dataframe}, from opal
  datasources.}
}
\value{
if no error message is thrown a list with the names of the
variables as call objects is returned.
}
\description{
This is an internal function required by the client
function \code{ds.glm} to carry checks.
}
\details{


the functions checks if all the variables in the lp formula
exist (i.e. are defined) on the server site and if any is
empty (contains NAs only). If any of the checks fails the
process stops and an error message is thrown.
}
\author{
Gaye, A.
}
\keyword{internal}

