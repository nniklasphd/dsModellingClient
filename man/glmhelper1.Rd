\name{glmhelper1}
\alias{glmhelper1}
\title{Checks if the elements in a regression formula are defined and not empty}
\usage{
glmhelper1(input, data, datasources)
}
\arguments{
  \item{input}{a character, a regression formula given as a
  string character}

  \item{data}{a character, the name of an optional data
  frame containing the variables in in the \code{formula}.}

  \item{datasources}{a list of opal object(s) obtained
  after login in to opal servers; these objects hold also
  the data assign to R, as \code{dataframe}, from opal
  datasources.}
}
\value{
an integer 0 if check was passed and 1 if failed
}
\description{
This is an internal function required by the client
function \code{ds.glm} to ensure all the variable in the LP
are defined and not empty, i.e. are not missing at
complete.
}
\details{

}
\author{
Gaye, A.
}
\keyword{internal}

