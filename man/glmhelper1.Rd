\name{glmhelper1}
\alias{glmhelper1}
\title{Extracts the elements of an expression of type call}
\usage{
glmhelper1(input)
}
\arguments{
  \item{input}{a call of the form 'x1+x2+...' or
  'x1*x2*...', typically this is the rightern side of a
  regression formula (i.e. the linear predictor) with the
  intercept.}
}
\value{
a list which contains the individual elements of the input
expression
}
\description{
This is an internal function required by the client
function \code{ds.glm} to extract variables from a
regression formula
}
\details{

}
\author{
Gaye, A.
}
\keyword{internal}

