\name{glmhelper2}
\alias{glmhelper2}
\title{Extracts object variables from a regression formula}
\usage{
glmhelper2(input)
}
\arguments{
  \item{input}{a regression formula of the form
  'y~x2+x2...' or 'y~x1*x2+x3...'}
}
\value{
a list which contains the individual objects of the linear
predictor but not the 'intercept' if there is one.
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

