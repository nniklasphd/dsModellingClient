\name{getBeta}
\alias{getBeta}
\title{Produces beta values and standard errors}
\usage{
getBeta(score, infoMatrix, Jmatrix, startCoeff)
}
\arguments{
  \item{score}{a numeric vector, a score vector}

  \item{infoMatrix}{a matrix, the information matrix}

  \item{Jmatrix}{a matrix, the j matrix}

  \item{startCoeff}{a numeric vector, the starting values
  for the beta coefficients.}
}
\value{
a list which contains the individual elements of the input
expression
}
\description{
This is an internal function required by the client
function \code{ds.gee}.
}
\details{

}
\author{
Gaye, A.; Jones EM.
}
\keyword{internal}

