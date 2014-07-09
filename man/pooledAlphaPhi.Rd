\name{pooledAlphaPhi}
\alias{pooledAlphaPhi}
\title{Computes the pooled parameters alpha and phi}
\usage{
pooledAlphaPhi(N, npara, M_study, alphaM, sum_p, corStructure = NULL)
}
\arguments{
  \item{N}{a numeric vector, the sample sizes of the
  studies}

  \item{npara}{the number of parameters/columns of the
  design matrix, for a regression model, generated from the
  table of individual level data.}

  \item{M_study}{}

  \item{alphaM}{}

  \item{sum_p}{}

  \item{corStructure}{an integer that set the correlation
  structure: 1 for 'ar1', 2 for 'exchangeable', 3 for
  'independence', 4 for for 'fixed' and 5 for an
  'unstructure' correlation structure.}
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

