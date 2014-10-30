\name{checkLevels}
\alias{checkLevels}
\title{Deals with differing factor levels between studies}
\usage{
checkLevels(form, variables, datasources)
}
\arguments{
  \item{form}{an object of type formula, the input formula}

  \item{variables}{a list of variable names given as calls}

  \item{datasources}{datasources a list of opal object(s)
  obtained after login to opal servers; these objects also
  hold the data assigned to R, as a \code{dataframe}, from
  opal datasources}
}
\value{
an object of type formula, an amended formula if number of
levels differ across studies for any one or more variables
and the initial formula if levels were the same across
studies.
}
\description{
This is an internal function.
}
\details{


The function goes through the variables in a formula and
looks for those of factor type. The number of levels are
then compared across all the studies. If the number of
levels differ, 'dummy' levels are added to the studies with
a lower number of levels and the formula is amended to
include the recoded variable(s).
}
\keyword{internal}

