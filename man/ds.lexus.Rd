\name{ds.lexus}
\alias{ds.lexus}
\title{Generates an expanded version of a dataset that contains survival data}
\usage{
ds.lexus(x = NULL, intervalWidth = NULL, idCol = NULL, entryCol = NULL,
  exitCol = NULL, statusCol = NULL, variables = NULL, newobj = NULL,
  datasources = NULL)
}
\arguments{
  \item{x}{a character, the name of the table that holds
  the original data, this is the data to be expanded.}

  \item{intervalWidth,}{a numeric vector which gives the
  chosen width of the intervals ('pieces'). This can be one
  value (in which case all the intervals have same width)
  or several different values. If no value(s) are provided
  a single default value is used. That default value is the
  set to be the 1/10th of the mean of the exit time values
  across all the studies.}

  \item{idCol}{a character the name of the column that
  holds the individual IDs of the subjects.}

  \item{entryCol}{a character, the name of the column that
  holds the entry times (i.e. start of follow up). If no
  name is provided the default is to set all the entry
  times to 0 in a column named "STARTTIME". A message is
  then printed to alert the user as this has serious
  consequences if the actual entry times are not 0 for all
  the subjects.}

  \item{exitCol}{a character, the name of the column that
  holds the exit times (i.e. end of follow up).}

  \item{statusCol}{a character, the name of the column that
  holds the 'failure' status of each subject, tells whether
  or not a subject has been censored.}

  \item{variables}{a character vector, the column names of
  the variables (covariates) to include in the final
  expanded table. The input table might have a large number
  of covariates and if only some of those variables are
  relevant for the sought analysis it make sense to only
  include those. By default (i.e. if no variables are
  indicated) all the covariates in the inout table are
  included and this will lengthen the run time of the
  function.}

  \item{newobj}{the name of the output expanded table. By
  default the name is the name of the input table with the
  suffixe "_expanded".}

  \item{datasources}{a list of opal object(s) obtained
  after login to opal servers; these objects also hold the
  data assigned to R, as a \code{data frame}, from opal
  datasources}
}
\value{
a dataframe, an expanded version of the input table.
}
\description{
This function is meant to be used as part of a piecewise
regression analysis.
}
\details{


It splits the survial interval time of subjects into
sub-intervals and reports the failure status of the
subjects at each sub-interval. Each of those sub-interval
is given an id e.g. if the overall interval of a subject is
split into 4 sub-interval, those sub-intervals have ids 1,
2, 3 and 4; so this is basically the count of periods for
each subject. The interval ids are held in a column named
"TIMEID". The entry and exit times in the input table are
used to compute the total survival time. By default all the
covariates in the input table are included in the expanded
output table but it is preferable to indicate the names of
the covariates to be included via the argument 'variables'.
}
\author{
Gaye, A.; Burton, P.
}

