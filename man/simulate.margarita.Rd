\name{simulate.margarita}
\alias{simulate.margarita}
\alias{simulate.margarita.prob}
\alias{simulate.margarita.rl}
\title{Simulate return levels and probabilities of threshold exceedences from a robust linear model and an extreme value model}
\usage{
\method{simulate}{margarita}(object, nsim = 1, seed = NULL, type = "rl",
  M = NULL, scale = "raw", ...)
}
\arguments{
  \item{object}{An object of class 'margarita'}

  \item{nsim}{Unused argument.}

  \item{seed}{Unused argument.}

  \item{type}{What type of prediction is required: either
  'rl' or 'prob'.  Defaults to \code{type = "rl"}.}

  \item{M}{The return level to be predicted. Defaults to
  \code{M=1000}. If \code{type="prob"}, M should be a
  vector containing the thresholds whose probabilities of
  exceedance the user is interested in, on the scale of the
  raw data.}

  \item{scale}{The type of prediction being made. Valid
  values are 'raw', 'proportional' and 'difference'. Due to
  how the calculations are performed, \code{scale} needs to
  be specified in the call to \code{simulate} if threshold
  exceedence probabilities are being predicted, or the call
  to \code{summary.simulate} if return levels are being
  predicted. If \code{scale='raw'} it is assumed that the
  values of \code{M} are on the same scale as the response
  variable in the GPD model. If
  \code{scale='proportional'}, the values of \code{M} are
  taken to be fold-changes (e.g. 2, 5, 10 times baseline
  for return levels, or P(yMax > 2, 5, 10 times yBase)). If
  \code{scale='difference'}, absolute changes from baseline
  are assumed.}

  \item{...}{Other agruments passed to \code{simulate}.
  Currently unused.}
}
\description{
Simulate return levels and probabilities of threshold
exceedences from a robust linear model and an extreme value
model. The procedure is specific to the situation of
extreme value modelling of residuals when predictions on
the scale of the original data are required.
}
\keyword{datagen}

