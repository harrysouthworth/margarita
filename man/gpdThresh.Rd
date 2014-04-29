\name{gpdThresh}
\alias{ggplot.gpdRangeFit}
\alias{ggplot.gpdThresh}
\alias{ggplot.mrl}
\alias{gpdThresh}
\title{Produce plots to aid threshold selection for GPD models}
\usage{
\method{ggplot}{gpdThresh}(x, umin = quantile(x, 0.05), umax = quantile(x,
  0.95), nint = 25, priorParameters = NULL, cov = "observed")
}
\arguments{
  \item{x}{A numeric vector.}

  \item{umin}{The minimum value of x to use as a threshold.
  Defaults to \code{umin=quantile(x, .05)}.}

  \item{umax}{The maximum value of x to use as a threshold.
  Defaults to \code{umin=quantile(x, .95)}.}

  \item{nint}{The number of values of x at which to compute
  the parameters. Defaults to \code{nint=25}.}

  \item{data}{An object returned by \code{gpdThresh}.}

  \item{xlab}{Label for the x axis.}

  \item{ylab}{Label fo rthe y axis.}

  \item{main}{Main title.}

  \item{fill}{Fill colour for polynomials representing
  confidence regions.}

  \item{col}{Colour for lines.}

  \item{addNexcesses}{Whether or not to annotate the mean
  residual life plot with the number of threshold excesses.
  Defaults to \code{TRUE}.}

  \item{textsize}{Font size for threshold excesses
  annotation.}

  \item{...}{Additional arguments passed to \code{ggplot}.}
}
\value{
A list of graphical objects created by \code{ggplot}.
}
\description{
Produce plots to aid threshold selection for GPD models
}
\details{
The only argument is the data, so there is no control over
other settings used by \code{mrl}, \code{gpdRangeFit} and
their plot functions.  If more control is needed, use those
functions directly.
}
\keyword{models}

