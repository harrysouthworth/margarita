\name{ggplot.evmOpt}
\alias{ggplot.evmOpt}
\alias{ggplot.evmOpt,}
\alias{ggplot.hist.evmOpt}
\alias{ggplot.ppevm}
\alias{ggplot.qqevm}
\alias{ggplotrl}
\title{Diagnostic plots for an evm object}
\usage{
\method{ggplot}{evmOpt}(data, which = 1:4, main = rep(NULL, 4),
  xlab = rep(NULL, 4), nsim = 1000, alpha = 0.05, ptcol = "blue",
  col = "light blue", fill = "orange", ...)
}
\arguments{
  \item{data}{An object of class \code{evm}.}

  \item{which}{Which plots to produce. Defaults to
  \code{which=1:4}.}

  \item{main}{Main titles. Should have length 4.}

  \item{xlab}{Labels for x-axes.}

  \item{ylab}{Labels for y-axes.}

  \item{nsim}{Number of simulated datasets to create to
  form tolerence regions.}

  \item{alpha}{Used to compute coverage of pointwise
  confidence intervals.}

  \item{ptcol}{Colour for points. Defaults to
  \code{ptcol="blue"}.}

  \item{col}{Colour for lines. Defaults to \code{col="light
  blue"}.}

  \item{fill}{Colour for confidence regions. Defaults to
  \code{fill="orange"}}

  \item{RetPeriodRange}{Range over which to compute return
  levels.}

  \item{...}{Other arguments passed through to underlying
  plot functions.}
}
\description{
Create and display diagnostic plots for an evm object.
}
\details{
The function attempts to arrange the plots nicely. If the
output isn't what was wanted, the function returns the
graphs to the user as a list so that the user can use
\code{grid.arrange} directly.
}
\keyword{hplot}

