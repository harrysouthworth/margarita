\name{ggplotrl}
\alias{ggplotrl}
\title{Diagnostic plots for an evm object}
\usage{
ggplotrl(data, alpha = 0.05, xlab, ylab, main, ptcol = "blue",
  col = "light blue", fill = "orange", RetPeriodRange = NULL, ...)
}
\arguments{
  \item{data}{An object of class 'evm'.}

  \item{alpha}{Used to compute coverage of pointwise
  confidence intervals.}

  \item{xlab}{Labels for x-axes.}

  \item{ylab}{Labels for y-axes.}

  \item{main}{Main titles. Should have length 4.}

  \item{ptcol}{Colour for points. Defaults to
  \code{ptcol="blue"}.}

  \item{col}{Colour for lines. Defaults to \code{col="light
  blue"}.}

  \item{fill}{Colour for confidence regions. Defaults to
  \code{fill="orange"}}

  \item{RetPeriodRange}{Range over which to compute return
  levels.}

  \item{...}{Additional arguments to \code{ggplot}.}
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

