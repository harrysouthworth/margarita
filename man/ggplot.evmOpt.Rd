\name{ggplot.evmOpt}
\alias{ggplot.evmOpt}
\alias{ggplot.hist.evmOpt}
\alias{ggplot.ppevm}
\alias{ggplot.qqevm}
\alias{ggplotrl}
\title{Diagnostic plots for an \code{evm} object}
\usage{

\method{ggplot}{evmOpt}(data, which = 1:4, main = rep(NULL, 4),
  xlab = rep(NULL, 4), nsim = 1000, alpha = 0.05, ptcol = "blue",
  col = "light blue", fill = "orange", ...)
}
\arguments{
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

  \item{...}{Other arguments passed through to underlying
  plot functions.}
}

