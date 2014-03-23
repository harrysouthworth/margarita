\name{boxplot.rlm}
\alias{boxplot.rlm}
\title{Boxplots of scaled residuals, by a factor in the model.}
\usage{
\method{boxplot}{rlm}(x, by, jitter.width = 0.1, xlab = "",
  ylab = "Scaled residuals", main = "", ptcol = "orange",
  boxcol = "blue", boxfill = "light blue", ...)
}
\arguments{
  \item{x}{An object of class 'rlm'.}

  \item{by}{A character string giving the name of a factor
  in \code{x}.}

  \item{jitter.width}{The amount of jittering to do.
  Defaults to \code{jitter.width = 0.1}.}

  \item{xlab}{Horizontal axis label. Defaults to
  \code{xlab=""}.}

  \item{ylab}{Vertical axis label. Defaults to
  \code{ylab="Scaled residuals"}.}

  \item{main}{Main title. Defaults to blank.}

  \item{ptcol}{The colour of the residuals on the plot.
  Defaults to \code{"orange"}.}

  \item{boxcol}{The colour of the outline of the boxplots.
  Defaults to \code{"blue"}.}

  \item{boxfill}{The fill colour of the boxplots. Defaults
  to \code{"light blue"}.}

  \item{...}{Additional arguments to \code{boxplot}. Not
  currently used.}
}
\description{
Boxplots of scaled residuals, split by a factor in the
model.
}

