\name{shiftplot}
\alias{shiftplot}
\title{Shiftplot}
\usage{
shiftplot(data, aes, by = NULL, ncol = NULL, trans = "log10",
  xlab = "Baseline", ylab = "Maximum", main = NULL, jitter.amount = 0,
  alpha = 1, shape = 16, ptcol = "blue", linecol = "orange", ...)
}
\arguments{
  \item{data}{A \code{data.frame}.}

  \item{aes}{An object created with \code{aes}.}

  \item{by}{A formula for passing to \code{facet_wrap}.}

  \item{ncol}{The number of columns to lay the plots out
  in.}

  \item{trans}{Character string defining the axis
  transformations. Defaults to \code{trans="log10"}. Use
  \code{trans="identity"} if no transformation is to take
  place.}

  \item{xlab}{The horizontal axis label.}

  \item{ylab}{The vertical axis label.}

  \item{main}{The main title.}

  \item{jitter.amount}{The amount by which to jitter the
  data. Defaults to no jittering.}

  \item{alpha}{Value of alpha for transparency. Defaults to
  no transparency.}

  \item{shape}{The plot symbol to use.}

  \item{ptcol}{The colour for points on the plots.}

  \item{linecol}{The colour for the reference line.}

  \item{...}{Additional arguments to \code{ggplot}.
  Currently unused.}
}
\description{
Shiftplots using ggplot with optional trellising
}

