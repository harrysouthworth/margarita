% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{stack.list}
\alias{stack.list}
\title{Stack a list of data.frames or matrices}
\usage{
\method{stack}{list}(x, ...)
}
\arguments{
\item{x}{A list containing data.frames or matrices with the same number of columns}

\item{...}{Additional arguments, currently unused.}
}
\description{
Stack a list of data.frames or matrices that have the same number of columns
}
\details{
A rudimentary check is performed to see if the objects in \code{x} have the
         same number of columns. If so, \code{rbind} is used to stack them.
}
