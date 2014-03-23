\name{lmr}
\alias{lmr}
\title{Robust regression using MM-estimation}
\usage{
lmr(fo, data, method = "MM", c = 3.44)
}
\arguments{
  \item{fo}{A formula describing a linear model.}

  \item{data}{An appropriate data frame.}

  \item{method}{The robust fitting method. Defaults to
  \code{method="MM"}.}

  \item{c}{Tuning parameter to the MM-algorithm. Defaults
  to \code{c=3.44} giving 85\% efficiency for Gaussian
  data.}
}
\value{
An object of class 'rlm', fit by the function in the MASS
package.  \code{lmr} is just a simple wrapper to
\code{rlm}. The returned object has an additional
component, \code{cov}.
}
\description{
Robust regression using MM-estimation with 85\% efficiency
for Gaussian data.
}
\keyword{models}

