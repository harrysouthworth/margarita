\name{egp3Thresh}
\alias{egp3Thresh}

\title{Test to find a threshold above which GPD modelling ought to be ok in the context of clinical trial lab data}
\description{
  Find a threshold above which the EGP3 power parameter is similar to one in a plausible
  range of the data
}
\usage{
egp3Thresh(data, umin=quantile(data, .5), umax=quantile(data, .75), nint = 20, 
            penalty = "gaussian", priorParameters = NULL, alpha=0.05)
}

\arguments{
  \item{data}{ The data vector to be modelled. }
  \item{umin}{ The minimum threshold above which to estimate the parameters. }
  \item{umax}{ The maximum threshold above which to estimate the parameters. }
  \item{nint}{ The number of thresholds at which to perform the estimation. }
  \item{penalty}{The type of penalty to be used in the maximum penalized
               likelihood estimation. Should be either "gaussian" or 
               "none". Defaults to "gaussian".}
  \item{priorParameters}{Parameters to be used for the penalty function.
                 See the help for \code{\link{evm}} for more informaiton.}
  \item{alpha}{100(1 - alpha)\% confidence intervals will be plotted with the 
               point estimates. Defaults to \code{alpha = 0.05}.}
}
\details{
  Internally, \code{egp3Thresh} uses \code{\link{egp3RangeFit}}, so see the help file for
  that function.
  
  In clinical trial applications, we usually deal with the residuals from a linear
  model that eliminates the baseline effect of the variable of interest. Such
  residuals will be centered near 0, so thresholds lower than 0 will be implausible.
  Also, clinical trial datasets tend not to be large, so thresholds greater than
  the upper quartile will unlikely leave enough data to do any meaninful extreme
  value modelling. Therefore, \code{egp3Thresh} considers thresholds between 0
  and the upper quartile and suggests 2 thresholds for consideration: the first
  is the lowest threshold for which the confidence interval for the power parameter
  contains 1; the second is the threshold for which the paramter itself is approximately
  1.
  
  It seems that no such automatic tests for threshold selection are foolproof, so
  standard threshold selection plots and standard diagnostic plots of the fitted
  model should always be examined.
}
\references{
I. Papastathopoulos and J. A. Tawn, Extended generalized Pareto modles for
tail estimation, Journal of Statistical Planning and Inference, 143, 131 -- 143,
2013
}
\author{Harry Southworth}
\seealso{ \code{\link{evm}}, \code{\link{gpdRangeFit}}, \code{\link{mrl}} }
\examples{
liver$r <- resid(rlm(log(ALT.M) ~ log(ALT.B) + as.numeric(dose), data=liver))
for (dose in LETTERS[1:4])
  print(egp3Thresh(liver$r[liver$dose == dose], umin=0))

\keyword{ models }

