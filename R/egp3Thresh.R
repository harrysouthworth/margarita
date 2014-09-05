#' Do a test to suggest two possible thresholds above which generalized Pareto modelling can be performed
#' The test is written in the context of modelling residual variation in clinical lab safety data and is unlikely to make sense outside of that context
#' @param data A vector of value to use, typically residuals from a robust linear model having eliminated the baselne effect and possible a treatment effect
#' @param umin The minumum value to consider as a threshold. Because the data are residuals and will be centered near zero, it defaults to zero. Note that negative values don't make sense since the data need to be positive
#' @param umax The maximum value to consider as a threshold. Clinical trials are usually quite small, so having this as a high quantile would often result in too few observations to meaningfully model. Therefore it defaults to being the upper quartile
#' @param nint The number of positions at which to fit the EGP3 model and plot the power parameter with its approximate confidence interval
#' @param penalty The penalty function to use, if any. By default, maximum likelihood is used with no penalty
#' @param priorParameters Parameters to use in the prior distribution if the likelihood is being penalized
#' @param alpha Determines the coverage of the approximate confidence interval on the power parameter
#' @details the EPG3 model of Papastathopoulos and Tawn is fit across a range of thresholds. The threshold at which the confidence interval first contains 1, and at which the point estimate is closest to 1 are returned as suggested thresholds for GPD modelling. The usual diagnostic checks should always be performed.
#' @export
egp3Thresh <- function(data, umin=0, umax=quantile(data, .75),
                       nint = 20, penalty="gaussian", priorParameters=NULL, alpha=.05){
  if (umin < 0) stop("EGP3 models can only be fit to non-negative values, so umin must be >= 0")
  wh <- egp3RangeFit(data=data, umin=umin, umax=umax, nint=nint, penalty=penalty,
                   priorParameters=priorParameters, alpha=alpha)

  shi <- spline(wh$th, wh$hi, n=200)
  slo <- spline(wh$th, wh$lo, n=200)
  spar <- spline(wh$th, wh$par, n=200)

  wh <- shi$y > 1 & slo$y < 1
  res1 <- min(slo$x[wh]) # A confidence limit touches 1
  res2 <- spar$x[abs(spar$y - 1) == min(abs(spar$y - 1))]
  c(res1, res2)
}