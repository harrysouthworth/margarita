egp3Thresh <- function(data, umin=quantile(data, .5), umax=quantile(data, .75),
                       nint = 20, penalty="gaussian", priorParameters=NULL, alpha=.05){
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