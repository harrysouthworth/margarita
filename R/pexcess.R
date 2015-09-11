#' Compute the probability of a threshold exceedance across a range of baseline values
#' 
#' @param object An object of class 'margarita'
#' @param M The threshold of interest. \code{M} must be of length 1 and is assumed to be
#'   on the scale of the raw data.
#' @param n The number of values to compute the probability at. Defaults to \code{n=200}. Using
#'   \code{n=50} would be faster, but the resulting graph would be less smooth.
#' @details Creates a sequency of \code{n} equally spaced values between the minimum and
#'   maximum of the baseline values in the robust regression component of the margarita
#'   object, repeats this sequence for each level of the treatment factor, gets the fitted
#'   values from the robust regression for each baseline and treatment, then adds the
#'   GP modelling threshold as obtained from the evmSim component of the margarita object.
#'   The function then uses \code{predict.evmSim} to compute the parameter matrix from
#'   the Markov chains and then computes the probability of exceeding \code{M} across
#'   the range of baseline values and treatment factors.
#' @note WARNING: the code has been tested with treatment as a factor and as linear numeric.
#'   More complicated structures, such as a combo of a factor and a linear effect are unliekly
#'   to work. There is a test function. Take a look at it.
#' @return A \code{data.frame} with columns for baseline, treatment, threshold (linear model
#'   fitted value plus GP threshold) and a summary of the distribution of the probabilities
#'   of exceeding \code{M}. The means and medians ought to be similar. Presumably a plot
#'   of the mean over the range of baseline values is desired. The expected probabilities
#'   could also be used to integrate over the observed baseline distribution to obtain an
#'   estimate of the population probability of exceeding \code{M}.
#' @export pExcessByBaseline
pExcessByBaseline <- function(object, M, n=200){
  if (class(object) != "margarita")
    stop("object should be of class 'margarita'")
  if (missing(M))
    stop("You need to provide the level above which you're interested in exceedances")

  lmod <- object[[1]]
  baseline <- object$baseline
  rawBaseline <- object$rawBaseline
  newdata <- object$newdata
  trans <- object$trans
  invtrans <- object$invtrans

  M <- trans(M)
  
  group <- names(newdata) # Should only be treatment group
  if (is.factor(newdata[, group])) nms <- levels(newdata[, group])
  else nms <- unique(newdata[, group])

  # Get sequence from minimum baseline to maximum
  x <- lmod$x
  x <- x[, colnames(x) == baseline]
  base <- seq(min(x, na.rm=TRUE), max(x, na.rm=TRUE), length.out=n)

  # Repeat it once for each treatment group
  base <- rep(base, length.out=nrow(newdata) * n)

  # Get fitted values given baselines and treatment groups
  newdata <- newdata[rep(1:nrow(newdata), each=n),, drop=FALSE]
  newdata[, rawBaseline] <- invtrans(base)

  # Add the GP modelling thresholds to newdata
  newdata$threshold <- suppressWarnings(predict.lm(lmod, newdata)) + object[[2]]$map$threshold

  # And split it on group
  snd <- split(newdata, newdata[, group])

  # Get the fitted parameters from the evmSim object. Need to get full posterior:
  # can't use ci.fit=TRUE because of relationships between parameters
  param <- predict(object[[2]], type = "lp", newdata=newdata, all=TRUE)$link
  # ^^^^^ XXX XXX XXX POSSIBLE BUG IN predict.evmSim. IF YOU PUT unique=FALSE IN THE ABOVE
  #                   THE OUTPUT IS WRONG! XXX

  # Reorder the elements of param to match snd. Need to be careful because the
  # order returned by predict.evmSim isn't necessarily what we'd expect
  
  # XXX DON'T unlist the line below or you lose a string with zero length
  wh <- lapply(param, function(x) colnames(x)[apply(x, 2, function(x) all(x == 1))])
  wh <- lapply(wh, function(X) gsub(group, "", X))
  names(param)[sapply(wh, function(x) length(x) == 0)] <- nms[1]
  names(param)[sapply(wh, function(x) length(x) > 0)] <- wh[sapply(wh, function(x) length(x) > 0)]
  param <- param[names(snd)]

  # Summarize P(x > u) for every value of baseline for each treatment
  sumfun <- function(x) c(mean(x), quantile(x, prob=c(.05, .25, .5, .75, .95)))
  p <- matrix(0, ncol=6, nrow=n)
  res <- vector("list", length(snd))
  
  theta <- object[[2]]$map$rate # P(x > u) for GP fitting threshold u
  r <- resid(object[[1]]) # Residuals from robust regression
  r <- r[r <= object[[2]]$map$threshold] # Resids beneath GP threshold

  bfun <- function(){
    r <- sample(r, length(fitted(object[[1]])), replace=TRUE) + u
    r[1:round(theta * length(fitted(object[[1]])))] <- u + 1
    mean(r > M)
  }

  for (i in 1:length(snd)){
    for (j in 1:n){
      u <- snd[[i]]$threshold[j]
      if (u <= M){
        p[j, ] <- sumfun(pgpd(M, exp(param[[i]][, 1]), param[[i]][, 2],
                              u=u, lower.tail=FALSE)) * theta
      } else { # u > M
        rp <- replicate(10000, bfun())
        p[j, ] <- sumfun(rp)
      }
    }
    p <- as.data.frame(p)
    names(p) <- c("mean", ".05", ".25", ".5", ".75", ".95")
 
    res[[i]] <- cbind(snd[[i]], p)
  }
  invisible(do.call("rbind", res))
}