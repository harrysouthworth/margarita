#' Compute the probability of a threshold exceedance across a range of baseline values
#' 
#' @param object An object of class 'margarita'
#' @param M The threshold of interest. \code{M} must be of length 1
#' @param n The number of values to compute the probability at. Defaults to \code{n=200}. Using
#'   \code{n=50} would be faster, but the resulting graph would be less smooth.
pExcessByBaseline <- function(object, M, n=200){
  if (class(object) != "margarita")
    stop("object should be of class 'margarita'")
  if (missing(M))
    stop("You need to provide the level above which you're interested in exceedances")
  
  lmod <- object[[1]]
  baseline <- object$baseline
  rawBaseline <- object$rawBaseline
  newdata <- object$newdata
  invtrans <- object$invtrans
  
  group <- names(newdata) # Should only be treatment group
  
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
  
  # Reorder the elements of param to match snd
  nms <- unlist(lapply(param, function(X) colnames(X)[3:ncol(X)][colSums(X[, 3:ncol(X)]) > 0]))
  nms <- gsub(group, "", nms)
  nms <- c(nms, names(snd)[!(names(snd) %in% nms)])
  names(param) <- nms
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