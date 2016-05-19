#' @method ggplot egp3RangeFit
#' @export
ggplot.egp3RangeFit <- function(data, xlab = "Threshold", ylab = NULL, main = NULL,
                               fill="cyan", col="orange", reflinecol="blue",
                               addNexcesses = TRUE, textsize=4, ...){
  if (is.null(ylab))
    ylab <- expression(hat(kappa))

  x <- data
  data <- data$data

#  yl <- range(x$hi, x$lo, 1)

  d <- data.frame(th=x$th, par=x$par)
  poly <- data.frame(x=c(x$th, rev(x$th)), y=c(x$lo, rev(x$hi)))

  p <- ggplot(poly, aes(x, y)) +
    geom_polygon(fill=fill, alpha=.5) +
    geom_line(data=d, aes(th, par), color=col) +
    geom_hline(yintercept=1, linetype=2, col=reflinecol) +
    scale_x_continuous(xlab) +
    scale_y_continuous(ylab) +
    theme(axis.title.y=element_text(angle=0)) +
    if (!missing(main)) ggtitle(main)

  if (addNexcesses)
    p <- addExcesses(p, poly$x, poly$y, data=data, u=u, textsize=textsize)
  p
}

#' Produce plots to aid threshold selection for GPD models
#'
#' @aliases gpdThresh, ggplot.gpdThresh ggplot.mrl ggplot.gpdRangeFit
#' @param x A numeric vector.
#' @param umin The minimum value of x to use as a threshold. Defaults to \code{umin=quantile(x, .05)}.
#' @param umax The maximum value of x to use as a threshold. Defaults to \code{umin=quantile(x, .95)}.
#' @param nint The number of values of x at which to compute the parameters. Defaults to \code{nint=25}.
#' @param which Which plots to produce. If \code{which} contains 1, the mean residual life plot is included; if it contains 2, the plots of the generalized Pareto parameters over increasing thresholds are included; if it includes 3, the plot of the EGP3 power parameter is also included. Defaults to including all 3.
#' @param priorParameters A list containing the parameters of the prior distribution, passed through to \code{evm}.
#' @param cov How to compute the covariance, passed through to \code{evm}. Defaults to \code{cov='observed'}, but you might want to use \code{cov='numeric'} if the computations are unstable.
#' @param ... Additional arguments passed to \code{ggplot}.
#' @return A list of graphical objects created by \code{ggplot}.
#' @details The only argument is the data, so there is no control over other
#'    settings used by \code{mrl}, \code{gpdRangeFit} and their plot functions.
#'    If more control is needed, use those functions directly.
#' @keywords models
#' @export gpdThresh
gpdThresh <- function(x, umin=quantile(x, .05),
                         umax=quantile(x, .95),
                         nint=25, which=1:3,
                         priorParameters=NULL, cov="observed"){
  if (class(x) != "numeric"){
    stop(paste("x has class ", class(x), ": should be numeric", sep = ""))
  }

  p <- list()

    wh <- x[x>=umin & x<=umax]
    nint <- min(nint, length(wh))

    if (1 %in% which)
      p <- ggplot(gpdRangeFit(x, umin=umin, umax=umax, nint=nint,
                                   priorParameters=priorParameters, cov=cov))
    if (2 %in% which)
      p <- c(p, mrl=list(ggplot(mrl(x, nint=length(x)))))

    if (3 %in% which)
      p <- c(p, egp3=list(ggplot(egp3RangeFit(x, umin=umin, umax=umax, nint=nint,
                                         priorParameters=priorParameters))))

    oldClass(p) <- 'gpdThresh'
    invisible(p)
}

#' @method ggplot gpdThresh
#' @export
ggplot.gpdThresh <- function(data, ...){
    do.call("grid.arrange", c(data, ncol=2))
    invisible()
}

