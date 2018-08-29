#' Shiftplot
#'
#' Shiftplots using ggplot with optional trellising
#' @param data A \code{data.frame}.
#' @param aes An object created with \code{aes}.
#' @param by A formula for passing to \code{facet_wrap}.
#' @param ncol The number of columns to lay the plots out in.
#' @param trans Character string defining the axis transformations. Defaults
#'              to \code{trans="log10"}. Use \code{trans="identity"} if no
#'              transformation is to take place.
#' @param xlab The horizontal axis label.
#' @param ylab The vertical axis label.
#' @param main The main title.
#' @param jitter.amount The amount by which to jitter the data. Defaults to no jittering.
#' @param alpha Value of alpha for transparency. Defaults to no transparency.
#' @param shape The plot symbol to use.
#' @param ptcol The colour for points on the plots.
#' @param linecol The colour for the reference line.
#' @param ... Additional arguments to \code{ggplot}. Currently unused.
#' @param theme An object produced by a call to \code{theme} to be added to the ggplot.
#' @export shiftplot
shiftplot <- function(data, aes, by=NULL, ncol=NULL, trans="identity",
                      xlab="Baseline", ylab="Maximum", main=NULL,
                      jitter.amount=0, alpha=1, shape=16,
                      ptcol="blue", linecol="orange", theme=NULL, ...){
  d <- data[, c(as.character(aes$x)[2], as.character(aes$y)[2])]
  limits <- range(d, na.rm=TRUE) + c(-jitter.amount, jitter.amount)

  ppts <-
    if (! "colour" %in% names(aes))
      geom_jitter(color=ptcol, alpha=alpha, shape=shape,
                  position=position_jitter(width=jitter.amount, height=jitter.amount))
  else
    geom_jitter(alpha=alpha, shape=shape,
                position=position_jitter(width=jitter.amount, height=jitter.amount))

  p <- ggplot(data, aes) + ppts + theme +
    geom_abline(color=linecol, intercept=0, slope=1) +
    scale_x_continuous(xlab, limit=limits, trans=trans) +
    scale_y_continuous(ylab, limit=limits, trans=trans) +
    coord_fixed() +
    ggtitle(main) +
    if (!is.null(by)) facet_wrap(by, ncol=ncol)
  p
}
