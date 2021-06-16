#' Boxplots of scaled residuals, by a factor in the model.
#'
#' Boxplots of scaled residuals, split by a factor in the model.
#' @param x An object of class 'rlm'.
#' @param by A character string giving the name of a factor in \code{x}.
#' @param jitter.width The amount of jittering to do. Defaults to \code{jitter.width = 0.1}.
#' @param box Whether to do the boxplot. If the amount of data is small, it might be better
#'        to just plot the jittered points. Defaults to \code{box=TRUE}.
#' @param xlab Horizontal axis label. Defaults to \code{xlab=""}.
#' @param ylab Vertical axis label. Defaults to \code{ylab="Scaled residuals"}.
#' @param main Main title. Defaults to blank.
#' @param size The size of the plot symbols. Defaults to 3.
#' @param ptcol The colour of the residuals on the plot. Defaults to \code{"orange"}.
#' @param boxcol The colour of the outline of the boxplots. Defaults to \code{"blue"}.
#' @param boxfill The fill colour of the boxplots. Defaults to \code{"light blue"}.
#' @param ... Additional arguments to \code{boxplot}. Not currently used.
#' @method boxplot lmr
#' @export boxplot.lmr
#' @importFrom graphics boxplot
boxplot.lmr <- function(x, by, jitter.width=.1, box=TRUE, xlab="", ylab="Scaled residuals",
                        main="", size=3, ptcol="orange", boxcol="blue", boxfill="light blue", ...){
    trt <- x$data[, by]
    data <- data.frame(sr = resid(x) / x$s, trt=as.factor(trt))
    if (box)
      bp <- stat_boxplot(color=boxcol, fill=boxfill, alpha=.5, outlier.shape="")
    else
      bp <- NULL
    p <- ggplot(data, aes(trt, sr)) +
           geom_jitter(color=ptcol, size=size, alpha=.7, position=position_jitter(width=jitter.width)) +
           scale_x_discrete(xlab) +
           scale_y_continuous(ylab) +
           ggtitle(main) +
           bp
    p
}

