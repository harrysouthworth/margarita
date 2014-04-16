#' Robust regression using MM-estimation
#'
#' Robust regression using MM-estimation with 85\% efficiency for Gaussian data.
#'
#' @param fo A formula describing a linear model.
#' @param data An appropriate data frame.
#' @param method The robust fitting method. Defaults to \code{method="MM"}.
#' @param c Tuning parameter to the MM-algorithm. Defaults to \code{c=3.44} giving 85\% efficiency for Gaussian data.
#' @return An object of class 'rlm', fit by the function in the MASS package.
#'         \code{lmr} is just a simple wrapper to \code{rlm}. The returned
#'         object has an additional component, \code{cov}.
#' @keywords models
#' @importFrom MASS rlm
#' @export lmr
lmr <- function(fo, data, method="MM", c=3.44){
    res <- rlm(fo, data, method=method, c=c)
    res$call$formula <- fo
    s <- summary(res)
    res$cov <- s$cov.unscaled * s$stddev^2
    res$data <- data # used by boxplot.rlm
    # Can't add simulated coefs here because we don't know how many
    # to simulate yet.
    res
}

#' QQ-plot for residuals from a model using ggplot2
#'
#' @param o A fitted model with a \code{resid} component.
#' @return a \code{ggplot2} plot representing the Gaussian QQ-plot of the residuals
#'           with a fitted line.
#; @keywords hplot
#' @export ggqqplot
ggqqplot <- function(o) {
    y <- quantile(o$resid[!is.na(o$resid)], c(0.25, 0.75))
    x <- qnorm(c(0.25, 0.75))
    slope <- diff(y)/diff(x)
    int <- y[1L] - slope * x[1L]

    d <- data.frame(r=resid(o))

    p <- ggplot(d, aes(sample = r)) +
           stat_qq(alpha = 0.5, color="blue") +
           geom_abline(slope = slope, intercept = int, color="orange")

    return(p)
}

#' ggplot for rlm results
#' @param data an object of class \code{rlm}.
#' @param ... Additional arguments passed to \code{ggplot}. Currently unused.
#' @method ggplot rlm
#' @importFrom gridExtra grid.arrange
#' @export
ggplot.rlm <- function(data=NULL, ...){
    d <- data.frame(r = resid(data), f=fitted(data), o=resid(data) + fitted(data))
    
    qq <- ggqqplot(data)
    
    hist <- ggplot() +
              geom_histogram(data=d, aes(r, ..density..), fill="blue" ) +
              geom_density(data=d, aes(r, ..density..), color="light blue" ) +
              geom_rug(data=d, aes(r), color="orange" ) +
              scale_x_continuous("Residuals")

    fr <- ggplot(d, aes(f, r)) +
            geom_point(color="blue") +
            geom_smooth(method="loess", color="orange") +
            scale_x_continuous("Fitted values") +
            scale_y_continuous("Residuals")

    fo <-  ggplot(d, aes(f, o)) +
             geom_point(color="blue") +
             geom_smooth(method="loess", color="orange") +
             scale_x_continuous("Fitted values") +
             scale_y_continuous("Observed values")

    grid.arrange(qq, hist, fr, fo)
}

#' Boxplots of scaled residuals, by a factor in the model.
#'
#' Boxplots of scaled residuals, split by a factor in the model.
#' @param x An object of class 'rlm'.
#' @param by A character string giving the name of a factor in \code{x}.
#' @param jitter.width The amount of jittering to do. Defaults to \code{jitter.width = 0.1}.
#' @param xlab Horizontal axis label. Defaults to \code{xlab=""}.
#' @param ylab Vertical axis label. Defaults to \code{ylab="Scaled residuals"}.
#' @param main Main title. Defaults to blank.
#' @param ptcol The colour of the residuals on the plot. Defaults to \code{"orange"}.
#' @param boxcol The colour of the outline of the boxplots. Defaults to \code{"blue"}.
#' @param boxfill The fill colour of the boxplots. Defaults to \code{"light blue"}.
#' @param ... Additional arguments to \code{boxplot}. Not currently used.
#' @method boxplot rlm
#' @export
#' @importFrom graphics boxplot
boxplot.rlm <- function(x, by, jitter.width=.1, xlab="", ylab="Scaled residuals",
                        main="", ptcol="orange", boxcol="blue", boxfill="light blue", ...){
    trt <- x$data[, by]
    data <- data.frame(sr = resid(x) / x$s, trt=as.factor(trt))
    p <- ggplot(data, aes(trt, sr)) +
           geom_jitter(color=ptcol, position=position_jitter(width=jitter.width)) +
           stat_boxplot(color=boxcol, fill=boxfill, alpha=.5, outlier.shape="") +
           scale_x_discrete(xlab) +
           scale_y_continuous(ylab) +
           ggtitle(main)
    p
}

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
#' @export shiftplot
shiftplot <- function(data, aes, by=NULL, ncol=NULL, trans="identity",
                      xlab="Baseline", ylab="Maximum", main=NULL,
                      jitter.amount=0, alpha=1, shape=16,
                      ptcol="blue", linecol="orange", ...){
  d <- data[, c(as.character(aes$x), as.character(aes$y))]
  limits <- range(d, na.rm=TRUE) + c(-jitter.amount, jitter.amount)

  p <- ggplot(data, aes) +
         geom_jitter(color=ptcol, alpha=alpha, shape=shape,
                     position=position_jitter(width=jitter.amount, height=jitter.amount)) +
         geom_abline(color=linecol, intercept=0, slope=1) +
         scale_x_continuous(xlab, limit=limits, trans=trans) +
         scale_y_continuous(ylab, limit=limits, trans=trans) +
         coord_fixed() +
         ggtitle(main) +
         if (!is.null(by)) facet_wrap(by, ncol=ncol)

  p
}

