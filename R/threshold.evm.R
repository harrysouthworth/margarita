ggplot.mrl <- function(data, xlab = "Threshold", ylab = "Mean excess", main=NULL,
                       fill="light blue", col="blue",
                       addNexcesses=TRUE, textsize=4, ...){
    x <- data
    data <- x$data
    x <- x$mrl
    d <- data.frame(th = x[, "threshold"],
                    mrl = x[, "MRL"],
                    xl = x[, "lo"],
                    xu = x[, "hi"])

    k <- !is.na(d$xl)
    poly <- data.frame(x=c(d$th, rev(d$th)), y=c(d$xl, rev(d$xu)))
    poly <- poly[c(k, rev(k)), ]

    p <- ggplot(poly, aes(x, y)) +
             geom_polygon(fill=fill) +
             geom_line(data=d, aes(th, mrl), color=col) +
             scale_x_continuous(xlab) +
             scale_y_continuous(ylab) +
             ggtitle(main)

    if (addNexcesses)
        p <- addExcesses(p, poly$x, poly$y, data=data, u=d$th,
                         textsize=textsize)

    p
}

ggplot.gpdRangeFit <- function(data, xlab = "Threshold", ylab = NULL, main = NULL,
                               fill="orange", col="blue",
                               addNexcesses = TRUE, textsize=4, ...){
    if (missing(ylab)) {
        ylab <- c(expression(hat(phi)[m]), expression(hat(xi)))
    }
    else if (length(ylab) != 2) {
        stop("length of ylab should be 2")
    }
    if (!missing(main) && length(main) != 2) {
        stop("length of main should be 2")
    }
    
    x <- data
    data <- data$data

    p <- vector("list", 2)
    
    for (i in 1:2) {
        yl <- range(x$hi[, i], x$lo[, i])
        
        d <- data.frame(th=x$th, par=x$par[, i])
        poly <- data.frame(x=c(x$th, rev(x$th)), y=c(x$lo[, i], rev(x$hi[, i])))
        
        p[[i]] <- ggplot(poly, aes(x, y)) +
                    geom_polygon(fill=fill) +
                    geom_line(data=d, aes(th, par), color=col) +
                    scale_x_continuous(xlab) +
                    scale_y_continuous(ylab[i]) +
                    theme(axis.title.y=element_text(angle=0)) +
                    if (!missing(main)) ggtitle(main[i])

        if (addNexcesses)
            p[[i]] <- addExcesses(p[[i]], poly$x, poly$y, data=data, u=u, textsize=textsize)
    } # Close for
    p
}

#' Produce plots to aid threshold selection for GPD models
#'
#' @aliases ggplot.gpdThresh ggplot.mrl ggplot.gpdRangeFit
#' @method ggplot gpdThresh
#' @method ggplot mrl
#' @method ggplot gpdRangeFit
#' @param x A numeric vector.
#' @param umin The minimum value of x to use as a threshold. Defaults to \code{umin=quantile(x, .05)}.
#' @param umax The maximum value of x to use as a threshold. Defaults to \code{umin=quantile(x, .95)}.
#' @param nint The number of values of x at which to compute the parameters. Defaults to \code{nint=25}.
#' @param data An object returned by \code{gpdThresh}.
#' @param xlab Label for the x axis.
#' @param ylab Label fo rthe y axis.
#' @param main Main title.
#' @param fill Fill colour for polynomials representing confidence regions.
#' @param col Colour for lines.
#' @param addNexcesses Whether or not to annotate the mean residual life plot with
#'        the number of threshold excesses. Defaults to \code{TRUE}.
#' @param textsize Font size for threshold excesses annotation.
#' @param ... Additional arguments passed to \code{ggplot}.
#' @return A list of graphical objects created by \code{ggplot}.
#' @details The only argument is the data, so there is no control over other
#'    settings used by \code{mrl}, \code{gpdRangeFit} and their plot functions.
#'    If more control is needed, use those functions directly.
#' @keywords models
#' @export gpdThresh
gpdThresh <- function(x, umin=quantile(x, .05),
                         umax=quantile(x, .95),
                         nint=25){
    m <- ggplot(mrl(x, nint=length(x)))
    g <- ggplot(gpdRangeFit(x, umin=umin, umax=umax, nint=nint))
    
    res <- list(g[[1]], g[[2]], m)
    oldClass(res) <- 'gpdThresh'
    invisible(res)
}

#' @method ggplot gpdThresh
#' @export
ggplot.gpdThresh <- function(data, ...){
    blankPanel <- grid.rect(gp=gpar(col="white"))    
    grid.arrange(data[[1]], data[[2]], data[[3]], blankPanel)
    invisible()
}

