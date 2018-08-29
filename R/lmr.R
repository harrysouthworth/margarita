#' Robust regression using MM-estimation
#'
#' Robust regression using MM-estimation with 85\% efficiency for Gaussian data.
#'
#' @param formula A formula describing a linear model.
#' @param data An appropriate data frame.
#' @param weights Not used. This is only here because \code{ggplot2::geom_smooth}
#'   appears to require any custom smoother to take the argument.
#' @param psi The psi function to use. The default is to use bisquare weight functions.
#'   Note that if \code{engine = "rlm"} then psi is set to \code{MASS::psi.bisquare}
#'   (i.e. a function) and if \code{engine = "lmrob"}, psi is set to \code{"bisqare"}
#'   (i.e. a character string). See the help files for \code{rlm} and \code{lmrob}
#'   for further information on available values for psi.
#' @param method The robust fitting method. Defaults to \code{method="MM"}.
#' @param c Tuning parameter to the MM-algorithm. Defaults to \code{c=3.443689} giving 85\% efficiency for Gaussian data.
#' @param engine Character string specifying either 'rlm' in which case \code{MASS::rlm} is used,
#'   or 'lmrob' in which case \code{robustbase::lmrob} is used. In the latter case, a robust
#'   version of R^2 is provided, but the default output produces p-values based on t-distributions
#'   that have no theoretical justification. Bootstrapping would be much better.
#' @param maxit The maximum number of iterations to perform. Defaults to \code{maxit = 40}
#' @param ... Other arguments to be passed to \code{rlm} (like contrasts)
#' @return An object of class 'rlm', fit by the function in the MASS package.
#'         \code{lmr} is just a simple wrapper to \code{rlm}. The returned
#'         object has an additional component, \code{cov}.
#' @return A list containing the same elements as an object of class 'lmr' but with
#'        additional elements containing the covariance matrix of the parameter
#'        estimates ('cov'), the data provided in the call to \code{lmr} ('data'),
#'        the tuning constant for the bisquare loss functions ('c'), and the
#'        residual degrees of freedom.
#' @details The tuning constant for the bisquare function defaults to
#'        \code{c=3.443689} providing 85\% efficiency for Gaussian data.
#'        Maronna et al suggest bisquare weight functions and 85\% efficiency
#'        with MM-estimation in Sections 5.9 and 11.2 of their book. In the setting of
#'        eliminating baseline effects from clinical trial data, the models
#'        considered are fairly simple and these defaults appear to work well.
#'        The value of 3.443689 is 'borrowed' from \code{lmRob} and its support
#'        functions in the 'robust' package. Rounded values for various Gaussian
#'        efficiencies appear in Seciton 2.2 of Maronna et al.
#'
#'        Note that \code{rlm} produces an object
#'        that inherits from \code{lm}, but \code{lmr} does not in order to avoid
#'        \code{lm} methods being used on it. Also, \code{lmr} attaches the
#'        residual degrees of freedom to the object; \code{rlm} deliberately sets
#'        it to NA to avoid erroneous application of \code{lm} methods.
#' @references Venables, W. N. & Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth Edition. Springer, New York. ISBN 0-387-95457-0
#'             Maronna, R. A, Martin, R. D and Yohai, V. J. (2006) Robust Statistics: Theory and Methods, Wiley
#' @keywords models
#' @importFrom MASS rlm
#' @export lmr
lmr <- function(formula, data, weights, psi=NULL, method="MM", c=3.443689, engine="rlm", maxit=40, ...){
    thecall <- match.call()

    if (engine == "rlm") {
      if (is.null(psi)) psi <- MASS::psi.bisquare

      res <- MASS::rlm(formula, data, psi=psi, method=method, c=c, maxit=maxit, ...)

      s <- summary(res)
      res$cov <- s$cov.unscaled * s$sigma^2
      res$data <- data # used by boxplot.rlm
      res$c <- c
      res$call <- thecall

      res$df.residual <- length(res$residuals) - length(res$coefficients)
      res$formula <- formula
    } else if (engine == "lmrob"){
      psi <- if (is.null(psi)) psi <- "bisquare"

      res <- robustbase::lmrob(formula, data, control=robustbase::lmrob.control(tuning.psi = c))
    } else {
      stop("engine should be either 'rlm' or 'lmrob'")
    }

    if (engine == "rlm") {
      class(res) <- c("lmr", "rlm") # drop "lm" because it can lead to errors
    } else {
      class(res) <- c("lmr", "lmrob")
    }
    res
}

#' @method summary lmr
#' @export summary.lmr
summary.lmr <- function(object, ...){
  if ("rlm" %in% class(object)){
    MASS:::summary.rlm(object, ...)
  } else if ("lmrob" %in% class(object)){
    robustbase:::summary.lmrob(object, showAlgo=FALSE, ...)
  } else {
    stop("the object class should include either 'rlm' or 'lmrob'")
  }
}


#' @method predict lmr
#' @export predict.lmr
predict.lmr <- function(object, newdata=NULL, interval="conf", level=0.95, ...){
  suppressWarnings(suppressMessages(stats::predict.lm(object, newdata=newdata, interval=interval, level=level, ...)))
}

#' QQ-plot for residuals from a model using ggplot2
#'
#' @param o A fitted model with a \code{resid} component.
#' @return a \code{ggplot2} plot representing the Gaussian QQ-plot of the residuals
#'           with a fitted line.
#; @keywords hplot
#' @export ggqqplot
ggqqplot <- function(o) {
    y <- quantile(resid(o)[!is.na(resid(o))], c(0.25, 0.75))
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
#' @param hist.scale A scaleing factor for bin widths in the histogram. The default bin width
#'        will be \code{range(resid(data))/hist.scale)}. Defaults to \code{hist.scale=10},
#'        which is smaler than the \code{ggplot} default of 30.
#' @param plot. Logical indicating whether to plot. If FALSE, the function
#'   returns a list of ggplot objects that can be passed to \code{gridExtra::grid.arrange}.
#' @param ... Additional arguments passed to \code{ggplot}. Currently unused.
#' @method ggplot lmr
#' @importFrom gridExtra grid.arrange
#' @export
ggplot.lmr <- function(data=NULL, hist.scale=10, plot.=TRUE, ...){
    d <- data.frame(r = resid(data), f=fitted(data), o=resid(data) + fitted(data))

    qq <- ggqqplot(data)

    bin <- diff(range(d$r, na.rm=TRUE) / hist.scale)
    hist <- ggplot() +
              geom_histogram(data=d, aes(r, ..density..), fill="blue", binwidth=bin) +
              geom_density(data=d, aes(r, ..density..), color="light blue" ) +
              geom_rug(data=d, aes(r), color="orange" ) +
              scale_x_continuous("Residuals")

    fr <- ggplot(d, aes(f, r)) +
            geom_point(color="blue") +
            stat_smooth(method="loess", color="orange",
                        method.args=list(span=2/3, family="symmetric", degree=1)) +
            scale_x_continuous("Fitted values") +
            scale_y_continuous("Residuals")

    fo <-  ggplot(d, aes(f, o)) +
             geom_point(color="blue") +
             stat_smooth(method="loess", color="orange",
                         method.args=list(span=2/3, family="symmetric", degree=1)) +
             scale_x_continuous("Fitted values") +
             scale_y_continuous("Observed values")

    if (plot.){
      gridExtra::grid.arrange(qq, hist, fr, fo)
    } else {
      list(qq, hist, fr, fo)
    }
}

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


