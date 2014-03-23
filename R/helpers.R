#' Get data for plotting segments representing interval estimates.
#' @param data A data frame with either 5 or 3 columns representing the point
#         estimates and interval estimates, and a column called 'group'.
getSegmentData <- function(data){
    if (ncol(data) > 5){
        seg1 <- data.frame(lo=data[, 1], hi=data[, 5], group=data$group, stringsAsFactors=FALSE)
        seg2 <- data.frame(lo=data[, 2], hi=data[, 4], group=data$group, stringsAsFactors=FALSE)
        if (any(seg1[, "lo"] > seg1[, "hi"]) | any(seg2[, "lo"] > seg2[, "hi"])){
            stop("CI segments with lo > hi")
        }
    }
    else if (ncol(data) > 3) {
        seg1 <- data.frame(lo=data[, 1], hi=data[, 3], group=data$group, stringsAsFactors=FALSE)
        seg2 <- NULL
        if (any(seg1[, "lo"] > seg1[, "hi"])){
            stop("CI segments with lo > hi")
        }
    }
    else {
        stop("segment data has 3 or fewer columns")
    }
    list(seg1, seg2)
}

#' Get quantiles representing confidence limits from a given test size.
#' @param alpha The size of the test. Values {alpha/2, 1 - alpha/2} will be
#'        returned. \code{alpha} can have length either 1 or 2.
getCIquantiles <- function(alpha){
    # Get quantiles for CI plus the median, in order
    if (length(alpha) > 2){ stop("alpha can have length 1 or 2") }
    if (any(alpha >= 1) | any(alpha <= 0)){ stop("alpha should be on the interval (0, 1)") }
    alpha <- sort(alpha)
    c(alpha/2, 0.50, 1-rev(alpha)/2)
}

#' Return period for GPD
#'
#' Compute the return period for a genralized Pareto distribution
#' @param X An indexing parameter used in a call to \code{sapply} or \code{lapply}.
#' @param xm The threshold above which we wish to estimated the exceedance rate.
#' @param u The threshold.
#' @param phi The log-scale parameter for the GPD model.
#' @param xi The shape parameter for the GPD model.
#' @param p The rate of threshold exceedance.
#' @param r The data to which the original model was fit.
margarita.rp <- function(X, xm, u, phi, xi, p, r) {
    xm <- xm[, X]
    res <- p * (1 + xi/exp(phi) * (xm - u))^(-1/xi)
    wh <- u > xm
    if (any(wh)){
         res[wh] <- sapply(1:sum(wh),
                              function(i, x, r, m, p) mean((r + x[i] - quantile(r,1-p)) > m[i]),
                              x=u[wh], r=r, m=xm[wh], p=p)
    }
    # Set P(x > upper limit) = 0 when upper limit exists (i.e. for any xi < 0).
    res[xi < 0 & xm > u - exp(phi)/xi] <- 0
    res
}

#' Get probabilities of threshold exceedance for a GPD model.
#'
#' Get probabilities of threshold exceedance for a genralized Pareto model.
#' @param X An index parameter used in a call to \code{sapply} or \code{lapply}.
#' @param par A numberic vector containin the parameters for the GPD model.
#' @param u The threshold.
#' @param p The rate of threshold exceedance.
#' @param r The original data to which the model was fit.
#' @param m The threshold above which we wish to estimated the exceedance rate.
margarita.getProbs <- function(X, par, u, p, r, m) {
    u <- u[[X]]
    phi <- par[[X]][, 1]
    xi <- par[[X]][, 2]
    lapply(1:ncol(m), margarita.rp, u=u, phi=phi, xi=xi, p=p, r=r, xm=m)
}

#' Construct a matrix of thresholds whose rate of exceedance we wish to estimate.
#'
#' @param M Numeric vector of thresholds.
#' @param scale Whether we are interested in the raw values of M or are treating
#'        them as fold-changes from baseline ('proportional') or absolute
#'        differences ('difference')
#' @param trans A function for transforming the data, the same as the function
#'        used to transform the data before fitting the original linear model.
#' @param d A list of data.frames
#' @param baseline The name of the baseline column in d
margarita.rp.matrix <- function(M, scale, trans, d, baseline){
    # Get baseline data - should be identical for each element of d, so use d[[1]]
    baseline <- d[[1]][, baseline]
    nr <- nrow(d[[1]])

    if (scale=="r"){ # raw
        m <- matrix(rep(trans(M), nr), ncol=length(M), byrow=TRUE)
    }
    else if (scale=="p") { # M is a multiple of baseline
        m <- matrix(rep(0, nr * length(M)), ncol=length(M))
        # d[[i]][, baseline] is the same for all i
        for (i in 1:length(M)){
            m[, i] <- trans(M[i] * baseline)
        }
    }
    else if (scale=="d"){ # difference
        m <- matrix(rep(0, nr * length(M)), ncol=length(M))
        for (i in 1:length(M)){
            m[, i] <- trans(M[i] + baseline)
        }
    }
    m
}

#' Annotate a threshold selection ggplot
#'
#' Annotate a threshold selection ggplot with the number
#' of exceedances of various thresholds.
#'
#' @param p An object produced by ggplot
#' @param x Horizontal axis data containing the full range.
#' @param y Verticle axis data containing the full range.
#' @param data The actual data being considered for GPD modelling.
#' @param u Thresholds above which we are interested.
#' @param textsize The size of the text in the annotations.
addExcesses <- function(p, x, y, data, u, textsize){
    x1 <- axisTicks(range(x), log=FALSE)
    yr <- range(y)
    delta <- abs(diff(yr)) * .1
    y1 <- rep(yr[2] + delta, length(x1))
    txt <- sapply(x1, function(u) sum(data > u))
    tx=data.frame(ex="Excesses:", x=min(x), y=y1[1] + delta)
    df <- data.frame(x=x1, y=y1, txt=txt)
    p <- p + geom_text(data=tx, aes(x,y,label=ex), size=textsize, hjust=0)
    p + geom_text(data=df, aes(x, y, label=txt), size=textsize)
}

#' Check an allowed value of the scale argument has been given and reduce it to
#' its first letter
#'
#' @param s A character string which should be 'raw', 'proportional' or 'difference'
#'         or an abbreviation of one of those. Only the first letter is returned.
margaritaScale <- function(s){
    if (length(s) > 1 | !is.character(s)){
        stop("scale should be a character string of length 1")
    }
    s <- substring(casefold(s), 1, 1)
    if (!is.element(s, c("r", "p", "d"))){
        stop("scale can be raw, proportional or difference (only the first letter is checked)")
    }
    else{
        s
    }
}
