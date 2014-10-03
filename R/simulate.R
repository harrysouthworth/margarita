#' Simulate baseline data used in construction of robust linear model and
#' predictions from extreme value model
#'
#' Simulate baseline data by resampling, inferring the number of samples from the input.
#'
#' @param lmod An object of class 'rlm' created by \code{lmr}.
#' @param gmod An object of class 'evmSim'.
#' @param baseline Character string identifying the baseline column of the
#'        design matrix used to fit the robust linear model.
#' @details The function decides how many baselines to simulate based on the
#'        number of simulated values in the \code{param} component of \code{gmod}.
#' @keywords datagen
#' @export simBase
simBase <- function(lmod, gmod, baseline="log(alt.b)"){
    if (class(gmod) != "evmSim"){
        stop(paste("I need an object of class evmSim, not", class(gmod)))
    }
    nsim <- nrow(gmod$param)
    x <- lmod$x
    x <- x[, colnames(x) == baseline]
    res <- sample(x, size=nsim, replace=TRUE)
    invisible(res)
}

#' Simulate expected responses from a robust linear model, accounting for
#' uncertainty in the parameter estimates.
#' @param lmod An object returned by \code{lmr}.
#' @param gmod An object of class 'evmSim'.
#' @param newdata A data.frame containing values of the covariates used
#'        to fit the robust model, excluding the baseline variable.
#' @param baseline Character string identifying the baseline column of the
#'        design matrix used to fit the robust linear model. Defaults to
#'        'log(b.alt)'.
#' @param rawBaseline Character string giving the name of the untransformed
#'        baseline values. If no transformation is performed, this will be
#'        the same value as for \code{baseline}. Defaults to
#'        \code{rawBaseline = "alt.b"}.
#' @param invtrans A function for back-transforming the baseline data if required.
#'        By default, the function assumes the user worked with log-transformed
#'        baseline data, so assumes \code{invtrans=exp}.
#' @keywords datagen
#' @importFrom mvtnorm rmvnorm
#' @export simLinear
simLinear <- function(lmod, gmod, newdata=NULL,
                      baseline="log(alt.b)", rawBaseline="alt.b", invtrans=exp){
    # Check each row of newdata is unique
    ur <- nrow(unique(newdata))
    if (ur != nrow(newdata)){
        stop("newdata should have unique rows")
    }
    # Get simulated baselines
    b <- simBase(lmod, gmod, baseline=baseline)
    nsim <- length(b)
    nrep <- nrow(newdata)

    # Repeat each row of newdata nsim times, add baseline
    newdata <- newdata[rep(1:nrep, each=nsim),, drop=FALSE]

    newdata[, rawBaseline] <- invtrans(b) # b gets repeated nrow times
    # Get design matrix and simulate parameters from linear model
    fo <- lmod$call$formula
    fo[[2]] <- NULL

    M <- model.matrix(fo, newdata, xlev=lmod$xlevels)
    M <- M[, names(coef(lmod)), drop=FALSE]

    # GPD parameters are reused for each observation in newdata,
    # and so are baselines, so do the same with the simulated
    # coefficients here.
    lco <- rmvnorm(nsim, coef(lmod), lmod$cov) -> co1
    if (nrep > 1)
      for (i in 1:(nrep - 1))
          lco <- rbind(lco, co1)

    res <- rowSums(M * lco)
    newdata$fitted <- res

    invisible(newdata)
}

#' Simulate return levels, probabilities of threshold exceedences or simple datasets from a robust linear model and an extreme value model
#'
#' Simulate return levels and probabilities of threshold exceedences from a
#' robust linear model and an extreme value model. The procedure is specific to
#' the situation of extreme value modelling of residuals when predictions on
#' the scale of the original data are required.
#'
#' @aliases simulate.margarita.rl simulate.margarita.prob simulate.margarita.simple
#' @aliases as.data.frame.margarita.sim.rl print.margarita.sim.rl head.margarita.sim.rl
#' @aliases as.data.frame.margarita.sim.prob print.margarita.sim.prob head.margarita.sim.prob
#' @aliases summary.margarita.sim.rl summary.margarita.sim.prob
#' @aliases ggplot.summary.margarita.sim.rl ggplot.summary.margarita.sim.prob
#' @aliases as.data.frame.summary.margarita.sim.rl print.summary.margarita.sim.rl
#' @aliases as.data.frame.summary.margarita.sim.prob print.summary.margarita.sim.prob
#' @aliases "/.summary.margarita.sim.rl" "*.summary.margarita.sim.rl" "-.summary.margarita.sim.rl" "+.summary.margarita.sim.rl"
#' @param object An object of class 'margarita'
#' @param nsim Unused argument.
#' @param seed Unused argument.
#' @param type What type of prediction is required: either 'rl', 'prob' or 'simple'.
#'        Defaults to \code{type = "rl"}. If \code{type='simple'} is used, the result
#'        is just a \code{data.frame} containing simulated baselines, fitted values
#'        (from the robust linear model) and values resulting from adding residuals
#'        to the fitted values. The residuals are sampled from GP distributions with
#'        parameters taken from the Markov chains, or are resampled from observed
#'        residuals, in proportion to the number of observations above and below
#'        the GP fitting threshold.
#' @param M The return level to be predicted. Defaults to \code{M=1000}. If
#'        \code{type="prob"}, M should be a vector containing the thresholds
#'        whose probabilities of exceedance the user is interested in, on the
#'        scale of the raw data.
#' @param scale The type of prediction being made. Valid values are 'raw',
#'        'proportional' and 'difference'. Due to how the calculations are
#'        performed, \code{scale} needs to be specified in the call to
#'        \code{simulate} if threshold exceedence probabilities are being
#'        predicted, or the call to \code{summary.simulate} if return levels
#'        are being predicted. If \code{scale='raw'} it is assumed that the
#'        values of \code{M} are on the same scale as the response variable
#'        in the GPD model. If \code{scale='proportional'}, the values of
#'        \code{M} are taken to be fold-changes (e.g. 2, 5, 10 times baseline
#'        for return levels, or P(yMax > 2, 5, 10 times yBase)). If
#'        \code{scale='difference'}, absolute changes from baseline are
#'        assumed.
#' @param Mlabels Labels to be used in the output. Defaults to \code{Mlabels=NULL}
#'        in which case the function guesses at meaningful labels.
#' @param ... Other agruments passed to \code{simulate}. Currently unused.
#' @keywords datagen
#' @method simulate margarita
#' @export
simulate.margarita <- function(object, nsim=1, seed=NULL,
                               type="rl", M=NULL, scale="raw", Mlabels=NULL, ...){
    if (missing(M) & type != "simple"){
        stop("You must provide a value for M")
    }

    res <- switch(type,
                  "rl" =, "return" =, "return level" =
                      simulate.margarita.rl(object, nsim=nsim, seed=seed, M=M, ...),
                  "prob" =, "probability"=, "excess probability" = 
                      simulate.margarita.prob(object, nsim=nsim, seed=seed, M=M, scale=scale, Mlabels=Mlabels, ...),
                  "simple"= simulate.margarita.simple(object, nsim=nsim, seed=seed, ...)
                  )
    attr(res, "baseline") <- object$rawBaseline
    attr(res, "trans") <- object$trans
    attr(res, "invtrans") <- object$invtrans
    
    if (type %in% c("rl", "return", "return level"))
      oldClass(res) <- "margarita.sim.rl"
    else if (type %in% c("prob", "probability", "excess probability"))
      oldClass(res) <- "margarita.sim.prob"
    # otherwise it's just a data.frame
    
    res
}

simulate.margarita.rl <- function(object, nsim=1, seed=NULL, M=NULL, ...){
    object <- unclass(object)

    res <- simLinear(object[[1]], object[[2]], newdata=object$newdata,
                     baseline=object$baseline, rawBaseline=object$rawBaseline,
                     invtrans=object$invtrans)
    fullres <- list()

    for (m in M){
      nmM <- paste("RL:", m)
      fullres[[nmM]] <- res
      p <- predict(object[[2]], newdata=object$newdata, all=TRUE, M=m, unique.=FALSE)
      p <- c(unclass(p)[[1]])
      
      fullres[[nmM]]$RLraw <- p
      fullres[[nmM]]$RLfull <- fullres[[nmM]]$RLraw + res$fitted
      if (object$minima) { fullres[[nmM]]$RLfull <- -fullres[[nmM]]$RLfull }
    }
    fullres
}

simulate.margarita.prob <- function(object, nsim=1, seed=NULL, M=NULL, scale="raw", 
                                    Mlabels=NULL, ...){
    object <- unclass(object)
    scale <- margaritaScale(scale)

    res <- simLinear(object[[1]], object[[2]], newdata=object$newdata,
                     baseline=object$baseline, rawBaseline=object$rawBaseline,
                     invtrans=object$invtrans)

    # Need res to be a list of length 4
    nr <- nrow(res)
    nn <- nrow(object$newdata)
    g <- rep(1:nn, each=nr/nn)

    res <- split(res, g)

    par <- predict(object[[2]], newdata=object$newdata, type="lp", unique.=FALSE, all=TRUE)[[1]]
    # par is a list with an entry for each row of newdata.
    # Each entry is a matrix containing the parameters in the first two columns.

    # Get thresholds, put them in u
    u <- lapply(res, function(x, u){ x$fitted + u }, u=object[[2]]$map$threshold)

    # Get matrix of M. Do it this way so that we can treat M as a multiple of baseline
    m <- margarita.rp.matrix(M, scale=scale, trans=object$trans,
                             d=res, baseline=object$rawBaseline)

    out <- lapply(1:nn, margarita.getProbs, u = u, par=par, m=m,
                                 r=object[[2]]$map$data$y, p = object[[2]]$map$rate)
    out <- lapply(out, do.call, what="cbind")
    # out is a list. Each element is a matrix with one column for each value of M

    onms <- apply(object$newdata, 2, as.character)
    # XXX NEXT LINE HAS HAD PROBLEMS
    # XXX Check the Git repo for previously attempted solutions if this fails again
    if (!is.matrix(onms)) onms <- as.data.frame(t(onms))
    onms <- apply(onms, 1, paste, collapse=" ")
    names(out) <- onms
    
    if (is.null(Mlabels)) Mlabels <- paste("RL =", format(M))

    out <- lapply(out, function(x, m){ colnames(x) <- m; x }, m=Mlabels)

    out
}

# Simulate pairs of baselines and on-treatment values
simulate.margarita.simple <- function(object, nsim=1, seed=NULL, ...){
  if (class(object) != "margarita") stop("object should be of class margarita")
  s <- simLinear(object[[1]], object[[2]], newdata=object$newdata,
                 baseline=object$baseline, rawBaseline=object$rawBaseline,
                 invtrans=object$invtrans)
  param <- object[[2]]$param # Posterior simulated parameter estimates
  th <- object[[2]]$map$threshold # Threshold used in GPD fit
  ru <- rgpd(nrow(param), exp(param[, 1]), param[, 2], u=th)
  rl <- resid(object[[1]])
  rl <- sample(rl[rl <= th], size=nrow(param), replace=TRUE)
  r <- sample(c(ru, rl))
  res <- s
    #data.frame(baseline=s[, object$rawBaseline], max=s$fitted + sample(r, size=nrow(s)))
  pe <- object[[2]]$map$rate # sample from ru with probability P(x > th)
  res$value <- res$fitted + sample(r, size=nrow(s), prob=rep(c(pe, 1-pe), each=nrow(param)))
  res$value <- object$invtrans(res$value)
  rownames(res) <- NULL
  invisible(res)
}


#' @method as.data.frame margarita.sim.rl
#' @export
as.data.frame.margarita.sim.rl <-
function(x, row.names=NULL, optional=FALSE, ...){
  x <- lapply(x, function(x) as.data.frame(unclass(x)))

  M <- rep(names(x), each=nrow(x[[1]]))
  x <- do.call("rbind", x)
  cat("sick\n")
  x <- as.data.frame(x)
  cat("scratch\n")
  x$M <- M
  x
}

#' @method print margarita.sim.rl
#' @export
print.margarita.sim.rl <-
function(x, ...){
    x <- as.data.frame(x)
    print(x, ...)
}

#' @method head margarita.sim.rl
#' @export
head.margarita.sim.rl <-
function(x, ...){
    x <- as.data.frame(x)
    head(x, ...)
}

#' @method summary margarita.sim.rl
#' @export
summary.margarita.sim.rl <- function(object, alpha=c(.1, .5), scale="raw", ...){
    baseline <- attr(object, "baseline")
    invtrans <- attr(object, "invtrans")
    scale <- margaritaScale(scale)
    
    # Get quantiles for credible interval
    qu <- getCIquantiles(alpha)
       
    sims <- c(baseline, "fitted", "RLraw", "RLfull")
    groups <- names(object[[1]])[!is.element(names(object[[1]]), sims)]

    groups <- object[[1]][, groups, drop=FALSE]
    groups <- apply(groups, 2, as.character)
    groups <- apply(groups, 1, paste, collapse=" ")

    getSummaries <- function(o){
      o <- as.data.frame(o)
      # Get response: the raw values, proportional or absolute change from baseline
      x <- switch(scale,
                  "r" = invtrans(o$RLfull),
                  "p" = invtrans(o$RLfull)/o[, baseline],
                  "d" = invtrans(o$RLfull) - o[, baseline])
      
      res <- tapply(x, groups, quantile, prob=qu)
      # Groups have been reordered alphabetically. Revert to original order
      res <- res[unique(groups)]
      res <- t(do.call("data.frame", res))
      colnames(res) <- paste("Q", qu*100, "%", sep = "")
      rownames(res) <- unique(groups)
      res
    }
    res <- lapply(object, getSummaries)

    oldClass(res) <- "summary.margarita.sim.rl"
    res
}

#' @method as.data.frame summary.margarita.sim.rl
#' @export
as.data.frame.summary.margarita.sim.rl <- function(x, row.names=NULL, optional=FALSE, ...){
  M <- rep(names(x), each=nrow(x[[1]]))
  groups <- rep(rownames(x[[1]]), length.out=length(x) * nrow(x[[1]]))

  x <- do.call("rbind", x)
  x <- as.data.frame(x)

  x$groups <- groups
  x$M <- M

  rownames(x) <- 1:nrow(x)
  x
}


#' @method print summary.margarita.sim.rl
#' @export
print.summary.margarita.sim.rl <- function(x, ...){ print(as.data.frame(x), ...) }

#' @method print margarita.sim.prob
#' @export
print.margarita.sim.prob <- function(x, ...){
    print(unclass(x), ...)
}

#' @method head margarita.sim.prob
#' @export
head.margarita.sim.prob <- function(x, ...){
    x <- unclass(x)
    lapply(x, head)
}

#' @method summary margarita.sim.prob
#' @export
summary.margarita.sim.prob <- function(object, alpha=c(.1, .5), ...){
    object <- unclass(object)

    qu <- getCIquantiles(alpha)    
    
    summat <- function(x, qu){
        t(apply(x, 2, quantile, prob=qu))
    }
    
    res <- lapply(object, summat, qu=qu)

    oldClass(res) <- "summary.margarita.sim.prob"
    res
}

#' @method print summary.margarita.sim.prob
#' @export
print.summary.margarita.sim.prob <- function(x, ...){
    for (i in 1:length(x)){
        cat(names(x)[i], "\n")
        print(x[[i]])
        cat("\n")
    }
    invisible()
}
