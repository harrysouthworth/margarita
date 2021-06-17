#' @importFrom stats simulate
#' @importFrom utils head
NULL

#' Create an object of class 'margarita'
#' @param rlm An object of class 'rlm' returned by \code{lmr}.
#' @param evmSim An object of class 'evmSim'.
#' @param newdata A \code{data.frame} containing the new data from which
#'        predictions are to be made. Defaults to \code{newdata=NULL} which will
#'        only work of there are no covariates in the linear model or the extreme
#'        value model (with the exception of the baseline term in the linear model)
#' @param trans A function matching that used to transform the response
#'        in the robust regression. Defaults to \code{trans=log}. Use
#'        \code{trans=I} if no transformation was made.
#' @param invtrans A function, the inverse of \code{trans}.
#' @param alpha Used to specify coverage of interval estimates. Defaults to
#'   \code{alpha = c(1, .5)}.
#' @param baseline Character string giving the name of the baseline variable.
#' @param arm Character string giving the name of the treatment group variable.
#'        Defaults to \code{arm = "arm"}. If there is no treatment factor,
#'        presumably there is a continuous exposure variable and the sample
#'        size is the number of rows in the data.
#' @param minima Are the extremes minima rather than maxima (i.e. the response
#'        was multiplied by -1 prior to all modelling). Defaults to
#'        \code{minima=FALSE}. Note that if you specify \code{minima=TRUE} the
#'        function will fail and you will be told to multiply by -1 at the very
#'        start of the analysis, and to then multiply M by -1 in the call to
#'        this function.
#' @details Returns a list with class 'margarita', to be used with
#'       \code{simulate.margarita}.
#' @keywords models
#' @export margarita
margarita <- function(rlm, evmSim, newdata=NULL,
                      trans=log, invtrans=exp,
                      alpha = c(.1, .5),
                      baseline="alt.b", arm="arm", minima=FALSE){
  if (minima){
    stop("If working with minima, multiply the response and baseline by -1 at the very start. Then call this function specifying -M, not M")
  }

  if (!inherits(rlm, c("rlm", "lmrob"))){ # class is likely c("lmr", "rlm")
    stop("object should have class 'rlm'")
  }
  else if (is.null(rlm$cov)){
    stop("object should be created by lmr, not rlm (I need the covariance matrix)")
  }
  if (class(evmSim) != "evmSim"){
    stop("object should have class 'evmSim'")
  }

  # Check each row of newdata is unique
  if (!is.null(newdata)){
    ur <- nrow(unique(newdata))
    if (ur != nrow(newdata))
      stop("newdata should have unique rows")
    if (ncol(newdata) > 1)
      stop("newdata should have only one column")
  }

  # Ensure factor levels are in the same order in newdata and rlm
  if (length(rlm$xlevels) > 0 && !all(levels(newdata[, 1]) == rlm$xlevels[[1]]))
      stop("Levels of the factor in newdata don't match those in the robust linear model (it might just be the ordering)")

  # Check trans and invtrans mirror each other
  b <- rlm$data[, baseline]
  bt <- try(trans(b), silent=TRUE)
  if ("try-error" %in% class(bt)){
    stop("trans(baseline) results in errors")
  }
  ib <- try(invtrans(bt))
  if ("try-error" %in% class(ib)){
    stop("itrans(trans(baseline)) results in errors")
  }

  cbt <- cor(b, invtrans(bt))

  if (is.na(cbt) | !all.equal(cbt, 1)){
    stop("trans and invtrans don't undo each other")
  }


  # Construct string for transformed baseline
  if (missing(trans)){
    ctrans <- "log"
  } else {
    ctrans <- as.character(as.list(match.call())$trans)
  }
  rawBaseline <- baseline
  baseline <- if (ctrans != "I") paste0(ctrans, "(", rawBaseline, ")")
  else baseline

  if (! baseline %in% colnames(rlm$x))
    stop(paste0("Baseline variable '", baseline, "' is not in the linear model"))

  if (is.null(newdata))
    if (length(coef(evmSim)) == 2) newdata <- data.frame(1)
  else stop("There are covariates in the model: you need to provide newdata")

  # Get sample sizes so that we can later adjust the distribution of the expected
  # proportions exceeding thresholds
  if (!(arm %in% names(rlm$data))){
    stop("arm not in rlm$data")
  } else if (inherits(rlm$data[, arm], c("factor", "character"))){
    n <- table(rlm$data$arm)
  } else {
    n <- nrow(rlm$data)
  }

  # Construct output object and return
  res <- list(rlm, evmSim, newdata=newdata,
              baseline=baseline, rawBaseline=rawBaseline,
              arm=arm, n=n,
              trans=trans, invtrans=invtrans,
              alpha = alpha,
              minima=minima)
  oldClass(res) <- "margarita"
  res
}

#' @method print margarita
#' @export
print.margarita <- function(x, ...){
  print(x[[1]])
  cat("\n")
  print(x[[2]])
}
