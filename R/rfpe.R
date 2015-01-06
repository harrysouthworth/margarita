#' Robust final prediction error for a linear model
#' @param x A robust linear model, fit by \code{lmr}.
#' @param scale The value of the scale to use. Defaults to \code{scale=NULL} and is take from \code{x}.
#' @return A number representing the robust final prediction error.
#' @details The definition of robust final prediction error is given in Section 5.12
#'     of Maronna et al. The function is generic, but so far only methods for objects
#'     fit by \code{lmr} and \code{lmRob} have been implemented (the latter purely for
#'     testing). It would be straightforward to implement
#'     for objects returned by \code{rlm} or \code{lmrob}, for example. Only bisquare
#'     weight functions are considered.
#' @references Maronna, R. A, Martin, R. D and Yohai, V. J. (2006) Robust Statistics: Theory and Methods, Wiley
#' @aliases RFPE.lmr RFPE.lmRob
#' @export RFPE
RFPE <- function(x, scale=NULL) NextMethod("RFPE", x)

#' @method RFPE lmr
#' @export RFPE.lmr
RFPE.lmr <- function (x, scale=NULL){
  if (is.null(scale)) scale <- x$s
  r <- resid(x) / scale
  if (any(is.na(r))) stop("missing values are not allowed")

  c <- x$c

  a <- mean(bisquare(r, c=c, d=0))
  A <- mean(bisquare(r, c=c, d=1)^2)
  B <- mean(bisquare(r, c=c, d=2))

  q <- length(x$coef) # q, page 151
  n <- length(r)

  if (B <= 0) return(NA)
  a + (q*A) / (n*B)
}

#' @method RFPE lmRob
#' @export RFPE.lmRob
RFPE.lmRob <- function(x, scale=NULL){
  if (is.null(scale)) scale <- x$scale
  x$s <- scale
  x$c <- x$yc
  RFPE.lmr(x, scale)
}

#' Bisquare weight function and its derivatives
#' @param x Scaled residuals from a linear model.
#' @param c The tuning constant to the bisquare function.
#' @param d The order of derivative of the bisquare function required. Defaulgs
#'        to \code{d=0}, the other allowed values being 1 and 2.
bisquare <- function(x, c=3.443689, d=0){
  if (d == 0){
    1 - (1 - pmin(1, abs(x/c))^2)^3
  } else if (d == 1){
    ifelse(abs(x) < c, (6*x / (c^6)) * (c^2 - x^2)^2, 0)
  } else if (d == 2){
    ifelse(abs(x) < c, 6/c^6 * (c^4 - 6*(c*x)^2 + 5*x^4), 0)
  } else {
    stop("d can take values 0, 1 or 2")
  }
}


#' Robust version of AIC for an MM-estimated linear model
#' 
#' @param object An object of class 'lmr'.
#' @param scale An optional scale parameter. If not provided, it is taken from \code{object}.
#' @param k The penalty per parameter to be used.
#' @param ... Not used.
#' @details The robust AIC correponds to equation (16) of Tharmaratnam and Claeskens.
#'     The 'penalty' term is not a simple function of the number of parameters in the
#'     model, so increaseing k does not necessarily result in simpler models being
#'     chosen.
#' @references Tharmaratnam, K. and Claeskens, G. A comparison of robust versions
#'      of the AIC based on M, S and MM-estimators, Technical Report KBI 1014,
#'      Katholieke Universiteit Leuven, 2010.
#'      https://lirias.kuleuven.be/bitstream/123456789/274771/1/KBI_1014.pdf
#' @method AIC lmr
#' @export AIC.lmr
AIC.lmr <- function(object, scale, ..., k=2){
  if (missing(scale)) scale <- object$s
  r <- resid(object) / scale
  n <- length(r)
  X <- object$x
 
  a <- bisquare(r, c=object$c, d=2)
  b <- bisquare(r, c=object$c, d=1)^2

  iJ <- solve(t(X) %*% diag(a) %*% X * (1/(scale^2)) / n)
  K <- (t(X) %*% diag(b) %*% X * (1/(scale^2))) / n

  2 * n * log(scale) + k * sum(diag(iJ %*% K))
}
