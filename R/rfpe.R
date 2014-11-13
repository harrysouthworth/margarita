#' Robust final prediction error for a linear model
#' @param x A robust linear model, fit by \code{lmr}.
#' @return A number representing the robust final prediction error.
#' @details The definition of robust final prediction error is given in Section 5.12
#'     of Maronna et al. The function is generic, but so far only methods for objects
#'     fit by \code{lmr} and \code{lmRob} have been implemented (the latter purely for
#'     testing). It would be straightforward to implement
#'     for objects returned by \code{rlm} or \code{lmrob}, for example. Only bisquare
#'     weight functions are considered.
#' @references Maronna, R. A, Martin, R. D and Yohai, V. J. (2006) Robust Statistics: Theory and Methods, Wiley
#' @export RFPE
RFPE <- function(x) NextMethod("RFPE", x)

#' @method RFPE lmr
#' @export RFPE.lmr
RFPE.lmr <- function (x){
  r <- resid(x) / x$s
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
RFPE.lmRob <- function(x){
  x$s <- x$scale
  x$c <- x$yc
  RFPE.lmr(x)
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

drop1.lmr <- function(object, scope, ...){
  
  
}


step.lmr <- function(object, scope, direction="backwards", target=RFPE){
  
  
  
}
