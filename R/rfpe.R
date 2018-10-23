#' Robust final prediction error for a linear model
#' @param x A robust linear model, fit by \code{margarita::lmr} or
#'   \code{robust::lmRob} or \code{robustbase::lmrob}, or a list containing
#'   multiple such models.
#' @param scale The value of the scale to use. If \code{x} is a list and \code{scale}
#'   is not provided by the user, the function takes the scale estimate to be the
#'   one associated with the largest model. Otherwise, no default is available. See details.
#' @return A number representing the robust final prediction error or, if \code{x}
#'   is a list, a \code{data.frame} giving the final prediction error and other
#'   details.
#' @details The definition of robust final prediction error is given in Section 5.12
#'     of Maronna et al. The function is generic, but so far only methods for objects
#'     fit by \code{lmr}, \code{lmrob}, \code{lmRob} and for a list of such
#'     models have been implemented. Only bisquare
#'     weight functions are considered. \strong{Note that} the purpose of RFPE
#'     is to \emph{compare} models -- for a single model, it's not very meaningful.
#'     To allow a valid comparison, the calculation must use the same scale
#'     estimate, usually the one taken from the largest model. For this reason,
#'     no default scale argument is provided, forcing the user to specify a value.
#' @references Maronna, R. A, Martin, R. D and Yohai, V. J. (2006) Robust Statistics: Theory and Methods, Wiley
#' @aliases RFPE.lmr RFPE.lmRob
#' @export RFPE
RFPE <- function(x, scale){
  UseMethod("RFPE", x)
}

#' @export
RFPE.lmr <- function (x, scale){
  r <- resid(x) / scale
  if (any(is.na(r))) stop("missing values are not allowed")

  if (class(x)[1] == "lmr"){
    c <- x$c
  } else if (class(x) == "lmrob"){
    c <- x$control$tuning.psi
  } else if (class(x) == "lmRob"){
    c <- x$yc
  }

  a <- mean(bisquare(r, c=c, d=0))
  A <- mean(bisquare(r, c=c, d=1)^2)
  B <- mean(bisquare(r, c=c, d=2))

  q <- length(x$coef) # q, page 151
  n <- length(r)

  if (B <= 0) return(NA)
  a + (q*A) / (n*B)
}

#' @export
RFPE.lmrob <- function(x, scale){
  if (x$control$psi != "bisquare"){
    stop("Only bisquare weight functions are implemented")
  }

  RFPE.lmr(x, scale)
}

#' @export
RFPE.lmRob <- function(x, scale){
  if (unique(x$robust.control$weight) != "bisquare"){
    stop("Only bisquare weight functions are implemented")
  }

  x$s <- scale
  x$c <- x$yc
  RFPE.lmr(x, scale)
}


#' @export
RFPE.list <- function(x, scale){
  # Check response vectors are all the same
  y <- lapply(x, function(X){
    X$residuals + X$fitted.values
  })
  tt <- sapply(y[2:length(y)], FUN = all.equal, y[[1]])

  if (class(tt) == "character"){
    stop("The models in x don't all have the same response vector")
  }

  if (is.null(names(x))){
    names(x) <- 1:length(x)
  }

  # Reverse order the models by rank
  rank <- sapply(x, function(X) X$rank)
  x <- x[rev(order(rank))]

  if (missing(scale)){
    if (class(x[[1]])[1] == "lmr"){
      scale <- x[[1]]$s
    } else {
      scale <- x[[1]]$scale
    }
  }

  # Check if maximal model is bigger than all others
  srank <- rev(sort(rank))
  if (rank[1] == rank[2]){
    message(paste("More than one model with maximum rank. Using scale =", scale))
  }

  rfpe <- sapply(x, RFPE.lmr, scale = scale)

  res <- data.frame(Model = names(x), Rank=rank[rev(order(rank))],
                    formula = as.character(sapply(x, function(X) X$call$formula)),
                    RFPE = rfpe, `Rel. RFPE` = rfpe / min(rfpe),
                    stringsAsFactors=FALSE, check.names=FALSE)
  res <- res[order(rfpe), ]

  attr(res, "scale") <- scale

  rownames(res) <- 1:nrow(res)

  res
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
