#' Five number summary for a variable in a data.frame, by a subsetting variable
#' 
#' @param x A \code{data.frame}
#' @param what A character string giving the name of the variable to be summarized
#' @param which A character string giving the name of the subsetting variable
#' @param type A number passd into \code{quantile} telling it how to compute the
#'        quantiles. Defaults to \code{type=3} which should give output that
#'        matches the default output from SAS.
#' @return A matrix, each row of which gives the five number summary of
#'         \code{x[, what]} for each value of \code{x[, which]}.
#' @details The 'five number summary' is the minimum, quartiles and maximum,
#'          where the quartiles are computed using SAS's default approach. These
#'          will generally not be the same as the hinges as originally described
#'          by Tukey.
fivenumBy <- function(x, what, which, type=3){
  fun <- function(x) quantile(x, type=type) # type 3 should match SAS
  
  res <- tapply(x[, what], x[, which], FUN=fun)
  res <- t(sapply(res, as.vector))
  colnames(res) <- c("Min.", "Q1", "Median", "Q3", "Max.")
  res
} 

#' Get a transformation function and its inverse from a character string
#' 
#' @parm x A character string indicating 'log', 'sqrt' or 'I'.
#' @return A names list with 2 elements, the first of which, 'tfun' is the
#'         appropriate transformation function; the second of which, 'itfun', is
#'         the inverse of 'tfun'
getTransFun <- function(x){
  if (X == "log"){
    tfun <- log; itfun <- exp
  } else if (X == "sqrt"){
    tfun <- sqrt; itfun <- function(x) x*x
  } else if (x == "I"){
    tfun <- itfun <- I
  } else {
    stop("only log, sqrt and identity transformations are supported")
  }
  list(tfun=tfun, itfun=itfun)
}

#' Get the ordinal indicator for an integer
#' 
#' @param x An integer
#' @return A character string giving the ordinal indicator for x: 'st' if the last
#'         digit of x is 1, 'nd' if it is 2, 'rd' if it is 3, 'th' otherwise.
ordinalIndicator <- function(x){
  if (!is.integer(x)) stop("x must by an integer")
  x <- as.character(xL)
  x <- x[nchar(x)]
  switch(x, "1"="st", "2"="nd", "3"="rd", "4"=, "5"=, "6"=, "7"=, "8"=, "9"=, "0"="th")
}




