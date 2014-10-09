readAutoInputs <- function(file){
  input <- readLines(file)
  ii <- substr(wh, 1, 1)
  input <- input[ii != "#"] # drop comment lines

  input <- strsplit(input, "#")

  # i is a list with 1 element for each line of the input file
  # i[[j]][1] gives the object name and value, anything else is junk

  input <- unlist(lapply(input, function(x) x[1]))

  output <- strsplit(input, ": ")
  nms <- unlist(lapply(output, function(x) x[1]))
  output <- unlist(lapply(output, function(x) x[2]))
  names(output) <- nms
  output
}


#' Read a data file, figuring out if it is a SAS file or something else
#' 
#' @param file The name, including the full path, of the data file
#' @param type The type, as indicated by the file extension. At present, must
#'             be either 'sas7bdat', 'csv' or NULL (the default). If NULL, the
#'             function uses the file extension to guess.
#
#' @details If the file appears to be a SAS file, the function uses read.sas7bdat
#'          from Matt Shotwell's sas7bdat package. Currently, compressed SAS
#'          data files are not supported (but see Shotwell's sas7bdat.parso
#'          package).
#' @importFrom tools file_ext
#' 
readData <- function(file, type=NULL){
  if (is.null(type)) type <- file_ext(file)
  if (type == "sas7bdat"){
    res <- read.sas7bdat(file)
  } else (type = ".csv") {
    res <- read.csv(file)
  } else {
    stop("File type (extension) must be sas7bdat or csv")
  }
  invisible(res)
}

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
#' @export
fivenumBy <- function(x, what, which, type=3){
  fun <- function(x) quantile(x, type=type) # type 3 should match SAS
  
  res <- tapply(x[, what], x[, which], FUN=fun)
  res <- t(sapply(res, as.vector))
  colnames(res) <- c("Min.", "Q1", "Median", "Q3", "Max.")
  res
} 

#' Get a transformation function and its inverse from a character string
#' 
#' @param x A character string indicating 'log', 'sqrt' or 'I'.
#' @return A names list with 2 elements, the first of which, 'tfun' is the
#'         appropriate transformation function; the second of which, 'itfun', is
#'         the inverse of 'tfun'
#' @export
getTransFun <- function(x){
  if (x == "log"){
    tfun <- log; itfun <- exp
  } else if (x == "sqrt"){
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
#' @export
ordinalIndicator <- function(x){
  if (!is.integer(x)) stop("x must by an integer")
  x <- as.character(x)
  x <- substring(x, nchar(x))
  switch(x, "1"="st", "2"="nd", "3"="rd", "4"=, "5"=, "6"=, "7"=, "8"=, "9"=, "0"="th")
}




