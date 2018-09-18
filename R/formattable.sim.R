#' Get a vector of grey colors
#' @param n The number of shades of grey to return.
#' @details The darkest value returned will be 'grey51' and the lightest will
#'   be 'grey100' which is white. (The darkest, black, is 'grey0')
shadesOfGrey <- shadesOfGray <-
  colorRampPalette(c("grey51", "grey100"))

#' Create an HTML table for outputing via Rmd.
#' @aliases formattable.summary.margarita.sim.prob
#' @param data The output of \code{summary.margarita.sim.rl}.
#' @param digits The number of digits to round numbers to. Defaults to
#'   \code{digits=2}.
#' @param backgroundColor Rows are colored according to treatment group. Defaults
#'   to using \code{shadesOfGrey} to select the correct number of shades of grey.
#' @param debug If \code{debug = TRUE}, the list of formatters is returned for
#'   inspection.
#' @note The returned object (if \code{debug = FALSE}) is a formattable object
#'   and, if desired, can be turned into a \code{DT::datatable} via
#'   \code{formattable::as.datatable}.
#' @method formattable summary.margarita.sim.rl
#' @export
formattable.summary.margarita.sim.rl <-
  function(data, digits=2, backgroundColor=NULL, debug=FALSE){
    data <- as.data.frame(data)

    area <- formattable::area
    formatter <- formattable::formatter
    formattable <- formattable::formattable

    ng <- length(unique(data$groups))
    nrl <- length(unique(data$M))

    fstring <- paste0("%.", digits, "f")

    if (is.null(backgroundColor)){
      backgroundColor <- shadesOfGrey(ng)
    } else if (length(backgroundColor) != ng){
      stop("backgroundColor should have length to match the number of groups")
    }

    data[, 1:5] <- apply(data[, 1:5], 2, function(X) sprintf(fstring, X))

    data <- mutate(data, M = factor(data$M, levels=unique(data$M))) %>%
      select(Group = groups, `#Subjects` = M, starts_with("Q"))


    areaList <- formatterList <- vector(mode = "list", length = nrl * 2)

    rowSeqs <- matrix(1:(ng * nrl), ncol = nrl, byrow=FALSE)
    rowSeqs <- apply(rowSeqs, 1, function(X) paste(X, collapse = ", "))

    for (j in 1:2){ # first is all columns, second is column 5 (Q50%)
      for (i in 1:ng){
        index <- i + (j - 1) * ng
        if (j == 1){
          a <- paste0("area(row=c(", rowSeqs[i], "))")
          ff <- paste0("formatter('span', style = x ~ style('background-color'='", backgroundColor[i], "'))")
          areaList[[index]] <- as.formula(paste(a, "~", ff))
        } else {
          a <- paste0("area(row=c(", rowSeqs[i], "), col=5)")
          ff <- paste0("formatter('span', style = x ~ style('background-color'='", backgroundColor[i], "', 'font-weight'='bold'))")
          areaList[[index]] <- as.formula(paste(a, "~", ff))      }
      }
    }

    if (debug){
      areaList
    } else {
      formattable(data, areaList)
    }
  }

#' @method formattable summary.margarita.sim.prob
#' @export
formattable.summary.margarita.sim.prob <-
  function(data, digits=2, backgroundColor=NULL, debug=FALSE){
    data <- as.data.frame(data)

    area <- formattable::area
    formatter <- formattable::formatter
    formattable <- formattable::formattable

    ng <- length(unique(data$groups))
    nrl <- length(unique(data$Exceedance))

    fstring <- paste0("%.", digits, "f")

    if (is.null(backgroundColor)){
      backgroundColor <- shadesOfGrey(ng)
    } else if (length(backgroundColor) != ng){
      stop("backgroundColor should have length to match the number of groups")
    }

    data[, 1:5] <- apply(data[, 1:5], 2, function(X) sprintf(fstring, X))

    data <- mutate(data, Exceedance = factor(data$Exceedance, levels=unique(data$Exceedance))) %>%
      select(Group = groups, `Threshold` = Exceedance, starts_with("Q"))

    areaList <- vector(mode = "list", length = ng * 2)

    rowSeqs <- matrix(1:(ng * nrl), ncol = nrl, byrow=FALSE)
    rowSeqs <- apply(rowSeqs, 1, function(X) paste(X, collapse = ", "))

    for (j in 1:2){ # first is all columns, second is column 5 (Q50%)
      for (i in 1:ng){
        index <- i + (j - 1) * ng
        if (j == 1){
          a <- paste0("area(row=c(", rowSeqs[i], "))")
          ff <- paste0("formatter('span', style = x ~ style('background-color'='", backgroundColor[i], "'))")
          areaList[[index]] <- as.formula(paste(a, "~", ff))
        } else {
          a <- paste0("area(row=c(", rowSeqs[i], "), col=5)")
          ff <- paste0("formatter('span', style = x ~ style('background-color'='", backgroundColor[i], "', 'font-weight'='bold'))")
          areaList[[index]] <- as.formula(paste(a, "~", ff))      }
      }
    }

    if (debug){
      areaList
    } else {
      formattable(data, areaList)
    }
  }
