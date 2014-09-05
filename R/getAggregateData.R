#' Take an SDTM type lab dataset and return a dataset ready for extreme value modelling
#' 
#' @param data A \code{data.frame}, usually the lab data, or vital signs data
#' @param subject The unique subject identifier. Defaults to \code{usubjid}
#' @param test The name of the column describing the lab or vital sign test in \code{data}. Defaults to \code{lbtest}
#' @param test.value The value of \code{test} on which to subset the data. Defaults to \code{ALT}
#' @param visit The name of the visit variable. This should be numeric rather than character. Defaults to \code{visitnum}
#' @param baseline.visit The value of \code{visit} corresponding to the baseline visit. Defaults to \code{visit = 1}
#' @param max.study.visit Th evalue of \code{visit} corresponding to the maximum study visit, allowing the user to discard codes representing dropout or follow-up visits of no interest. Defaults to \code{max.study.visit = 20}
#' @param baseline The name of the column containing the baseline values, if any. Defaults to \code{base}
#' @param value The name of the column containing the test values. Defaults to \code{aval}
#' @param aggregate.fun The function to do the aggregation. Defaults to \code{max}
#' @details Subjects with no post-baseline values are excluded from the output data. Also, no treatment codes are included in the output; these should be spliced in from another dataset. If the baseline column does not exist, the function attempts to use the last observed value at the stated baseline visit.\n The code was commissioned by AstraZeneca and ought to work with their implementation of SDTM: other implementations may differ.
#' @export
getAggregateData <- function(data, subject="usubjid", test="lbtest", test.value="ALT (SGOT)",
                         visit="visitnum", baseline.visit=1, max.study.visit=90,
                         baseline="base", value="aval", aggregate.fun=max){

  data <- data[data[, test] == test.value, ]

  # Get baseline data
  b <- data[data[, visit] == baseline.visit, ]

  if (baseline %in% names(data))
    b <- b[, c(subject, baseline)]
  else{ # baseline is NULL
    b <- b[order(b[, visit]), ]
    b <- b[order(b[, subject]), ]
    b <- b[cumsum(rle(b[, subject])$lengths), c(subject, value)]
    names(b)[names(b) == value] <- baseline
  }

  # Get maximum post-baseline values
  data <- data[data[, visit] > baseline.visit & data[, visit] <= max.study.visit, ]
  dr <- na.omit(data[, c(subject, value)]) # avoid warnings in next line due to pats with no data
  m <- aggregate(dr[, value], by=list(dr[, subject]), FUN=aggregate.fun, na.rm=TRUE)

  b <- b[b[, subject] %in% m$Group.1, ]
  b[, value] <- m$x[match(b[, subject], m$Group.1)]
  b <- b[b[, value] > -Inf, ]
  invisible(b)
}