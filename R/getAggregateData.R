#' Take an SDTM type lab dataset and return a dataset ready for extreme value modelling
#' 
#' @param data A \code{data.frame}, usually the lab data, or vital signs data
#' @param subject The unique subject identifier. Defaults to \code{usubjid}
#' @param test The name of the column describing the lab or vital sign test in \code{data}. Defaults to \code{lbtest}
#' @param test.value The value of \code{test} on which to subset the data. Defaults to \code{ALT}
#' @param visit The name of the visit variable. This should be numeric rather than character. Defaults to \code{visitnum}
#' @param baseline.visit The value of \code{visit} corresponding to the baseline visit. Defaults to \code{visit = 1}
#' @param max.study.visit Th evalue of \code{visit} corresponding to the maximum study visit, allowing the user to discard codes representing dropout or follow-up visits of no interest. Defaults to \code{max.study.visit = Inf}
#' @param baseline The name of the column containing the baseline values, if any. Defaults to \code{base}
#' @param value The name of the column containing the test values. Defaults to \code{aval}
#' @param aggregate.fun The function to do the aggregation. Defaults to \code{max}
#' @details Subjects with no post-baseline values are excluded from the output data.
#'          Also, no treatment codes are included in the output; these should be spliced
#'          in from another dataset. If the baseline column does not exist, the function
#'          attempts to use the last observed value at the stated baseline visit.
#'          Also, if a patient has a missing baseline value and no non-missing post-baseline
#'          values, they are silently dropped from the returned data object.
#'          The code was commissioned by AstraZeneca and ought to work with their implementation of SDTM: other implementations may differ.
#' @export
getAggregateData <- function(data, subject="usubjid", test="lbtest", test.value="ALT",
                         visit="visitnum", baseline.visit=1, max.study.visit=Inf,
                         baseline="base", value="aval", aggregate.fun=max){

  data <- data[data[, test] == test.value, ]
  if (nrow(data) == 0) stop("Subsetting on test.value leaves no data")

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
  invisible(na.omit(b))
}

#' Add variables from another data set, matching by patient identifier.
#' @param data The data to have one or more new variables added
#' @param additional.data The dataset containing the new variables to be added to \code{data}
#' @param subject The unique subject identifier. Defaults to \code{subject="usubjid"}
#' @param vars Character vector with the names of the varialbes in \code{additional.data} to be added to \code{data}
#' @param The function adds the first value of \code{var} as matched by \code{subject}. If that's not unique, you want to subset \code{additional.data} to make it uniqure, or sort it appropriately first.
#' @export
addVariables <- function(data, additional.data, subject="usubjid", vars = "trt"){
  for (var in vars){
    data[, var] <- additional.data[match(data[, subject], additional.data[, subject]), var]
  }
  invisible(data)
}

#' Compute study day (days since randomization, screening, or whatever), from date fields
#' @param data A data.frame
#' @param datecol The name of the column containing the dates of the study visits. Defaults to \code{datecol="vsdt"}
#' @param id The name of the column containing the unique subject IDs. Defaults to \code{id="usubjid"}
#' @param visit The name of the columne containing the visit identifiers. Defaults to \code{visit="visit"}
#' @param base.visit The name of the visit taken to be day 0. Defaults to \code{base.visit="Screening"}, but could be numeric
#' @details If any patient has more than one value at the baseline visit, one is taken almost arbitrarily
#' @return The same data.frame but with a new column called "day0" and another called "studyday"
#' @export
getStudyDay <- function(data, datecol="vsdt", id="usubjid", visit="visit", base.visit="Screening"){
  data[, datecol] <- as.Date(data[, datecol], format="%d%b%Y")
  s <- data[data[, visit] == base.visit, ]
  s <- s[order(s[, datecol]), ]
  s <- s[order(s[, id]), c(id, datecol)]
  s <- s[cumsum(rle(s[, id])$lengths), ]
  data$day0 <- s[match(data[, id], s[, id]), datecol]

  data$studyday <- as.numeric(data[, datecol] - data$day0)

  invisible(data)
}

#' Add baseline values of the analysis variable to an analysis dataset
#' @param data A \code{data.frame}
#' @param id The unique patient identifier. Defaults to \code{id = "usubid"}
#' @param visit The name of the visit column. Defaultso to \code{visit="visit"}, but a days on study column could be used instead
#' @param base.visit The value of the visit column at the time the baselne measurements are made
#' @param value The column containing the data values to be analysed
#' @details If multiple measurements are available for some patients at baseline, one is chosen almost arbitrarily
#' @export
getBaselines <- function(data, id="usubjid", visit="visit", base.visit =1, value="aval"){
  data <- data[!is.na(data[, value]), ]
  b <- data[data[, visit] == base.visit, ]
  b <- b[order(b[, id]), c(id, value)]
  b <- b[cumsum(rle(b[, id])$lengths), ]
  
  data$base <- b[match(data[, id], b[, id]), value]
  invisible(data)
}