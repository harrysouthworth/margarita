#' Take an SDTM type lab dataset and return a dataset ready for extreme value modelling
#' 
#' @param data A \code{data.frame}, usually the lab data, or vital signs data
#' @param subject The unique subject identifier. Defaults to \code{usubjid}
#' @param visit The name of the visit variable. This should be numeric rather than character. Defaults to \code{visitnum}
#' @param baseline.visit The value of \code{visit} corresponding to the baseline visit. Defaults to \code{visit = 1}
#' @param max.study.visit Th evalue of \code{visit} corresponding to the maximum study visit, allowing the user to discard codes representing dropout or follow-up visits of no interest. Defaults to \code{max.study.visit = Inf}
#' @param baseline The name of the column containing the baseline values, if any. Defaults to \code{base}
#' @param value The name of the column containing the test values. Defaults to \code{aval}
#' @param aggregate.fun The function to do the aggregation. Defaults to \code{max}
#' @details Subjects with no post-baseline values are excluded from the output data.
#'          Also, no treatment codes are included in the output; these should be spliced
#'          in from another dataset.
#'          Also, if a patient has a missing baseline value or no non-missing post-baseline
#'          values, they are silently dropped from the returned data object.
#'          The code was commissioned by AstraZeneca and ought to work with their
#'          implementation of SDTM: other implementations may differ.
#' @examples alt <- getBaselines(lb[lb$lbtestcd == "ALT", ], visit="visitnum")
#'           alt <- getAggregateData(alt, visit="visitnum")
#'           alt <- addVariables(alt, dm, vars="armcd")
#' @export
getAggregateData <- function(data, subject="usubjid",
                         visit="visitnum", baseline.visit=1, max.study.visit=Inf,
                         baseline="baseline", value="aval", aggregate.fun=max){
  if (nrow(data) == 0) stop("data has 0 rows")
  if (baseline %in% names(data)){
    b <- data[, c(subject, baseline)]
    b <- b[!(is.na(b[, subject]) | is.na(b[, baseline])), ]
    b <- b[order(b[, subject]), ]
    b <- b[cumsum(rle(as.character(b[, subject]))$lengths), ]
  } else stop ("baseline variable is not in the data")

  # Get maximum post-baseline values
  data <- data[data[, visit] > baseline.visit & data[, visit] <= max.study.visit, ]
  dr <- na.omit(data[, c(subject, value)]) # avoid warnings in next line due to pats with no data
  m <- aggregate(dr[, value], by=list(dr[, subject]), FUN=aggregate.fun, na.rm=TRUE)

  b <- b[b[, subject] %in% m$Group.1, ] # Dropping cases where there is no maximum
  b[, value] <- m$x[match(b[, subject], m$Group.1)]
  b <- b[b[, value] > -Inf, ]
  invisible(na.omit(b))
}

#' Add variables from another data set, matching by patient identifier.
#' @param data The data to have one or more new variables added
#' @param additional.data The dataset containing the new variables to be added to \code{data}
#' @param subject The unique subject identifier. Defaults to \code{subject="usubjid"}
#' @param vars Character vector with the names of the varialbes in \code{additional.data} to be added to \code{data}
#' @details The function adds the first value of \code{var} as matched by \code{subject}. If
#'          that's not unique, you want to subset \code{additional.data} to make it unique,
#'          or sort it appropriately first.
#'          An alternative is to use the \code{join} function in the \code{plyr}
#'          package. NOTE that the call should
#'          be like \code{mydata <- addVariables(mydata, otherdata, vars="trt")}:
#'          that is, the returned object is
#'          a \code{data.frame}, NOT a new column to be appended to \code{myhdata}.
#' @examples lb <- addVariables(lb, dm, vars="armcd")
#' @export
addVariables <- function(data, additional.data, subject="usubjid", vars = "trt"){
  m <- match(data[, subject], additional.data[, subject])
  for (var in vars){
    if (var %in% names(data)){
      if (nrow(data) > sum(is.na(data[, var]))){
        warning(paste(var, "already in data and it isn't empty"))
      } else {
        data[, var] <- additional.data[m, var]
      }
    } else {
      data[, var] <- additional.data[m, var]
    }
  }
  invisible(data)
}

#' Compute study day (days since randomization, screening, or whatever), from date fields
#' @param data A data.frame
#' @param datecol The name of the column containing the dates of the study visits. Defaults to \code{datecol="vsdt"}
#' @param subject The name of the column containing the unique subject IDs. Defaults to \code{id="usubjid"}
#' @param visit The name of the columne containing the visit identifiers. Defaults to \code{visit="visit"}
#' @param baseline.visit The name of the visit taken to be day 0. Defaults to \code{baseline.visit=1}, but could be numeric
#' @details If any patient has more than one value at the baseline visit, one is taken almost arbitrarily
#' @return The same data.frame but with a new column called "day0" and another called "studyday"
#' @export
getStudyDay <- function(data, datecol="vsdt", subject="usubjid", visit="visit", baseline.visit=1){
  data[, datecol] <- as.Date(data[, datecol], format="%d%b%Y")
  s <- data[data[, visit] == baseline.visit, ]
  s <- s[order(s[, datecol]), ]
  s <- s[order(s[, subject]), c(subject, datecol)]
  s <- s[cumsum(rle(s[, subject])$lengths), ]
  data$day0 <- s[match(data[, subject], s[, subject]), datecol]

  data$studyday <- as.numeric(data[, datecol] - data$day0)

  invisible(data)
}

#' Add baseline values of the analysis variable to an analysis dataset
#' @param data A \code{data.frame}
#' @param subject The unique patient identifier. Defaults to \code{id = "usubid"}
#' @param visit The name of the visit column. Defaultso to \code{visit="visit"}, but
#'        a days on study column could be used instead
#' @param baseline.visit The value of the visit column at the time the baselne
#'        measurements are made
#' @param baseline The name of the baseline column. If it is already in the data and not
#'        empty, the function will exit with a warning. Otherwise, this is the
#'        name of the baseline variable in the output dataset.
#' @param value The column containing the data values to be analysed
#' @details If multiple measurements are available for some patients at baseline,
#'          one is chosen almost arbitrarily.
#' @export
getBaselines <- function(data, subject="usubjid", visit="visit", baseline.visit =1,
                         value="aval", baseline="baseline"){
  if (nrow(data) == 0) stop("data has 0 rows")
  if (baseline %in% names(data)){
    if (nrow(data) > sum(is.na(data[, baseline]))){
      warning("data contains the baseline column already, and it isn't empty")
      invisible(data)
    }
  }

  b <- data[data[, visit] == baseline.visit, ]
#  b <- b[!is.na(b[, value]), ]
  b <- b[order(b[, subject]), c(subject, value)]
  b <- b[cumsum(rle(as.character(b[, subject]))$lengths), ]

  data[, baseline] <- b[match(data[, subject], b[, subject]), value]
  invisible(data)
}
