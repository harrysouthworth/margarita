
#' @method drop1 lmr
#' @export drop1.lmr
drop1.lmr <- 
function (object, scope, ...){
  x <- model.matrix.lm(object)
  asgn <- attr(x, "assign")
  term.labels <- attr(object$terms, "term.labels")
  dfs <- table(asgn[asgn > 0])
  names(dfs) <- term.labels

  if (missing(scope)) 
    scope <- drop.scope(object)
  else {
    if (!is.character(scope)) 
      scope <- attr(terms(update.formula(object, scope)), 
                    "term.labels")
    if (!all(match(scope, term.labels, FALSE))) 
      stop("scope is not a subset of term labels")
  }
  dfs <- dfs[scope]
  k <- length(scope)
cat("meow\n")

  rfpe <- double(k)

  for (i in 1:k) {
    curfrm <- as.formula(paste(".~.-", scope[[i]]))
cat("purr\n")
#browser()
    curobj <- update(object, formula=curfrm)
cat("piss\n")
    rfpe[i] <- RFPE(curobj)
browser()
    if (length(keep)) {
      value[i, 1] <- list(curobj$coefficients)
      value[i, 2] <- list(curobj$fitted)
      value[i, 3] <- list(curobj$residuals)
    }
  }
  
  scope <- c("<none>", scope)
  dfs <- c(0, dfs)
  rfpe <- c(RFPE(object), rfpe)
  dfs[1] <- NA
  aod <- data.frame(Df = dfs, RFPE = rfpe, row.names = scope, 
                    check.names = FALSE)
  head <- c("\nSingle term deletions", "\nModel:", deparse(as.vector(formula(object))))
  if (!missing(scale)) 
    head <- c(head, paste("\nscale: ", format(scale), "\n"))
  oldClass(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  if (length(keep)) {
    value <- value[, keep, drop = FALSE]
    oldClass(value) <- "matrix"
    list(anova = aod, keep = value)
  }
  else aod
}
