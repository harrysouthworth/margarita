#' @method drop1 lmr
#' @export drop1.lmr
drop1.lmr <- function (object, scope, scale=NULL, target="RFPE", k=2, cores=NULL, ...){
  tnm <- target
  aick <- k

  if (target == "RFPE") target <- RFPE.lmr
  else if (target == "AIC") target <- function(x) AIC.lmr(x, k=aick)
  else stop("available targets are 'RFPE' and 'AIC'")
  
  if (is.null(scale)) scale <- object$s
  x <- model.matrix.lm(object)
  asgn <- attr(x, "assign")
  term.labels <- attr(object$terms, "term.labels")
  dfs <- table(asgn[asgn > 0])
  names(dfs) <- term.labels

  if (missing(scope))
    scope <- drop.scope(object)
  else {
    if (!is.character(scope)) 
      scope <- attr(terms(update.formula(object, scope)), "term.labels")
    if (!all(match(scope, term.labels, FALSE)))
      stop("scope is not a subset of term labels")
  }

  dfs <- dfs[scope]
  k <- length(scope)

  fun <- function(X, scope, object, tnm, target){
    curfrm <- as.formula(paste(".~.-", scope[[X]]))
    curobj <- update(object, curfrm, data=object$data, method=object$method, c=object$c)
    if (tnm == "RFPE") target(curobj, scale)
    else target(curobj)
  }

  suppressWarnings(res <- sapply(X=1:k, fun, scope=scope, object=object, tnm=tnm, target=target))

  scope <- c("<none>", scope)
  dfs <- c(0, dfs)
  res <- c(target(object), res) # <------------------------------------------------------------------
  dfs[1] <- NA
  aod <- data.frame(Df = dfs, res = res, row.names = scope, check.names = FALSE)
  names(aod)[2] <- tnm
  head <- c("\nSingle term deletions", "\nModel:", deparse(as.vector(formula(object))))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head

  aod
}

#' Stepwise model selection for a robust linear model
#' @param object An object fit by \code{lmr}.
#' @param scope The lowest model to consider.
#' @param target Character string giving the method for selecting models. Defaults to "RFPE",
#'   the only current alternative being "AIC". 
#' @param direction What direction to take steps in. Only "backwards" is implemented.
#' @param trace Whether or not to report progress. Defaults to \code{trace=TRUE}.
#' @param steps The maximum number of steps to take.
#' @param cores The number of cores to use when running \code{drop1.lmr} in parallel.
#'    Defaults to \code{cores=NULL} and the function will try to guess how many cores
#'    to use.
#' @details The function uses robust finite prediction error to decide when to stop
#'        the selection process. In principle, other approaches could be implemented.
#' @aliases drop1.lmr
#' @export step.lmr
step.lmr <- function (object, scope, target="RFPE", direction = "backward", trace = TRUE, steps = 1000, k=2, cores=NULL, ...){
  if (missing(direction)) direction <- "backward"
  if (direction != "backward") 
    stop("Presently step.lmr only supports backward model selection.")

  tnm <- target
  if (target == "RFPE") target <- RFPE.lmr
  else if (target == "AIC") target <- AIC.lmr
  else stop("target should be either 'RFPE' or 'AIC'")
  
  make.step <- function(models, fit, tnm, object) {
    change <- sapply(models, "[[", "change")
    rdf <- sapply(models, "[[", "df.resid")
    ddf <- c(NA, diff(rdf))
    RES <- sapply(models, "[[", "tnm")
    heading <- c("Stepwise Model Path \nAnalysis of Deviance Table", 
                 "\nInitial Model:", deparse(as.vector(formula(object))), 
                 "\nFinal Model:", deparse(as.vector(formula(fit))), 
                 "\n")
    aod <- data.frame(Step = change, Df = ddf, `Resid. Df` = rdf, RES = RES, check.names=FALSE)
    names(aod)[4] <- tnm
    attr(aod, "heading") <- heading
    class(aod) <- c("anova", "data.frame")
    fit$anova <- aod
    fit
  }

  if (missing(scope)) {
    fdrop <- numeric(0)
    fadd <- NULL
  }
  else {
    if (is.list(scope)) {
      fdrop <- if (!is.null(fdrop <- scope$lower)) 
        attr(terms(update.formula(object, fdrop)), "factor")
      else numeric(0)
      fadd <- if (!is.null(fadd <- scope$upper)) 
        attr(terms(update.formula(object, fadd)), "factor")
    }
    else {
      fadd <- if (!is.null(fadd <- scope)) 
        attr(terms(update.formula(object, scope)), "factor")
      fdrop <- numeric(0)
    }
  }

  m <- model.frame(object)
  obconts <- object$contrasts
  objectcall <- object$call

  Terms <- object$terms
  x <- model.matrix(Terms, m, contrasts = obconts)

  Asgn <- attr(x, "assign")
  term.labels <- attr(Terms, "term.labels")
  a <- attributes(m)
  y <- model.extract(m, "response")
  w <- model.extract(m, "weights")
  if (is.null(w)) w <- rep(1, nrow(m))
  models <- vector("list", steps)

  n <- length(object$fitted)
  scale <- object$s
  fit <- object
  bRES <- target(fit)
  nm <- 1
  Terms <- fit$terms
  if (trace) 
    cat("Start:  RFPE=", format(round(bRES, 4)), "\n", deparse(as.vector(formula(fit))), "\n\n")
  models[[nm]] <- list(df.resid = fit$df.residual, change = "", tnm = bRES)
  RES <- bRES + 1
  
  while (bRES < RES & steps > 0) {
    steps <- steps - 1
    RES <- bRES
    bfit <- fit
    ffac <- attr(Terms, "factor")
    scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
    aod <- NULL
    change <- NULL

    if (ndrop <- length(scope$drop)) {
      aod <- drop1.lmr(fit, scope$drop, scale, target=tnm, k=k, cores=cores)
      if (trace) print(aod)
      change <- rep("-", ndrop + 1)
    }

    if (is.null(aod)) break
    o <- order(aod[, tnm])[1]
    if (o[1] == 1) break

    change <- paste(change[o], dimnames(aod)[[1]][o])
    Terms <- terms(update(formula(fit), eval(parse(text = paste("~ .", change)))))
    attr(Terms, "formula") <- new.formula <- formula(Terms)
    newfit <- lmr(new.formula, data = object$data, c=object$c, method="MM")

    bRES <- aod[, tnm][o]

    if (trace) cat("\nStep:  RFPE =", format(round(bRES, 4)), "\n", deparse(as.vector(formula(Terms))), "\n\n")
    if (bRES >= RES) break

    nm <- nm + 1
    models[[nm]] <- list(df.resid = newfit$df.resid, change = change, tnm = bRES)
    fit <- c(newfit, list(formula = new.formula))
    oc <- objectcall
    oc$formula <- as.vector(fit$formula)
    fit$call <- oc
    class(fit) <- class(object)
  }
  make.step(models = models[seq(nm)], fit, tnm, object)
}
