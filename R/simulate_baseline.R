#' Get range of baseline values and their predicted values
#'
#' @param object An object of class 'margarita'.
#' @param n Number of points in the sequence over baseline.
#'
#' @details The sequence of baseline values is taken on the transformed scale. You
#'   might want to add the option to do it on the observed scale. The sequence of
#'   expected values is on whatever scale was used to fit the regression model.
#'   So baseline and expected might be on different scales.
getLMrange <- function(object, n){

  ## Get baseline data grid, expected values on trans scale
  od <- object[[1]]$data
  b <- exp(seq(min(object$trans(od[, object$rawBaseline])), max(object$trans(od[, object$rawBaseline])), length.out = n))
  newdata <- expand.grid(b, object$newdata[, 1])
  names(newdata) <- c(object$rawBaseline, names(object$newdata))

  p <- suppressWarnings(predict.lm(object[[1]], newdata = newdata))

  cbind(newdata, expected = p)
}


#' Get return levels over observed range of baseline values
#'
#' @param object An object of class 'margarita'.
#' @param baseline String identifying the baseline variable.
simulate.margarita.baseline.rl <- function(object, M, nsim = 1, seed = NULL, grid.n = 25, ...){

  ## Get contents of (), if they exist
  ##if (substring(wh, nchar(wh)) == ")"){
  ##  baseline <- stringr::str_extract_all(object$baseline, "(?<=\\().+?(?=\\))")[[1]]
  ##}
  baseline <- object$rawBaseline

  alpha <- object$alpha
  rlm <- object[[1]]
  evmSim <- object[[2]]
  trans <- object$trans
  itrans <- object$invtrans


  out <- getLMrange(object, n = grid.n)

  p <- predict(evmSim, M = M, ci.fit = TRUE, alpha = alpha)$obj
  nms <- gsub("%", "", paste0("Q", colnames(p[[1]])))

  p <- unlist(p) %>%
    matrix(nrow = length(M), byrow = TRUE) %>%
    as.data.frame() %>%
    setNames(c("Mean", nms[-1])) %>%
    mutate(M = M)

  res <- list()
  for (m in M){
    om <- apply(p[p$M == m, !(names(p) == "M")], 2, function(X) itrans(X + out$expected)) %>%
      as.data.frame() %>%
      bind_cols(out) %>%
      mutate(M = paste0(m, "-subject return level"),
             tooltip = paste("Baseline:", round(baseline), "\nReturn level:", round(Q50)))

    res[[paste("M =", m)]] <- om
  }

  bind_rows(res) %>%
    mutate(M = factor(M, levels = unique(M)))
}

simulate.margarita.baseline.prob <- function(object, nsim = 1, seed = NULL, M,
                                             Mlabels = NULL, grid.n = 25, ...){

  rlm <- object[[1]]
  evmSim <- object[[2]]
  baseline <- object$rawBaseline
  trans <- object$trans
  itrans <- object$invtrans
  alpha <- object$alpha
  newdata <- object$newdata

  if (is.null(Mlabels)){
    Mlabels <- paste0(M / min(M), "xULN")
  }

  alpha <- sort(c(alpha / 2, .5, 1 - alpha / 2))
  family <- evmSim$map$family
  thresholds <- trans(M)

  d <- getLMrange(object, n = grid.n)
  par <- predict(evmSim, type = "lp", all = TRUE)$obj$link[[1]]

  rate <- evmSim$map$rate
  u <- d$expected + evmSim$map$threshold

  m <- matrix(rep(thresholds, nrow(d)), ncol = length(thresholds), byrow = TRUE)
  r <- resid(rlm)

  # The functions below is ripped off from margarita:::margarita.getProbs and
  # margarita:::rp because those function requires par to be a list, expects
  # multiple treatment groups and we need to be able to debug
  getProbs <- function (X, par, u, p, r, m, family, th){
    lapply(1:ncol(m), rp, u = u[X], par = par, p = p,
           r = r, xm = m, family = family, th = th)
  }

  rp <- function (X, xm, u, par, p, r, family, th){
    xm <- xm[, X]
    res <- p * (1 - family$prob(xm, par, list(threshold = u)))
    wh <- u > xm

    if (any(wh)) {
      r <- r[r < quantile(r, 1 - p)]
      res[wh] <- mean((r + u) > m[1])
      res[wh] <- p + (1 - p) * res[wh]
    }
    res
  }

  sumfun <- function(x){
    res <- c(mean(x), quantile(x, probs = alpha))
    names(res) <- c("Mean", paste0("Q", 100 * alpha))
    res
  }

  out <- lapply(1:nrow(d), getProbs,
                u = u, par = par,
                m = m, r = r, p = rate, family = evmSim$map$family,
                th = evmSim$map$threshold) %>%
    lapply(function(X){
      matrix(unlist(X), ncol = length(thresholds), byrow=FALSE)
    })

  # This next block could be piped from the previous one, but it's kept
  # separate to enable debugging at this point.
  out <-  lapply(out, function(X){
    apply(X, 2, sumfun) %>% t() %>%
      as.data.frame() %>%
      mutate(threshold = Mlabels,
             threshold = factor(threshold, levels = Mlabels))
  }) %>% bind_rows() %>%
    mutate(baseline = rep(d[, baseline], each = length(thresholds)),
           ..group.. = rep(newdata[, 1], length.out = n()))

  names(out)[names(out) == "..group.."] <- object$arm

  invisible(out)
}

#' @method as.data.frame margarita.sim.baseline.prob
#' @export
as.data.frame.margarita.sim.baseline.prob <- function(x, row.names = NULL, optional = FALSE, ...){
  as.data.frame(unclass(x))
}
