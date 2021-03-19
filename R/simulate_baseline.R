#' Get range of baseline values and their predicted values
#'
#' @param object An object of class 'margarita'.
#' @param n Number of points in the sequence over baseline.
#'
#' @details The sequence of baseline values is taken on the transformed scale. You
#'   might want to add the option to do it on the observed scale. The sequence of
#'   expected values is on whatever scale was used to fit the regression model.
#'   So baseline and expected might be on different scales.
getLMrange <- function(object, n, range = NULL){

  ## Get baseline data grid, expected values on trans scale
  od <- object[[1]]$data
  baseline <- object$rawBaseline
  newdata <- object$newdata
  
  if (is.null(range)){
    range <- object$trans(range(od[, baseline], na.rm = TRUE))
  } else {
    range <- object$trans(range)
  }
  
  b <- object$invtrans(seq(range[1], range[2], length.out = n))

  newdata <- expand.grid(b, newdata[, 1])
  names(newdata) <- c(baseline, names(object$newdata))

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
                                             Mlabels = NULL, grid.n = 25,
                                             baseline.range = NULL, ...){
  rlm <- object[[1]]
  evmSim <- object[[2]]
  baseline <- object$rawBaseline
  trans <- object$trans
  itrans <- object$invtrans
  summary_alpha <- object$alpha
  newdata <- object$newdata

  if (is.null(Mlabels)){
    Mlabels <- paste0(M / min(M), "xULN")
  }

  summary_alpha <- sort(c(summary_alpha / 2, .5, 1 - summary_alpha / 2))
  family <- evmSim$map$family
  thresholds <- trans(M)

  d <- getLMrange(object, n = grid.n, range = baseline.range)
  par <- predict(evmSim, newdata = d, type = "lp", all = TRUE)$obj$link

  rate <- evmSim$map$rate
  u <- d$expected + evmSim$map$threshold

  ##m <- matrix(rep(thresholds, nrow(d)), ncol = length(thresholds), byrow = TRUE)
  m <- thresholds
  r <- resid(rlm)
  
  out <- lapply(1:grid.n, simBaseline_getProbs,
                u = u, par = par,
                m = m, r = r, p = rate, family = evmSim$map$family,
                th = evmSim$map$threshold)

  out <- lapply(out, function(X){
    o <- lapply(X, function(Z){
      Z <- as.data.frame(apply(Z, 2, simBaseline_sumfun, probs = summary_alpha))
      names(Z) <- Mlabels
      Z <- as.data.frame(t(Z))
      Z$level <- rownames(Z)
      Z
    }) %>% bind_rows() %>% 
      mutate(..group.. = rep(newdata[, 1], each = n() / nrow(newdata)))
  }) %>% bind_rows() %>% 
    mutate(baseline = rep(unique(d[, baseline]), each = length(m) * nrow(newdata)))
  
  names(out)[names(out) == "..group.."] <- object$arm
  rownames(out) <- 1:nrow(out)

  invisible(out)
}

#' @method as.data.frame margarita.sim.baseline.prob
#' @export
as.data.frame.margarita.sim.baseline.prob <- function(x, row.names = NULL, optional = FALSE, ...){
  as.data.frame(unclass(x))
}


# The functions below is ripped off from margarita:::margarita.getProbs and
# margarita:::rp because those functions require par to be a list, expects
# multiple treatment groups and we need to be able to debug
simBaseline_getProbs <- function (X, par, u, p, r, m, family, th){
  ##lapply(1:ncol(m), rp, u = u[X], par = par, p = p,
  lapply(1:length(m), simBaseline_rp, u = u[X], par = par, p = p,
         r = r, xm = m, family = family, th = th)
}

simBaseline_rp <- function (X, xm, u, par, p, r, family, th){
  xm <- xm[X]
  
  res <- sapply(1:length(par), function(ZZ){
    o <- p * (1 - family$prob(xm, par[[ZZ]], list(threshold = u)))
    wh <- u > xm

    if (any(wh)) {
      r <- r[r < quantile(r, 1 - p)]
      o[wh] <- mean((r + u) > xm)
      o[wh] <- p + (1 - p) * o[wh]
    }
    o
  })

  # res is a matrix with length(par), with unique(sapply(par, nrow))
  # estimated probabilities in each column
  res
}

simBaseline_sumfun <- function(x, probs){
  res <- c(mean(x), quantile(x, probs = probs))
  names(res) <- c("Mean", paste0("Q", 100 * probs))
  res
}
