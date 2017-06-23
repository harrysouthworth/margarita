context("P(excess)")

test_that("pExcessByBaseline output behaves as expected", {
  d <- liver

  ##############################################################################
  ### TREATMENT AS A FACTOR
  ###
  rmod <- lmr(log(ALT.M) ~ log(ALT.B) + dose, data=d)
  d$r <- resid(rmod)
  th <- quantile(d$r, .5)
  gmod <- evm(r, data=d, th=th, xi=~dose, method="sim", verbose=FALSE)
  nd <- data.frame(dose=LETTERS[1:4])
  m <- margarita(rmod, gmod, newdata=nd, baseline="ALT.B")

  p3 <- pExcessByBaseline(m, M=3*36, n=200)
  p3.25 <- pExcessByBaseline(m, M=3*36, n=25)

  expect_that(nrow(p3), equals(nrow(nd) * 200), label="pExcessByBaseline: output data size, n=200")
  expect_that(nrow(p3.25), equals(nrow(nd) * 25), label="pExcessByBaseline: output data size, n=25")

  expect_that(range(d$ALT.B), equals(range(p3$ALT.B)), label="pExcessByBaseline: baseline ranges, n=200")
  expect_that(range(d$ALT.B), equals(range(p3.25$ALT.B)), label="pExcessByBaseline: baseline ranges ok, n=25")

  expect_more_than(cor(sort(fitted(rmod)), sort(sample(p3$threshold, size=nrow(d), replace=TRUE))), .90,
                   label="pExcessByBaseline: correlation of fitted values with thresholds, n=200")
  expect_more_than(cor(sort(fitted(rmod)), sort(sample(p3.25$threshold, size=nrow(d), replace=TRUE))), .90,
                   label="pExcessByBaseline: correlation of fitted values with thresholds, n=25")

  expect_more_than(cor(p3$mean, p3$".5"), .95,
                   label="pExcessByBaseline: means and medians correlated, n=200")
  expect_more_than(cor(p3.25$mean, p3.25$".5"), .95,
                   label="pExcessByBaseline: means and medians correlated, n=25")

  expect_true(min(p3[, -c(1:3)]) >= 0, label="pExcessByBaseline: probabilities >= 0, n=200")
  expect_true(min(p3.25[, -c(1:3)]) >= 0, label="pExcessByBaseline: probabilities >= 0, n=25")

  expect_true(max(p3[, -c(1:3)]) <= 1, label="pExcessByBaseline: probabilities <= 1, n=200")
  expect_true(max(p3.25[, -c(1:3)]) <= 1, label="pExcessByBaseline: probabilities <= 1, n=25")

  sp3 <- split(p3, p3$dose)
  sp3.25 <- split(p3.25, p3.25$dose)

  # THE FOLLOWING WILL ONLY WORK FOR HIGH THREHSOLDS. If simulation error comes into it, would expect this to fail
  for (i in 1:length(sp3)){
    expect_true(min(diff(sp3[[i]]$mean)) >= 0, label="pExcessByBaseline: mean probabilities are monotonic increasing with baseline, n=200")
    expect_true(min(diff(sp3[[i]]$".5")) >= 0, label="pExcessByBaseline: median probabilities are monotonic increasing with baseline, n=200")

    expect_true(min(diff(sp3.25[[i]]$mean)) >= 0, label="pExcessByBaseline: mean probabilities are monotonic increasing with baseline, n=25")
    expect_true(min(diff(sp3.25[[i]]$".5")) >= 0, label="pExcessByBaseline: median probabilities are monotonic increasing with baseline, n=25")
  }

  mp3 <- as.matrix(p3[, 5:9])      # Stops apply coercing to character
  mp3.25 <- as.matrix(p3[, 5:9])
  wh <- apply(mp3, 1, function(X) expect_true(min(diff(X)) >= 0, label="pExcessByBaseline: quantiles are monotonic increasing, n=200"))
  wh <- apply(mp3.25, 1, function(X) expect_true(min(diff(X)) >= 0, label="pExcessByBaseline: quantiles are monotonic increasing, n=25"))

  ##############################################################################
  ### TREATMENT AS LINEAR
  ###
  d$ndose <- as.numeric(d$dose)
  rmod <- lmr(log(ALT.M) ~ log(ALT.B) + ndose, data=d)
  d$r <- resid(rmod)
  th <- quantile(d$r, .5)
  gmod <- evm(r, data=d, th=th, xi=~ndose, method="sim", verbose=FALSE)
  nd <- data.frame(ndose=1:4)
  m <- margarita(rmod, gmod, newdata=nd, baseline="ALT.B")
  
  p3 <- pExcessByBaseline(m, M=3*36, n=200)
  p3.25 <- pExcessByBaseline(m, M=3*36, n=25)

  expect_that(nrow(p3), equals(nrow(nd) * 200), label="pExcessByBaseline: output data size, n=200, linear dose")
  expect_that(nrow(p3.25), equals(nrow(nd) * 25), label="pExcessByBaseline: output data size, n=25, linear dose")

  expect_that(range(d$ALT.B), equals(range(p3$ALT.B)), label="pExcessByBaseline: baseline ranges, n=200, linear dose")
  expect_that(range(d$ALT.B), equals(range(p3.25$ALT.B)), label="pExcessByBaseline: baseline ranges ok, n=25, linear dose")

  expect_more_than(cor(sort(fitted(rmod)), sort(sample(p3$threshold, size=nrow(d), replace=TRUE))), .90,
                   label="pExcessByBaseline: correlation of fitted values with thresholds, n=200, linear dose")
  expect_more_than(cor(sort(fitted(rmod)), sort(sample(p3.25$threshold, size=nrow(d), replace=TRUE))), .90,
                   label="pExcessByBaseline: correlation of fitted values with thresholds, n=25, linear dose")

  expect_more_than(cor(p3$mean, p3$".5"), .95,
                   label="pExcessByBaseline: means and medians correlated, n=200, linear dose")
  expect_more_than(cor(p3.25$mean, p3.25$".5"), .95,
                   label="pExcessByBaseline: means and medians correlated, n=25, linear dose")

  expect_true(min(p3[, -c(1:3)]) >= 0, label="pExcessByBaseline: probabilities >= 0, n=20, linear dose0")
  expect_true(min(p3.25[, -c(1:3)]) >= 0, label="pExcessByBaseline: probabilities >= 0, n=25, linear dose")

  expect_true(max(p3[, -c(1:3)]) <= 1, label="pExcessByBaseline: probabilities <= 1, n=200, linear dose")
  expect_true(max(p3.25[, -c(1:3)]) <= 1, label="pExcessByBaseline: probabilities <= 1, n=25, linear dose")

  sp3 <- split(p3, p3$ndose)
  sp3.25 <- split(p3.25, p3.25$ndose)

  # THE FOLLOWING WILL ONLY WORK FOR HIGH THREHSOLDS. If simulation error comes into it, would expect this to fail
  for (i in 1:length(sp3)){
    expect_true(min(diff(sp3[[i]]$mean)) >= 0, label="pExcessByBaseline: mean probabilities are monotonic increasing with baseline, n=200, linear dose")
    expect_true(min(diff(sp3[[i]]$".5")) >= 0, label="pExcessByBaseline: median probabilities are monotonic increasing with baseline, n=200, linear dose")

    expect_true(min(diff(sp3.25[[i]]$mean)) >= 0, label="pExcessByBaseline: mean probabilities are monotonic increasing with baseline, n=25, linear dose")
    expect_true(min(diff(sp3.25[[i]]$".5")) >= 0, label="pExcessByBaseline: median probabilities are monotonic increasing with baseline, n=25, linear dose")
  }

  mp3 <- as.matrix(p3[, 5:9])      # Stops apply coercing to character
  mp3.25 <- as.matrix(p3[, 5:9])
  wh <- apply(mp3, 1, function(X) expect_true(min(diff(X)) >= 0, label="pExcessByBaseline: quantiles are monotonic increasing, n=200, linear dose"))
  wh <- apply(mp3.25, 1, function(X) expect_true(min(diff(X)) >= 0, label="pExcessByBaseline: quantiles are monotonic increasing, n=25, linear dose"))
  })
