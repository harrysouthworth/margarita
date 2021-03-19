library(tidyverse)
library(robustbase)
library(robust)
library(margarita)

test_that("getLMrange behaves as expected", {
  liver$ndose <- as.numeric(liver$dose)

  rmod <- lmr(log(ALT.M) ~ log(ALT.B) + ndose, data = liver)
  liver$r <- resid(rmod)
  bmod <- evm(r, data = liver, qu = .7, xi = ~ndose, method = "sim")

  mar <- margarita(rmod, bmod, newdata = data.frame(ndose = 1:4),
                   arm = "ndose", baseline = "ALT.B")

  for (i in 1:5){
    n <- sample(10:25, size = 1)

    o <- getLMrange(mar, n = n)

    expect_equal(unique(table(o$ndose)), n,
                 label = "getLMrange: ndose repeats correctly")

    expect_equal(nrow(o), n * length(unique(liver$ndose)))

    expect_equal(min(liver$ALT.B), o$ALT.B[1],
                 label = "getLMrange: min ALT.B matches data")
    expect_equal(max(liver$ALT.B), o$ALT.B[nrow(o)],
                 label = "getLMrange: max ALT.B matches data")

    br <- c(sample(4:10, size = 1), sample(30:60, size = 1))

    o <- getLMrange(mar, n = n, range = br)

    expect_equal(unique(table(o$ndose)), n,
                 label = "getLMrange: specifying range results in ndose repeats correctly")
    expect_equal(nrow(o), n * length(unique(liver$ndose)),
                 label = "getLMrange: specifying range results in correct size output")
    expect_equal(br[1], o$ALT.B[1],
                 label = "getLMrange: specifying range results in correct minimum ALT.B")
    expect_equal(br[2], o$ALT.B[nrow(o)],
                 label = "getLMrange: specifying range results in correct maximum ALT.B")


    co <- coef(rmod)

    ealt <- co[1] + co[2] * log(o$ALT.B) + co[3] * o$ndose
    expect_equal(cor(ealt, o$expected), 1,
                 label = "getLMrange: expected values are correct")
  }
})


test_that("simBaseline_rp does what it ought", {
  liver$ndose <- as.numeric(liver$dose)

  rmod <- lmr(log(ALT.M) ~ log(ALT.B) + ndose, data = liver)
  liver$r <- resid(rmod)
  bmod <- evm(r, data = liver, qu = .7, xi = ~ndose, method = "sim")

  mar <- margarita(rmod, bmod, newdata = data.frame(ndose = 1:4),
                   arm = "ndose", baseline = "ALT.B")


  n <- sample(10:25, size = 1)

  br <- c(sample(4:10, size = 1), sample(30:60, size = 1))

  o <- getLMrange(mar, n = n, range = br)

  r <- resid(rmod)
  u <- o$expected + bmod$map$threshold

  # pd <- cbind(o, u)
  # ggplot(pd, aes(log(ALT.B), expected)) +
  #   geom_point() +
  #   facet_wrap(~ndose) +
  #   geom_point(aes(y = u, color = "blue"))

  par <- predict(bmod, newdata = o, type = "lp", all = TRUE)$obj$link

  simBaseline_rp(1, log(3 * 36), )




})

