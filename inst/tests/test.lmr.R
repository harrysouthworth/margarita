context("lmr functions")

test_that("lmr and lmRob give similar results", {
  wh <- try(library(robust), silent=TRUE)
  if (class(wh) == "try-error") return()

  r1 <- lmr(log(ALT.M) ~ log(ALT.B) + as.numeric(dose), data=liver)
  r2 <- lmRob(log(ALT.M) ~ log(ALT.B) + as.numeric(dose), data=liver,
              control=lmRob.control(weight=c("Bisquare", "Bisquare"),
                                    efficiency=.85))

  co1 <- coef(summary(r1))
  co2 <- coef(summary(r2))

  expect_equal(co1[, 1], co2[, 1], tol=.001, label="lmr: point estimates")
  expect_equal(co1[, 2], co2[, 2], tol=.05, label="lmr: standard errors")
  expect_equal(r1$c, r2$yc, label="lmr: tuning constant")
  expect_equal(resid(r1), resid(r2), tol=.001, label="lmr: residuals")
})

test_that("RFPE works", {
  rfpe1 <- RFPE.lmr(r1)
  rfpe2 <- RFPE.lmRob(r2)
  expect_equal(rfpe1, rfpe2, tol=.001, lable="lmr: RFPE for lmr and lmRob")
  expect_equal(rfpe2 * nrow(liver), lmRob.RFPE(r2), tol=.001, lable="lmr: RFPE for lmRob")
  
  fun <- function(){
    d <- liver[sample(1:nrow(liver), size=nrow(liver), replace=TRUE), ]
    r1 <- lmr(log(ALT.M) ~ log(ALT.B) + as.numeric(dose), data=d)
    r2 <- lmRob(log(ALT.M) ~ log(ALT.B) + as.numeric(dose), data=d,
                     control=lmRob.control(weight=c("Bisquare", "Bisquare"),
                                           efficiency=.85))
    c(RFPE.lmr(r1), RFPE.lmRob(r2))
  }
  bs <- t(replicate(100, fun()))
  expect_more_than(cor(bs[, 1], bs[, 2]), .99, label="lmr: lmr vs lmRob with RFPE")
})

test_that("bisquare works", {
  # Since RFPE is ok, bisquare must be. Test it directly, anyway
  x <- rnorm(1000)
  expect_equal(cor(rho.weight(x, 2, 3.443689), margarita:::bisquare(x, d=0)), 1, label="lmr: bisquare rho")
  expect_equal(cor(psi.weight(x, 2, 3.443689), margarita:::bisquare(x, d=1)), 1, label="lmr: bisquare psi")
  expect_equal(cor(psp.weight(x, 2, 3.443689), margarita:::bisquare(x, d=2)), 1, label="lmr: bisquare dpsi")
  # Note MASS implementation doesn't multiply by input
  expect_equal(cor(psi.bisquare(x, c=3.443689)*x, margarita:::bisquare(x, d=1)), 1, label="lmr: bisquare MASS")
  expect_equal(cor(psi.bisquare(x, c=3.443689, deriv=1), margarita:::bisquare(x, d=2)), 1, label="lmr: bisquare MASS")
})

test_that("drop1.lmr works", {
  m <- lmr(log(ALT.M) ~ log(ALT.B) + as.numeric(dose), data=liver)
  d <- drop1(m)
  
  s <- m$s
  rf <- c(RFPE(m), RFPE(lmr(log(ALT.M) ~ as.numeric(dose), data=liver), scale=s),
          RFPE(lmr(log(ALT.M) ~ log(ALT.B), data=liver), scale=s))

  expect_equal(d$RFPE, rf)
})

test_that("step.lmr works", {
  d <- liver
  d <- cbind(d, data.frame(rmvnorm(nrow(d), mean=rep(0, 4))))
  m <- lmr(log(ALT.M) ~ log(ALT.B) + as.numeric(dose) + X1 + X2 + X3 + X4, data=d)
  step.lmr(m)
  
})