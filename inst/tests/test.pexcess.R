context("P(excess)")

test_that("Baseline values are produced correctly", {
  d <- liver
  rmod <- lmr(log(ALT.M) ~ log(ALT.B) + dose, data=d)
  d$r <- resid(rmod)
  gmod <- evm(r, data=d, qu=.5, xi=~dose, method="sim", verbose=FALSE)
  
  
  
})