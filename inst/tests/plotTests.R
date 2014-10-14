context("plot functions")

test_that("ggplot.evmOpt output matches plot.evmOpt output", {
  m <- lmr(log(ALT.M) ~ log(ALT.B) + as.numeric(dose), data=liver)
  liver$r <- resid(m)
  em <- evm(r, data=liver, qu=.5)
  par(mfrow=c(2, 2))
  plot(em)
  ggplot(em)
  
  em <- evm(r, data=liver, qu=.5, xi=~as.numeric(dose))
  par(mfrow=c(2, 2))
  plot(em)
  ggplot(em)
})

test_that("ggplot.gpdThresh output matches plot.gpdThresh output", {
  m <- lmr(log(ALT.M) ~ log(ALT.B) + as.numeric(dose), data=liver)
  liver$r <- resid(m)
  em <- evm(r, data=liver, qu=.5)
  par(mfrow=c(2, 2))
  plot(em)
  ggplot(em)
})


