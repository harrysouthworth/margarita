context("plot functions")

test_that("ggplot.lmr output matches manually created plots", {
  m <- lmr(log(ALT.M) ~ log(ALT.B) + as.numeric(dose), data=liver)
  r <- resid(m); f <- fitted(m); o <- log(liver$ALT.M)
  
  par(mfrow=c(2, 2))
  qqnorm(r); qqline(r)
  hist(r, prob=TRUE, ylim=c(0, 1.5)); rug(r); lines(density(r))
  plot(f, r); lines(lowess(f, r))
  plot(f, o); lines(lowess(f, o))
  ggplot(m)
})

test_that("ggplot.gpdThresh output matches plot.gpdThresh output", {
  m <- lmr(log(ALT.M) ~ log(ALT.B) + as.numeric(dose), data=liver)
  liver$r <- resid(m)
  thresh <- gpdThresh(liver$r)
  par(mfrow=c(2, 2))
  plot(gpdRangeFit(liver$r))
  plot(mrl(liver$r))
  plot(egp3RangeFit(liver$r))
  ggplot(thresh)
})

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

test_that("ggplot.evmSim output matches plot.evmSim output", {
  m <- lmr(log(ALT.M) ~ log(ALT.B) + as.numeric(dose), data=liver)
  liver$r <- resid(m)
  em <- evm(r, data=liver, qu=.5, xi=~as.numeric(dose), method="sim", iter=10500, verbose=FALSE)
  par(mfrow=c(3, 3))
  plot(em)
  ggplot(em)  
})
