if (FALSE){
test_that("getCIquantiles behaves as expected", {
  expect_that(getCIquantiles(.05), equals(c(0.025, 0.50, 0.975)))
  expect_that(getCIquantiles(.1), equals(c(0.050, 0.50, 0.950)))
  expect_that(getCIquantiles(c(.1, .5)), equals(c(.050, .250, .50, .750, .950)))
  expect_that(getCIquantiles(c(.5, .1)), equals(sort(getCIquantiles(c(.5, .1)))))
  expect_that(getCIquantiles(0), throws_error())
  expect_that(getCIquantiles(1), throws_error())
  expect_that(getCIquantiles(c(.05, .1, .5)), throws_error())
})

test_that("Segments are returned correctly", {
  data <- data.frame(rmvnorm(100, rep(0, 3)))

  expect_that(getSegmentData(data), throws_error())

  data$group <- rep(LETTERS[1:2], each=50)

  expect_that(getSegmentsData(data), throws_error())
  
  data[, 1:3] <- t(apply(data[, 1:3], 1, sort))

  expect_that(getSegmentData(data)[[1]], is_equivalent_to(data[, c("X1", "X3", "group")]))
  expect_that(getSegmentData(data)[[2]], equals(NULL))

  data <- data.frame(rmvnorm(100, rep(0, 7)))
  data$group <- rep(LETTERS[1:2], each=50)
  data[, 1:7] <- t(apply(data[, 1:7], 1, sort))
  
  expect_that(getSegmentData(data)[[1]], is_equivalent_to(data[, c("X1", "X5", "group")]))
  expect_that(getSegmentData(data)[[2]], is_equivalent_to(data[, c("X2", "X4", "group")]))

  data <- data.frame(rmvnorm(100, rep(0, 2)))
  data$group <- rep(LETTERS[1:2], each=50)
  data[, 1:2] <- t(apply(data[, 1:2], 1, sort))
  
  expect_that(getSegmentData(data), throws_error())
})

test_that("simBase returns a sample as expected", {
  ll <- liver
  ll$ndose <- as.numeric(ll$dose)
  mm <- lmr(log(ALT.M) ~ log(ALT.B) + ndose, data=ll)
  ll$r <- resid(mm)
  nsim <- 10500; burn <- 500
  gm <- evm(r, data=ll, qu=.5, xi=~ndose, method="sim",
            iter=nsim, burn=burn, thin=1, verbose=FALSE)
  b <- simBase(mm, gm, baseline="log(ALT.B)")

  # Check length is the same as the number of simulated
  # parameters in the GPD model
  expect_that(length(b), equals(nsim - burn))
  
  # Check the simualted baselines are from the same distribution
  # as the original. Don't want to be warned about ties in the data,
  # so suppress warnings.
  suppressWarnings(k <- ks.test(log(ll$ALT.B), b)$p.value)
  expect_that(k > .2, is_true())
})

test_that("Distribution of expected values from robust regression is ok", {
  ll <- liver
  ll$ndose <- as.numeric(ll$dose)
  mm <- lmr(log(ALT.M) ~ log(ALT.B) + ndose, data=ll)
  ll$r <- resid(mm)
  nsim <- 10500; burn <- 500
  gm <- evm(r, data=ll, qu=.5, xi=~ndose, method="sim",
            iter=nsim, burn=burn, thin=1, verbose=FALSE)
  nd <- data.frame(ndose=1:4)
  sl <- simLinear(mm, gm, nd,
                  baseline="log(ALT.B)", rawBaseline="ALT.B",
                  invtrans=exp)

  # Check number of simulations is same as number of GPD
  # simualted parameters times number of rows in newdata
  expect_that(nrow(sl), equals(nrow(nd) * (nsim - burn)))
  
  # Check the (nsim-burn) baseline values got replicated correctly
  expect_that(log(sl$ALT.B[1:(nsim-burn)]), equals(log(sl$ALT.B[((nsim-burn)+1):(2*(nsim-burn))])))
  expect_that(log(sl$ALT.B[1:(nsim-burn)]), equals(log(sl$ALT.B[(2*(nsim-burn)+1):(3*(nsim-burn))])))
  expect_that(log(sl$ALT.B[1:(nsim-burn)]), equals(log(sl$ALT.B[(3*(nsim-burn)+1):nrow(sl)])))
  
  # Check that coefficients match those from the original robust regression.
  # The data have not been centred, so omit the intercept.
  ols <- lm(fitted ~ log(ALT.B) + ndose, data=sl)
  # Data not centred, so intercepts could differ
  expect_that(coef(ols)[-1], equals(coef(mm)[-1], tol=.01))

  # Check that the difference in fitted values is the same from one dose to
  # the next (which checks the simulated regression coefficients have been
  # replicated in the same way as the baselines).
  ssl <- split(sl, sl$ndose)
  diff1 <- ssl[[2]]$fitted - ssl[[1]]$fitted
  diff2 <- ssl[[3]]$fitted - ssl[[2]]$fitted
  diff3 <- ssl[[4]]$fitted - ssl[[3]]$fitted
  expect_that(diff1, equals(diff2))
  expect_that(diff2, equals(diff3))
  
  # Check the differences correspond to the dose effect
  expect_that(mean(diff1), equals(as.numeric(coef(mm)[3]), tol=.01))
  co <- coef(summary(mm))
  # tol is mean relative diff tol. Simulation suggests 0.025 should be fairly safe
  expect_that(sd(diff1), equals(co[3,2], tol=.03))
})

test_that("Probabilities of threshold exceedance are calculated correctly", {
    # margarita.rp and margarita.getProbs are just wrappers that allow for
    # dealing with matrices and lists of parmeters and return levels.
    for (i in 1:10){
        sigma <- abs(rnorm(1, 10))
        xi <- rnorm(1, 0, .5)
        phi <- log(sigma)
        m <- 2
        mrp <- margarita.rp(1, matrix(m), u=0, phi, xi, 1, NULL)
        p <- 1-pgpd(2, u=0, sigma=sigma, xi=xi)
        mgp <- margarita.getProbs(1, list(matrix(c(phi, xi), ncol=2)), u=0, p=1, r=NULL, m=matrix(2))
        expect_that(mrp, equals(p))
        expect_that(unlist(mgp), equals(p))
    }
})

test_that("margaritaScale traps erroneous input and returns correct substrig", {

    expect_that(margaritaScale("robots"), matches("r"))
    expect_that(margaritaScale("pingpong"), matches("p"))
    expect_that(margaritaScale("dumbdumb"), matches("d"))
    expect_that(margaritaScale(2), throws_error())
    expect_that(margaritaScale(letters), throws_error())
    f <- function(x){
            res <- try(margaritaScale(x), silent=TRUE)
            if (class(res) == "try-error"){ FALSE }
            else { TRUE }
         }
    s <- sapply(letters, f)  
    expect_that(letters[s], matches(c("d", "p", "r")))
})
}
