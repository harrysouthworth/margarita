context("margarita helpers")

test_that("margaritaScale traps erroneous input and returns correct substring", {
  margaritaScale <- margarita:::margaritaScale
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
  for (i in 1:sum(s)) # s is a logical vector
    expect_that(letters[s][i], matches(c("d", "p", "r")[i]))
})

test_that("getCIquantiles behaves as expected", {
  getCIquantiles <- margarita:::getCIquantiles
  expect_that(getCIquantiles(.05), equals(c(0.025, 0.50, 0.975)))
  expect_that(getCIquantiles(.1), equals(c(0.050, 0.50, 0.950)))
  expect_that(getCIquantiles(c(.1, .5)), equals(c(.050, .250, .50, .750, .950)))
  expect_that(getCIquantiles(c(.5, .1)), equals(sort(getCIquantiles(c(.5, .1)))))
  expect_that(getCIquantiles(0), throws_error())
  expect_that(getCIquantiles(1), throws_error())
  expect_that(getCIquantiles(c(.05, .1, .5)), throws_error())
})

test_that("Segments are returned correctly", {
  getSegmentData <- margarita:::getSegmentData
  data <- data.frame(rmvnorm(100, rep(0, 3)))

  expect_that(getSegmentData(data), throws_error())

  data$groups <- rep(LETTERS[1:2], each=50)

  expect_that(getSegmentsData(data), throws_error())
  
  data[, 1:3] <- t(apply(data[, 1:3], 1, sort))

  expect_that(getSegmentData(data)[[1]], is_equivalent_to(data[, c("X1", "X3", "groups")]))
  expect_that(getSegmentData(data)[[2]], equals(NULL))

  data <- data.frame(rmvnorm(100, rep(0, 7)))
  data$groups <- rep(LETTERS[1:2], each=50)
  data[, 1:7] <- t(apply(data[, 1:7], 1, sort))

  expect_that(getSegmentData(data)[[1]], is_equivalent_to(data[, c("X1", "X5", "groups")]))
  expect_that(getSegmentData(data)[[2]], is_equivalent_to(data[, c("X2", "X4", "groups")]))

  data <- data.frame(rmvnorm(100, rep(0, 2)))
  data$groups <- rep(LETTERS[1:2], each=50)
  data[, 1:2] <- t(apply(data[, 1:2], 1, sort))

  expect_that(getSegmentData(data), throws_error())
})

test_that("summary.margarita arithmetic works", {
  # Implicity tests as.data.frame.summary.margarita, too
  liver$ndose <- as.numeric(liver$dose)
  mm <- lmr(log(ALT.M) ~ log(ALT.B) + ndose, data=liver)
  liver$r <- resid(mm)
  em <- evm(r, data=liver, qu=.5, xi=~ndose, method="sim", iter=10500, verbose=FALSE)
  mar <- margarita(mm, em, newdata=data.frame(ndose=1:4), baseline="ALT.B")


  rl <- simulate(mar, M=c(100, 500, 1000))
  srl <- summary(rl)

  # division
  d <- srl / 25
  whd <- lapply(srl, function(x) x/25)
  whd <- do.call("rbind", whd)
  expect_that(as.data.frame(whd), is_equivalent_to(as.data.frame(d)[, 1:5]))

  # multiplication
  d <- srl * 25
  whd <- lapply(srl, function(x) x*25)
  whd <- do.call("rbind", whd)
  expect_that(as.data.frame(whd), is_equivalent_to(as.data.frame(d)[, 1:5]))

  # subtraction
  d <- srl - 25
  whd <- lapply(srl, function(x) x-25)
  whd <- do.call("rbind", whd)
  expect_that(as.data.frame(whd), is_equivalent_to(as.data.frame(d)[, 1:5]))
  
  # addition
  d <- srl + 25
  whd <- lapply(srl, function(x) x+25)
  whd <- do.call("rbind", whd)
  expect_that(as.data.frame(whd), is_equivalent_to(as.data.frame(d)[, 1:5]))
  })



