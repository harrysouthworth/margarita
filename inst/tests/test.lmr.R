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

  # A bug in step.lmRob means we need to take logs and create numeric dose
  ll <- liver
  ll$logALT.M <- log(ll$ALT.M)
  ll$logALT.B <- log(ll$ALT.B)
  ll$ndose <- as.numeric(ll$dose)
  r1 <- lmRob(logALT.M ~ logALT.B + ndose + TBL.B + TBL.M, data=ll,
              control=lmRob.control(eff=.85, weight=c("Bisquare", "Bisquare")))
  d1 <- drop1(r1)
  
  r2 <- lmr(log(ALT.M) ~ log(ALT.B) + as.numeric(dose) + TBL.B + TBL.M, data=liver)
  d2 <- drop1(r2)
  
  expect_equal(d2[, 2], d1[, 2]/nrow(liver), tol=.001)
})

test_that("step.lmr works", {
  # A bug in step.lmRob means we need to take logs and create numeric dose
  ll <- liver
  ll$logALT.M <- log(ll$ALT.M)
  ll$logALT.B <- log(ll$ALT.B)
  ll$ndose <- as.numeric(ll$dose)
  x <- rmvnorm(nrow(ll), mean=rep(0, 10))
  ll <- cbind(ll, data.frame(x))
  r1 <- lmRob(logALT.M ~ logALT.B + ndose + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data=ll,
              control=lmRob.control(eff=.85, weight=c("Bisquare", "Bisquare")))
  s1 <- step.lmRob(r1, trace=FALSE)
  
  r2 <- lmr(logALT.M ~ logALT.B + ndose + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data=ll)
  s2 <- step.lmr(r2, trace=FALSE)
  
  expect_equal(coef(s2), coef(s1), tol=.001)

  for (i in 1:10){
    ll <- liver
    x <- rmvnorm(nrow(ll), mean=rep(0, 10))
    ll <- cbind(ll, data.frame(x))
    r1 <- lmr(log(ALT.M) ~ log(ALT.B) + as.numeric(dose) + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, data=ll)
    s1 <- step.lmr(r1, trace=FALSE)
    expect_less_than(length(coef(s1)), length(coef(r1)))
  }
})


test_that("AIC.lmr does what it ought", {
  # Define some functions using code found at 
  #https://www.google.co.uk/url?sa=t&rct=j&q=&esrc=s&source=web&cd=2&cad=rja&uact=8&ved=0CCoQFjAB&url=https%3A%2F%2Ffeb.kuleuven.be%2Fpublic%2Fndbaf45%2Fpapers%2FTharmaratnamClaeskens.pdf&ei=cBdmVJmjHMPmaLOigfAO&usg=AFQjCNFGLJHy7L5oaMjmWX_0NmrKGJEFFg&sig2=sXjQyR-Tyz9UHEDeCZD5sQ&bvm=bv.79400599,d.d2s
  Rho = function(x, cc){
    U = x/cc; U1 = 3 * U^2 - 3 * U^4 + U^6
    U1[abs(U) > 1] = 1; return(U1)
  }
  Psi = function(x, cc){
    U = x/cc; U1 = 6/cc * U * (1 - U^2)^2
    U1[abs(U) > 1] = 0; return(U1)
  }
  
  dPsi = function(x, cc){
    U = x/cc; U1 = (6/(cc^2)) *(1- 6* U^2+ 5* U^4)
    U1[abs(U) > 1] = 0; return(U1)
  }
  AIC.S <- function(y, X, beta.s, scale.s, cc){
    n <- length(y)
    
    U=matrix(ncol=1,nrow=n)
    UU=matrix(ncol=1,nrow=n)
    UUU=matrix(ncol=1,nrow=n)
    for(i in 1:n){
      U[i,]=dPsi((y[i,]-X[i,] %*% beta.s)/scale.s ,cc)
      UU[i,]=(Psi((y[i,]-X[i,] %*% beta.s)/scale.s ,cc))^2
    }
    J= (t(X) %*% diag(as.vector(U)) %*% X * (1/(scale.s^2)))/n
    inv.J<- solve(J)
    K= (t(X) %*% diag(as.vector(UU)) %*% X * (1/(scale.s^2)))/n
    AIC =2*n*(log(scale.s))+ 2* sum(diag(inv.J %*%(K)))
    return(AIC)
  } # close AIC.S
  
  r1 <- lmr(log(ALT.M) ~ log(ALT.B) + as.numeric(dose), data=liver)
  a1 <- AIC.lmr(r1)
  a2 <- AIC.S(y=as.matrix(log(liver$ALT.M), ncol=1), r1$x, coef(r1), r1$s, r1$c)
  expect_equal(a1, a2)
})

