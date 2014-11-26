## ----system, echo=FALSE--------------------------------------------------
path <- "~/Work/repos/github/margarita/inst/auto/"
source(file.path(path, "system.R"))


## ----topmatter, echo=FALSE-----------------------------------------------
infile <- "userinput.txt"
#source(file.path(path, "topmatter.R"))
readAutoInputs(file.path(infile))


## ----data, echo=FALSE----------------------------------------------------
source(file.path(path, "data.R"))


## ----summary, echo=FALSE, results="asis"---------------------------------
print(xtable(fivenumBy(lab, baseline, arm), label="tab:sb",
             caption="Summary of baseline data by treatment group."))

agg <- if (minmax == "max") "maxima" else "minima"
print(xtable(fivenumBy(lab, value, arm), label="tab:sm",
             caption=paste("Summary of post-baseline", agg, "by treatment group.")))


## ----shiftplots, echo=FALSE, fig.cap=""----------------------------------
shiftplot(lab, aes(base, aval), by=~arm)


## ----robust, echo=FALSE, results="asis", fig.cap="Diagnostic plots for the robust linear model."----
rfo <- as.formula(paste0("tfun(", value, ") ~ ", "tfun(", baseline, ") + ", arm))
rmod <- lmr(rfo, data=lab)
rmod <- step.lmr(rmod, trace=FALSE)

wh <- coef(summary(rmod))
rownames(wh) <- tidyCoefNames(rmod)
rownames(wh) <- gsub("tfun", trans, rownames(wh))

print(xtable(wh, label="tab:robust",
             caption="Summary of the robust linear model used to eliminate effect of baseline and the effect of treatment on the location of the distribution."))

ggplot(rmod)
boxplot(rmod, by="arm")


## ----thresh, echo=FALSE, fig.cap="Threshold selection plots. The parameter estimates for the GP distribution should be constant above an appropriate threshold, the mean residual life plot should be  linear, and the confidence interval for the power parameter should include 1."----
lab$r <- resid(rmod)
ggplot(gpdThresh(lab$r))
th <- min(egp3Thresh(lab$r))
qu <- mean(lab$r <= th) # checked and it's <=

cqu <- ordinalIndicator(as.integer(round(qu*100)))
cstr <- paste0("$", as.integer(round(qu*100)), "^{", cqu, "}$")


## ----gpddiag, echo=FALSE, results="asis", fig.cap="Diagnostic plots for the selected GP model."----

# XXX The sane way of constructing the models list currently escapes me...
# XXX Need to try MLE, then try penalized if it fails. Not obvious if that can be
#     automated in general...

mods <- list()
mods[[1]] <- list()
mods[[1]][[1]] <- evm(r, data=lab, qu=qu) # null model

for (i in 1:length(models)){
  fo <- as.formula(models[i]) # do here to avoid scoping problems
  mods[[i+1]] <- vector("list", length=3)
  mods[[i+1]][[1]] <- evm(r, data=lab, qu=qu, phi=fo)
  mods[[i+1]][[2]] <- evm(r, data=lab, qu=qu, xi=fo)
  mods[[i+1]][[3]] <- evm(r, data=lab, qu=qu, phi=fo, xi=fo)
}

mods <- unlist(mods, recursive=FALSE)

print(AICtable(mods), include.rownames=TRUE)

prefmod <- which.min(sapply(mods, AIC))
gmod <- mods[[prefmod]]

if (prefmod == 1){
    prefstring <- paste0(prefmod, ", the null model with no terms for treatment.")
} else {
  prefstring <- paste0(prefmod, ".")
}

suppressMessages(ggplot(gmod))


## ----mcmcdiag, echo=FALSE, fig.cap="Plots of the Markov chains used to simulate from the posterior distribution of the parameters in the selected GP model."----
bmod <- evm(r, data=lab, qu=qu, phi=gmod$formulae[[1]], xi=gmod$formulae[[2]],
            method="sim", thin=thin, burn=burn, iter=iter, verbose=FALSE)
ggplot(bmod)


## ----returnLevels, echo=FALSE, fig.cap="Predicted return levels in terms of multiples of ULN. The heavy lines represent 50\\% interval estimates and the lighter lines represent 90\\% interval estimates.", results="asis"----
if (prefmod != 0){ # need to construct newdata to pass into margarita
  # Add any transformed variables to the data
  for (i in 1:2){
    tp <- transInFormula(gmod$formulae[[i]])
    if (length(tp) > 0){
      m <- model.matrix(gmod$formulae[[i]], lab)
      m <- m[, tp, drop=FALSE]
      m <- m[, ! colnames(m) %in% names(lab), drop=FALSE]
      lab <- cbind(lab, m)
    }
  }

  # Get names of all variables used in the models and reduce the data
  ap <- unique(varsInFormula(gmod$formula[[1]]),
               varsInFormula(gmod$formula[[2]]))
  nd <- unique(lab[, ap, drop=FALSE])

#  nd <- data.frame(sort(unique(lab[, arm])))
#  names(nd) <- arm
} else {
  nd <- NULL
}

marg <- margarita(rmod, bmod, newdata=nd, trans=tfun, invtrans=itfun,
                  baseline=baseline, minima=minmax=="min")

rl <- simulate(marg, M=returnLevels)
srl <- summary(rl)
srl <- srl/ULN
ggplot(srl, ncol=1, as.table=FALSE)

srl <- as.data.frame(srl)[, c(7:6, 1:5)]
srl[, 3:7] <- round(srl[, 3:7], 1)
ng <- length(unique(lab[, arm]))
names(srl)[2] <- "Group"
print(xtable(srl, label="tab:rl", digits=1,
             caption="Predicted M-subject return levels in terms of multiples of ULN."),
      include.rownames=FALSE, hline.after=c(-1, 0, 4*c(1:length(returnLevels))))


## ----multiplesOfULN, echo=FALSE, fig.cap="Predicted probability of any particular subject having a value exceed various multiples of the ULN.", results="asis"----
ep <- simulate(marg, type="prob", M=multiplesOfULN * ULN,
               Mlabels=paste0(multiplesOfULN, "xULN"))
sep <- summary(ep)
ggplot(sep, ncol=1, as.table=FALSE, scales="fixed")

sep <- as.data.frame(sep)
sep[, 3:ncol(sep)] <- sep[, 3:ncol(sep)] * 100
names(sep)[2] <- "Group"
print(xtable(sep, label="tab:ep", digits=2,
             caption="Predicted probability (\\%) of any particular subject having a value exceed various multiples of the ULN."),
      hline.after=c(-1, 0, (1:length(multiplesOfULN))*length(unique(lab[, arm]))),
      include.rownames=FALSE)


## ----sessionInfo, echo=FALSE, results="verbatim"-------------------------
sessionInfo()


## ----getenv, echo=FALSE, results="verbatim"------------------------------
as.matrix(Sys.getenv())


