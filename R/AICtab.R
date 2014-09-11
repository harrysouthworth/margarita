
#' Get an xtable object that describes change in deviance in terms of
#' strength of evidence
#' 
#' @import xtable
#' @export
evidence <- function(){
    dev <- cbind(c("< 2", "2 -- 6", "6 -- 10", "> 10"),
                 c("Weak", "Modest", "Strong", "Very strong"))
    colnames(dev) <- c("Change in deviance", "Strength of evidence")
    rownames(dev) <- c("", " ", "  ", "   ")
    xtable(dev, label="tab:deviance", align="llr",
           caption="Interpretation of change in deviance in terms of strength of evidence.")
}

#' Take a list of evmOpt objects and return an xtable
#' 
#' @param x A list, each element of which is an object of class "evm". It is
#'          assumed that the first element of \code{x} is the null model.
#' @param digits The number of significant digits to display.
#' @param label The label to use for the table. Defaults to \code{label='tab:aic'}.
#' @param caption The caption for the table.
#' @details Use \code{print(res, sanitize.text.function=function(x) x, include.rownames=FALSE)}.
#' @export
AICtable <- function(x, digits=3, label="tab:aic",
                     caption="Comparison of various GPD models in terms of number of coefficients, log-likelihood, AIC and change in deviance from the null model."){
    aic <- sapply(x, function(x){ AIC(x) })
    ll <- sapply(x, function(x){ x$loglik })
    dev <- as.character(signif(2*(ll - ll[1]), digits))
    dev[1] <- "-"
    np <- sapply(x, function(x) length(coef(x)))
    
    getfo <- function(wh){
        fo <- wh$formulae # <----------------------------------- XXX XXX XXX XXX
        fo <- lapply(fo, function(x) unlist(strsplit(as.character(x)[2], " + ", fixed=TRUE)))
        fo <- sapply(fo, function(x) paste0("f(", paste(x, collapse=", "), ")"))
        fo <- paste0("$", paste(paste0("\\", names(fo), "=", fo), collapse=", "), "$")
        fo
    }
    fo <- sapply(x, getfo)
    
    res <- cbind(fo, np, signif(ll, digits), signif(aic, digits), dev)
    colnames(res) <- c("Model", "\\#Coef.", "Loglik.", "AIC", "$\\Delta$ Dev.")
    # Even if print.xtable is told to ignore rownames, the align argument needs to assume
    # they're there, which is why it appears to have the wrong length.
    res <- xtable(res, label=label, caption=caption, align=c("l", "l", rep("r", 4)))
    oldClass(res) <- "AICtable"
    res
}

#' @method print AICtable
print.AICtable <- function(x, ...){
  print.xtable(x, include.rownames=FALSE, sanitize.colnames.function=function(x) x)
  invisible()
}