#' Get an xtable object that describes change in deviance in terms of
#' strength of evidence
#' 
#' @export
evidence <- function(){
    dev <- cbind(c("< 2", "2 -- 6", "6 -- 10", "> 10"),
                 c("Weak", "Modest", "Strong", "Very strong"))
    colnames(dev) <- c("Change in deviance", "Strength of evidence")
    rownames(dev) <- c("", " ", "  ", "   ")
    xtable(dev, label="tab:deviance",
           caption="Interpretation of change in deviance in terms of strength of evidence.")
}

#' Take a list of evmOpt objects and return an xtable
#' 
#' @param x A list, each element of which is an object of class "evm".
#' 
#' @export
AICtable <- function(x){
    aic <- sapply(x, function(x){ AIC(x) })
    ll <- sapply(x, function(x){ x$loglik })

    getfo <- function(wh){
        fo <- wh$formulae # <----------------------------------- XXX XXX XXX XXX
        fo <- lapply(fo, function(x) unlist(strsplit(as.character(x)[2], " + ", fixed=TRUE)))
        fo <- sapply(fo, function(x) paste0("f(", paste(x, collapse=", "), ")"))
        fo <- paste(paste0("\\", names(fo), "=", fo), collapse=", ")
        fo
    }
    fo <- sapply(x, getfo)

    res <- cbind(fo, format(ll, digits=4), format(aic, digits=4))
    colnames(res) <- c("Model", "Log-likelihood", "AIC")
    res
}