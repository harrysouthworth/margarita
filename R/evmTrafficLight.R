#' Plot predicted return levels of a safety lab variable with some guidance as to
#' whether it seems safe, potentially suspect, or dangerous
#' @param data An object created by \code{summary.simulate.margarita.prob}
#' @param ptcol The colour of points on the plot. Defaults to \code{ptcol="blue"}
#' @param loref The lower reference level beneath which a return level is considered
#'        to be acceptably low
#' @param hiref The upper reference level above which a return level is considered to
#'        be unacceptably high
#' @param refcol Background colours for the acceptable, intermediate, and unaccepable
#'        return levels
#' @param linecol Colour of the lines representing the credible intervals on the
#'        predicted return levels
#' @param ptsize The size of the point estimate (median) predicted return level
#' @param linesize The weight of the lines representing the credible intervals
#' @param scales Whether to allow trellis panel scales to match each other or be free
#' @param ncol The number of columns in the output plots
#' @param xlab The x-axis label
#' @param ylab The y-axis label
#' @param main The main title
#' @param ... Other arguments to the plotting functions. Currently unused
#' @details The function is a fairly simple modification of \code{plot.summary.simulate.margarita.prob}
#' @export
evmTrafficLight <- function(data=NULL, ptcol="blue",
                            loref=1/1000, hiref=1/200,
                            refcol=c("green", "orange", "red"),
                            linecol=c("blue", "blue"),
                            ptsize=4, linesize=c(.5, 1.5),
                            scales="free", ncol=NULL,
                            xlab="P( > RL)", ylab="", M, main=NULL,
                            ...){
  g <- names(data)
  data <- unclass(data)
  nM <- nrow(data[[1]])
  g <- rep(g, each=nM)

  # Add M to each data.frame
  if (missing(M))
    M <- factor(rownames(data[[1]]), levels=rownames(data[[1]]))

  data <- lapply(1:length(data), function(x, data, M) {
                                   data <- as.data.frame(data[[x]])
                                   data$M <- M
                                   data },
                 data=data, M=M)

  # Make groups to trellis on
  data <- do.call("rbind", data)
  data$groups <- factor(g, levels=unique(g))
  if (ncol(data) == 7){
    names(data)[3] <- "mid"
  } else if (ncol(data) == 5){
    names(data)[2] <- "mid"
  } else {
    stop("data object has wrong number of columns (should be 5 or 7)")
  }

  seg <- getSegmentData(data)
  seg <- lapply(seg, function(x, M){ if (!is.null(x)){ x$M <- M }; x }, M=data$M)

  p <-
  ggplot(data, aes(mid, groups)) +
    annotate("rect", ymin=-Inf, ymax=Inf, xmin=-Inf, xmax=loref, fill=refcol[1], alpha=0.5) +
    annotate("rect", ymin=-Inf, ymax=Inf, xmin=loref, xmax=hiref, fill=refcol[2], alpha=0.5) +
    annotate("rect", ymin=-Inf, ymax=Inf, xmin=hiref, xmax=Inf, fill=refcol[3], alpha=0.5) +
    geom_point(size=ptsize, color=ptcol) +
    facet_wrap(~M, scales=scales, ncol=ncol) +
    scale_x_continuous(xlab) +
    scale_y_discrete(ylab) +
    ggtitle(main) +
    geom_segment(data=seg[[1]], aes(x=lo, xend=hi, y=group, yend=group),
                 size=linesize[1], color=linecol[1]) +
    if (!is.null(seg[[2]])){
      geom_segment(data=seg[[2]], aes(x=lo, xend=hi, y=group, yend=group),
                   size=linesize[2], color=linecol[2])
    }
  p
}
