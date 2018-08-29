#' @method ggplot summary.margarita.sim.rl
#' @export
#' @importFrom scales comma
ggplot.summary.margarita.sim.rl <- function(data=NULL, trans="log10", labels=comma,
                                         xlab="Return level", ylab="", main=NULL,
                                         xbreaks = waiver(),
                                         ptcol="blue", linecol=c("blue", "blue"),
                                         ptsize=4, linesize=c(.5, 1.5),
                                         ncol=1, as.table=TRUE,
                                         ...){
    data <- as.data.frame(data)
    data$M <- factor(data$M, levels=unique(data$M))

    ng <- length(unique(data$groups))
    if (ng == 1) data$groups <- data$M # <------------------ Redundant now???

    nint <- ncol(data)/2 - .5 # Number of intervals

    names(data)[(ncol(data)-2)/2 + .5] <- "median" # Middle column (could be mean or median or something else)
#    data$group <- factor(rownames(data), levels=rownames(data))

    seg <- getSegmentData(data)
    seg[[1]]$M <- seg[[2]]$M <- data$M

    if (ng > 1){
      p <- ggplot(data=data, aes(median, groups)) +
             geom_point(size=ptsize, color=ptcol) +
             facet_wrap(~M, ncol=ncol, as.table=as.table) +
             scale_x_continuous(xlab, trans=trans, labels=labels, breaks=xbreaks) +
             scale_y_discrete(ylab) +
             ggtitle(main) +
             geom_segment(data=seg[[1]], aes(x=lo, xend=hi, y=group, yend=group),
                          size=linesize[1], color=linecol[1]) +
             if (!is.null(seg[[2]])){
               geom_segment(data=seg[[2]], aes(x=lo, xend=hi, y=group, yend=group),
                            size=linesize[2], color=linecol[2])
             } # Close if
    } # Close if ng > 1
    else{
      p <- ggplot(data=data, aes(median, groups)) +
             geom_point(size=ptsize, color=ptcol) +
             scale_x_continuous(xlab, trans=trans, labels=labels, breaks=xbreaks) +
             scale_y_discrete(ylab) +
             ggtitle(main) +
             geom_segment(data=seg[[1]], aes(x=lo, xend=hi, y=group, yend=group),
                          size=linesize[1], color=linecol[1]) +
             if (!is.null(seg[[2]])){
               geom_segment(data=seg[[2]], aes(x=lo, xend=hi, y=group, yend=group),
                            size=linesize[2], color=linecol[2])
             }
    }
    p
}

#' @method ggplot summary.margarita.sim.prob
#' @export
ggplot.summary.margarita.sim.prob <- function(data=NULL, ptcol="blue",
                                              linecol=c("blue", "blue"),
                                              ptsize=4, linesize=c(.5, 1.5),
                                              scales="free", ncol=NULL, as.table=TRUE,
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
                             data }, data=data, M=M)

    # Make groups to trellis on
    data <- do.call("rbind", data)
    data$groups <- factor(g, levels=unique(g))

    if (ncol(data) == 7){
        names(data)[3] <- "mid"
    }
    else if (ncol(data) == 5){
        names(data)[2] <- "mid"
    }
    else {
        stop("data object has wrong number of columns")
    }

    seg <- getSegmentData(data)
    seg <- lapply(seg, function(x, M){ if (!is.null(x)){ x$M <- M }; x }, M=data$M)

    p <- ggplot(data, aes(mid, groups)) +
             geom_point(size=ptsize, color=ptcol) +
             facet_wrap(~M, scales=scales, ncol=ncol, as.table=as.table) +
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
