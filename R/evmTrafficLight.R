probplot <- function(data=NULL, ptcol="blue",
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
  if (missing(M)){
    M <- factor(rownames(data[[1]]), levels=rownames(data[[1]]))
  }
  
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
    stop("data object has wrong number of columns (should be 5 or 7)")
  }
  
  seg <- margarita:::getSegmentData(data)
  seg <- lapply(seg, function(x, M){ if (!is.null(x)){ x$M <- M }; x }, M=data$M)
  
  p <- ggplot(data, aes(mid, groups)) +
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
