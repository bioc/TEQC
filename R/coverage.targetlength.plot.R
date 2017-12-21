coverage.targetlength.plot <-
function(targets, plotcolumn, linecol=2, xlab, ylab, lwd, pch, cex, ...){

  mctargets <- mcols(targets)
  if(ncol(mctargets) == 0)
    stop("'targets' does not include values to plot")

  if(is.character(plotcolumn) & !(plotcolumn %in% names(mctargets)))
    stop("selected 'plotcolumn' does not exist")

  targetlen <- width(targets)
  y <- mctargets[,plotcolumn]

  # set graphical parameters
  if(is.numeric(plotcolumn))
    plotcolumn <- names(mctargets)[plotcolumn]
  if(missing(xlab)) xlab <- "Target length (bp)"
  if(missing(ylab)) ylab <- plotcolumn
  if(missing(lwd)) lwd <- 3
  if(missing(cex)) cex <- 2
  if(missing(pch)) pch <- "."

  # scatter plot
  plot(targetlen, y, xlab=xlab, ylab=ylab, pch=pch, cex=cex, ...)

  # loess curve
  if(length(unique(targetlen)) > 3)  # otherwise smooth.spline doesn't work
    lines(smooth.spline(x=targetlen, y=y), col=linecol, lwd=lwd)
}

