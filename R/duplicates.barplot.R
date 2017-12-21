duplicates.barplot <-
function(reads, targets, returnDups=FALSE, truncateX, col=c("red","lightblue"), xlab, ylab, ylim, ...){

  # in case 'reads' is output of 'reads2pairs' and contains also 'singleReads'
  if(is.list(reads) & ("readpairs" %in% names(reads)))
    reads <- reads$readpairs

  # which reads are on target
  on.target <- reads %over% targets
  reads.on <- reads[on.target,]
  reads.off <- reads[!on.target,]

  reads.onu <- unique(reads.on)
  reads.offu <- unique(reads.off)
  
  # multiplicities
  m.on <- countOverlaps(reads.onu, reads.on, type="equal")
  m.off <- countOverlaps(reads.offu, reads.off, type="equal")
  T1 <- table(m.on)
  T2 <- table(m.off)
  
  # prepare data for barplot
  N <- union(names(T1), names(T2))
  l <- length(N)
  barmat0 <- matrix(0, nrow=2, ncol=l)
  colnames(barmat0) <- N
  barmat0[1, names(T1)] <- T1
  barmat0[2, names(T2)] <- T2

  # plot fractions instead of absolute numbers (separate for on- and off-target)
  barmat <- barmat0 / rowSums(barmat0)

  # graphical parameters
  if(missing(xlab)) xlab <- "Read multiplicity"
  if(missing(ylab)) ylab <- "Fraction of reads"
  if(missing(ylim)) ylim <- c(0, max(barmat) + 0.05)

  if(!missing(truncateX))
    l <- truncateX
  B <- barplot(barmat[,1:l], beside=TRUE, col=col, xlab=xlab, ylab=ylab, ylim=ylim, ...)
  legend("topright", c("on target", "off target"), fill=col)

  # add absolute numbers (in millions) on top of the bars
  nr <- round(as.vector(barmat0[,1:l]) / 1e6, 2)
  text(x=as.vector(B), y=as.vector(barmat[,1:l])+ 0.02, labels=nr, cex=.8)

  if(returnDups){
    rownames(barmat) <- rownames(barmat0) <- c("on target", "off target")
    return(list(absolute=barmat0, relative=barmat))
  }
}
