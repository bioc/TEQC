get.targets <-
function(targetsfile, chrcol=1, startcol=2, endcol=3, zerobased=TRUE, sep="\t", skip=1, header=FALSE, ...){

  # read only the required columns
  n <- max(count.fields(targetsfile, sep=sep))
  colclasses <- rep("NULL", n)
  colclasses[chrcol] <- "character"
  colclasses[c(startcol, endcol)] <- "integer"
  
  dat <- read.delim(targetsfile, colClasses=colclasses, sep=sep, skip=skip, header=header, ...)

  # make IRanges object
  ir <- IRanges(start=dat[,startcol], end=dat[,endcol])

  # shift start position forward by 1 to go from 0-based to 1-based system
  if(zerobased)
    start(ir) <- start(ir) + 1

  # make GRanges object
  gr <- GRanges(seqnames=dat[,chrcol], ranges=ir)

  # print number of targets and, if applicable, number of non-overlapping targets
  n.targ <- length(gr)
  gr <- reduce(gr)
  n.targ.no <- length(gr)
  if(n.targ == n.targ.no)
    print(paste("read", n.targ, "(non-overlapping) target regions"))
  else
    print(paste("read", n.targ, "target regions in total, which are collapsed to", n.targ.no, "non-overlapping target regions"))

  # sort by chromosome and position
  gr <- sortSeqlevels(gr)
  gr <- sort(gr)
  
  # return non-overlapping targets
  return(gr)
}

