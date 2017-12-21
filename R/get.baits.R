get.baits <-
function(baitsfile, chrcol=1, startcol=2, endcol=3, seqcol=4, zerobased=TRUE, sep="\t", header=TRUE, ...){

  dat <- read.delim(baitsfile, sep=sep, as.is=TRUE, header=header, ...)

  # sort baits
  o <- order(dat[,chrcol], dat[,startcol], dat[,endcol])
  if(!identical(o, 1:nrow(dat)))
    dat <- dat[o,]

  ir <- IRanges(start=dat[,startcol], end=dat[,endcol])

  # shift start position forward by 1 to go from 0-based to 1-based system
  if(zerobased)
    start(ir) <- start(ir) + 1

  # make GRanges object
  gr <- GRanges(seqnames=dat[,chrcol], ranges=ir, sequence=dat[,seqcol])
  gr <- sortSeqlevels(gr)
  gr <- sort(gr)
  n.bait <- length(gr)
  print(paste("read", n.bait, "hybridization probes"))
  return(gr)
}

