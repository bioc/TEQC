get.reads <-
function(readsfile, filetype=c("bed", "bam"), chrcol=1, startcol=2, endcol=3, idcol, zerobased=TRUE, sep="\t", skip=1, header=FALSE, ...){

  filetype <- match.arg(filetype)
  
  if(length(grep(".bam$", readsfile) == 0) & filetype == "bed"){
    filetype <- "bam"
    message("'filetype' was changed to 'bam'")
  }

  if(filetype == "bam"){
    # read BAM file
# !! add flag to read only unique (primary) alignments 
    param <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=FALSE, isSecondaryAlignment=FALSE), what=c("qname", "pos", "qwidth", "rname"))
    aln <- scanBam(readsfile, param=param)[[1]]
# !!
    # create GRanges object
    #rd <- with(aln, RangedData(IRanges(pos, width=qwidth), ID=qname,  space=rname))
    gr <- with(aln, GRanges(seqnames=rname, ranges=IRanges(pos, width=qwidth, names=qname)))
  }

  else {
    # read only the required columns
    n <- max(count.fields(readsfile, sep=sep))
    colclasses <- rep("NULL", n)
    colclasses[chrcol] <- "character"
    colclasses[c(startcol, endcol)] <- "integer"
    if(!missing(idcol))
      colclasses[idcol] <- "character"

    dat <- read.delim(readsfile, colClasses=colclasses, sep=sep, skip=skip, header=header, ...)

    # sort reads (better for example when calling 'findOverlaps()' several times)
    #o <- order(dat[,chrcol], dat[,startcol])
    #if(!identical(o, 1:nrow(dat)))
    #  dat <- dat[o,]

    # make IRanges object
    if(missing(idcol))
      ir <- IRanges(start=dat[,startcol], end=dat[,endcol])
    else
      ir <- IRanges(start=dat[,startcol], end=dat[,endcol], names=dat[,idcol])
    
    # shift start position forward by 1 to go from 0-based to 1-based system
    if(zerobased)
      start(ir) <- start(ir) + 1

    # make GRanges object
    gr <- GRanges(seqnames=dat[,chrcol], ranges=ir)
  }

  # check for Illumina read pair IDs - #0/1 and #0/2 have to be removed
  #if(length(colnames(rd)) > 0){
    illu <- grep("#0/", names(gr))
    if(length(illu) > 0)
      names(gr) <- gsub("#0/[[:digit:]]", "", names(gr))
  #}
  
  # sort reads
  gr <- sortSeqlevels(gr)
  gr <- sort(gr)

  print(paste("read", length(gr), "sequenced reads"))
  return(gr)
}


