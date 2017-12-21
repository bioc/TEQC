
reads2pairs <-
#function(reads){
function(reads, max.distance){

  if(is.null(names(reads)))
    stop("read pair IDs are needed as names of the 'reads' GRanges table")

  # check if there are chromosomes without any mapped reads and remove those
  #nreadsperchrom <- table(seqnames(reads))
  #zerochr <- names(nreadsperchrom)[nreadsperchrom == 0]
  #reads <- reads[!(seqnames(reads) %in% zerochr)]
  
  # check if there are as many 1st as 2nd reads on each chromosome
  # -> single reads or reads where pairs matching different chromosomes will be removed
  ID <- names(reads)
  dups <- tapply(ID, seqnames(reads), duplicated) 
  singlereads <- any(sapply(dups, function(x) sum(x) != length(x)/2))    
  
  # give back single reads separately
  if(singlereads){
    ID2 <- ID[unlist(dups)]     
    
    # check if there are any pairs at all
    if(length(ID2) == 0)
      stop("there are no read pairs")

    # check if some IDs are not unique to one pair  
    if(any(duplicated(ID2)))                       
      stop("read pair IDs do not seem to be unique") 
    
    nondups <- Rle(!(ID %in% ID2))     
    print(paste("there were", sum(nondups), "reads found without matching second read, or whose second read matches to a different chromosome"))
    res <- list(singleReads=reads[which(nondups),,drop=TRUE])
    reads <- reads[which(!nondups),,drop=TRUE]
  }
  
  # merge reads of each pair
  r <- sort(reads, by = ~ seqnames + names + start + end)
  
  # split ranges table into one for 1st and one for 2nd read
  x.split <- split(r, rep(1:2, length.out=length(r)))  
  
  # merge the tables such that each line starts with start position of first read
  #   of the pair and ends with end position of second read
  ir <- IRanges(start=start(x.split[[1]]), end=end(x.split[[2]]))
  
  res2 <- GRanges(seqnames=seqnames(x.split[[1]]), ranges=ir)
  names(res2) <- names(x.split[[1]])
  
  # remove read pairs that are more than max.distance apart (within same chromosome)
  if(!missing(max.distance)){
    toofar <- which(width(res2) > max.distance)
    if(length(toofar) > 0){
      print(paste("there were", length(toofar), "read pairs found with distance between the two reads exceeding max.distance =", max.distance))
      id.toofar <- names(res2)[toofar]
      reads.toofar <- reads[names(reads) %in% id.toofar,,drop=T]
      res2 <- res2[-toofar,,drop=T]
      
      # add them to 'singleReads' output
      if(singlereads)
        res$singleReads <- c(res$singleReads, reads.toofar)
      else
        res <- list(singleReads=reads.toofar)
      res$singleReads <- sort(res$singleReads)
    }
  }

  # sort again according to position
  res2 <- sort(res2)

  if(singlereads | !missing(max.distance))
    res <- c(res, readpairs=res2)
  else
    res <- res2
  return(res)
}
