myreduce <- function(x, by=character(), drop.empty.ranges=FALSE,
                   min.gapwidth=1L, with.inframe.attrib=FALSE)
          {
            if (!isTRUEorFALSE(drop.empty.ranges))
              stop("'drop.empty.ranges' must be TRUE or FALSE")
            if (!isSingleNumber(min.gapwidth))
              stop("'min.gapwidth' must be a single integer")
            if (!is.integer(min.gapwidth))
              min.gapwidth <- as.integer(min.gapwidth)
            if (min.gapwidth < 0L)
              stop("'min.gapwidth' must be non-negative")
            FUN <- function(y) {
              name <- names(y)
              ranges <- ranges(y)[[1L]]
              values <- values(y)[[1L]]
              inds <-
                unname(split(seq_len(nrow(values)), lapply(values, as.vector)))
              rlist <-
                lapply(inds, function(i) {
                  rngs <-
                    reduce(ranges[i],
                           drop.empty.ranges=drop.empty.ranges,
                           min.gapwidth=min.gapwidth,
                           with.inframe.attrib=with.inframe.attrib)
                  list(ranges = rngs,
                       values =
                         values[rep(i, length.out = length(rngs)), ,
                                drop=FALSE])
                })
              ranges <-
                IRangesList(do.call(c, lapply(rlist, "[[", "ranges")))
              names(ranges) <- name
              values <-
                SplitDataFrameList(do.call(rbind,
                                           lapply(rlist, "[[", "values")))
              names(values) <- name
              new2(class(y), ranges = ranges, values = values, check = FALSE)
            }
            if (ncol(x) == 0 || length(by) == 0) {
              ranges <-
                reduce(ranges(x),
                       drop.empty.ranges = drop.empty.ranges,
                       min.gapwidth = min.gapwidth,
                       with.inframe.attrib = with.inframe.attrib)
              listData <- new2("DataFrame", nrows=sum(elementNROWS(ranges)),
                               check=FALSE)
              end <- cumsum(elementNROWS(ranges))
              names(end) <- names(ranges)
              partitioning <- PartitioningByEnd(end)
              initialize(x,
                         ranges = ranges,
                         values = relist(listData, partitioning))
            } else {
              endoapply(x[,by], FUN)
            }
          }
