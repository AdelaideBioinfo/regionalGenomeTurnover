rmSeqGapsFromNetOutput <- function(netOutput, seqGaps){
  # removes seq gaps from a net output file by taking the setdiff
  sDiff <- GenomicRanges::setdiff(netOutput, seqGaps)
  if(!all(overlapsAny(sDiff,netOutput))){
    stop("some netOutput ranges have been completely removed, check if net output file is correct")
  }
  ol <- findOverlaps(sDiff, netOutput)
  # this intersection makes sure no regions have been unnecesarily joined
  sDiff <- pintersect(sDiff[queryHits(ol)], netOutput[subjectHits(ol)])
  ol <- findOverlaps(sDiff, netOutput)
  if(!all(sDiff == sDiff[queryHits(ol)])){
    stop("not a single layer net output and orger has not been maintained, cannot correctly assign mcols")
  }
  mcols(sDiff) <- mcols(netOutput[subjectHits(ol)])
  return(sDiff)
}

switchGenome <- function(queGenome.gr){
  
  cName <- colnames(mcols(queGenome.gr))[colnames(mcols(queGenome.gr)) != "queRanges"]
  
  switched <- GRanges(mcols(queGenome.gr)$queRanges, 
          queRanges = granges(queGenome.gr, use.mcols = FALSE)
          )
  mcols(switched)[,cName] <- mcols(queGenome.gr)[,cName]
  return(switched)
}



genoExpandBreak <- function(x.gr, synthGenome, expandedSeqlengths){
  seqlengths(x.gr) <- seqlengths(synthGenome)
  ol <- findOverlaps(x.gr, synthGenome)
  pInt <- pintersect(x.gr[queryHits(ol)], synthGenome[subjectHits(ol)], drop.nohit.ranges=TRUE)
  seqlengths(pInt) <- expandedSeqlengths
  expanded.gr <- shift(pInt, shift = synthGenome$shift[subjectHits(ol)])
  genome(expanded.gr) <- paste("stretched.", genome(synthGenome), sep = "")
  return(expanded.gr)
}

genoExpandStretch <- function(x.gr, synthGenome, expandedSeqlengths){
  x.gr <- sort(sortSeqlevels(x.gr))
  seqlengths(x.gr) <- seqlengths(synthGenome)
  olStart <- findOverlaps(x.gr, synthGenome, select = "first")
  olEnd <- findOverlaps(x.gr, synthGenome, select = "last")
  expanded.gr <- GRanges(seqnames = seqnames(x.gr), 
                         ranges = IRanges(start = start(x.gr) + synthGenome[olStart]$shift,
                                          end = end(x.gr) + synthGenome[olEnd]$shift))
  mcols(expanded.gr) <- mcols(x.gr)
  seqlengths(expanded.gr) <- expandedSeqlengths
  genome(expanded.gr) <- paste("stretched.", genome(synthGenome), sep = "")
  return(expanded.gr)
}
