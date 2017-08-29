#!/usr/bin/env Rscript

# stretching chromosomes

# could write this as a scripts



rm(list = ls())

library(optparse)
option_list = list(
  make_option(c("-i", "--input"), default=NA, type='character',
              help="input Rdata file containting output from refGapFillStats"),
  make_option(c("-o", "--output"), default=NA, type='character',
              help="stretched genome output, svaed as RData"),
  make_option(c("-s", "--shiftGenome"), default=NA, type='character',
              help="newSynthRefShift genome file, requried coordiantes for genome expansion,saved as RData"),
  make_option(c("-R", "--useRepeats"), default=FALSE, type='logical', action = "store_true",
              help="option to use repeats to assign gaps"),
  make_option(c("-r", "--refRep"), default=NULL, type='character',
              help="reference genome new repeats, svaed as RData"),
  make_option(c("-q", "--queRep"), default=NULL, type='character',
              help="query genome new repeats, svaed as RData")
)
opt = parse_args(OptionParser(option_list=option_list))

if(any(is.na(opt))){
  stop("missing options")
}



library(dplyr)
library(zoo)
library(GenomicRanges)
devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/comparativeGenomics/netScripts/netDataFunctions.R")


load(opt$input)

if(opt$useRepeats){
  print("using repeats!!!!!!!!!")
  load(opt$refRep)
  refRep.gr <- newRep.gr
  refAncDna.gr <- gaps(refRep.gr)
  refAncDna.gr <- refAncDna.gr[strand(refAncDna.gr) == "*"]
  
  
  load(opt$queRep)
  queRep.gr <- newRep.gr
  queAncDna.gr <- gaps(queRep.gr)
  queAncDna.gr <- queAncDna.gr[strand(queAncDna.gr) == "*"]
  
}



## get the loss bases for the query genome
ol <- findOverlaps(queGap.gr,queAncDna.gr)
queGapAnc.gr <- pintersect(queGap.gr[queryHits(ol)], queAncDna.gr[subjectHits(ol)], drop.nohit.ranges=TRUE)

# get the gain bases for the query
ol <- findOverlaps(queGap.gr,
                   GenomicRanges::setdiff(queGap.gr,queAncDna.gr, ignore.strand = TRUE))
queGapNonAnc.gr <- pintersect(queGap.gr[queryHits(ol)], 
                              GenomicRanges::setdiff(queGap.gr,queAncDna.gr,ignore.strand = TRUE)[subjectHits(ol)], 
                              drop.nohit.ranges=TRUE)

# check seperation of ancestral and non ancestral
if(!all(!overlapsAny(queGapAnc.gr, queGapNonAnc.gr))){
  stop("ancestral and non andcetral did not seperate correctly for query")
}


## get the loss bases for the reference genome
ol <- findOverlaps(refGap.gr,refAncDna.gr)
refGapAnc.gr <- pintersect(refGap.gr[queryHits(ol)], refAncDna.gr[subjectHits(ol)], drop.nohit.ranges=TRUE)

# get the gain bases for the ref
ol <- findOverlaps(refGap.gr,
                   GenomicRanges::setdiff(refGap.gr,refAncDna.gr, ignore.strand = TRUE))
refGapNonAnc.gr <- pintersect(refGap.gr[queryHits(ol)], 
                              GenomicRanges::setdiff(refGap.gr,refAncDna.gr, ignore.strand = TRUE)[subjectHits(ol)], 
                              drop.nohit.ranges=TRUE)

# check seperation of ancestral and non ancestral
if(!all(!overlapsAny(refGapAnc.gr, refGapNonAnc.gr))){
  stop("ancestral and non andcetral did not seperate correctly for reference")
}



# set up genomes
queDel.gr <- refGapAnc.gr
queDel.gr$type = "queDel"
queDel.gr <- sort(queDel.gr)
queDel.gr$featureID = 1:length(queDel.gr)

refIns.gr <- refGapNonAnc.gr
refIns.gr$type = "refIns"
refIns.gr <- sort(refIns.gr)
refIns.gr$featureID = 1:length(refIns.gr)

queIns.gr <- switchGenome(queGapNonAnc.gr)
queIns.gr$type = "queIns"
queIns.gr <- sort(queIns.gr)
queIns.gr$featureID = 1:length(queIns.gr)

refDel.gr <- switchGenome(queGapAnc.gr)
refDel.gr$type = "refDel"
refDel.gr <- sort(refDel.gr)
refDel.gr$featureID = 1:length(refDel.gr)

mapped.gr <- c(refIns.gr,queDel.gr)
mapped.gr <- sort(sortSeqlevels(mapped.gr))

remapped.gr <- c(refDel.gr, queIns.gr)
remapped.gr <- sort(sortSeqlevels(remapped.gr))
# set widths for remapped bases
end(remapped.gr) = start(remapped.gr) + width(remapped.gr$queRanges)  -1

# need to make room for memmory
rm(refDel.gr, queIns.gr, queDel.gr, refIns.gr)
rm(refGapNonAnc.gr, refGapAnc.gr, queGapNonAnc.gr, queGapAnc.gr)
rm(queFill.gr, queGap.gr, refFill.gr, refGap.gr, refFillGaps.gr, queFillGaps.gr)

#mapped.gr <- mapped.gr[seqnames(mapped.gr) == "chr1"]
#remapped.gr <- remapped.gr[seqnames(remapped.gr) == "chr1"]

mapped.gr$refRanges <- granges(mapped.gr, use.mcols = FALSE)
remapped.gr$refRanges <- granges(remapped.gr, use.mcols = FALSE)


# calculate the total amount of shift required for each gap
agg <- aggregate(x = width(remapped.gr$queRanges), by = list(as.character(granges(remapped.gr, use.mcols = FALSE))), FUN = sum)


synthRefShift <- GRanges(agg$Group.1, shift = agg$x)
rm(agg)
seqlevels(synthRefShift) <- seqlevels(remapped.gr)
genome(synthRefShift) <- genome(remapped.gr)
seqlengths(synthRefShift) <- seqlengths(remapped.gr)
synthRefShift <- sort(sortSeqlevels(synthRefShift))
synthRefShift <- resize(synthRefShift, width = 1,fix = "start")

newDNA <- NULL
newSynthRefShift <- NULL
for(chr in seqlevels(synthRefShift)){
  newSynthRefShift0 <- synthRefShift[seqnames(synthRefShift) == chr]
  if(length(newSynthRefShift0) == 0){
    newSynthRefShift0 <- GRanges(seqnames = chr, ranges = IRanges(start = 1, width = 1))
  }else{
    startRange <- GRanges(seqnames = chr, 
                          ranges = IRanges(start = 1, end = start(newSynthRefShift0[1]) - 1), shift = 0)
    newSynthRefShift0 <- c(startRange, newSynthRefShift0)
  }
  newSynthRefShift0$shift[1] <- 0
  if(length(newSynthRefShift0) > 1){
    end(newSynthRefShift0) <- c(start(newSynthRefShift0)[2:length(newSynthRefShift0)] - 1, 
                                seqlengths(newSynthRefShift0)[chr])
  }else{
    end(newSynthRefShift0) <- seqlengths(synthRefShift)[chr]
  }
  newSynthRefShift0$shift <- cumsum(newSynthRefShift0$shift)
  newDNA <- c(newDNA,newSynthRefShift0$shift[length(newSynthRefShift0)] )
  newSynthRefShift <- c(newSynthRefShift, newSynthRefShift0)
}
names(newDNA) <- seqlevels(synthRefShift)
newSynthRefShift <- unlist(GRangesList(newSynthRefShift))

save(newSynthRefShift, file = opt$shiftGenome)


# newSynthRefShift spans the entire reference genome containg intervals between gaps, 
# intervals with a specifc betweenGap regions are are shifted by that regions shift value
# This makes room for the gaped areas of the genome


# shift our mapped ref indels to make way for the query indels
ol <- findOverlaps(mapped.gr, newSynthRefShift)
newMapped.gr <- pintersect(mapped.gr[queryHits(ol)], newSynthRefShift[subjectHits(ol)],drop.nohit.ranges=TRUE)
seqlengths(newMapped.gr) <- seqlengths(newMapped.gr) + newDNA
newMapped.gr <- shift(newMapped.gr, shift = newSynthRefShift$shift[subjectHits(ol)])


# shift query regions into gap spaces made in reference
# increase seqlengths to allow for shift
seqlengths(newSynthRefShift) <- seqlengths(newSynthRefShift) + 1
seqlengths(remapped.gr) <- seqlengths(remapped.gr) + 1
ol <- findOverlaps(resize(remapped.gr, width = 1, fix = "start"), shift(newSynthRefShift, 1))
# decrease agian back to normal
seqlengths(newSynthRefShift) <- seqlengths(newSynthRefShift) - 1
seqlengths(remapped.gr) <- seqlengths(remapped.gr) - 1

newRemapped.gr <- remapped.gr
seqlengths(newRemapped.gr) <- seqlengths(newMapped.gr)
newRemapped.gr <- shift(newRemapped.gr[queryHits(ol)], newSynthRefShift$shift[subjectHits(ol)])

# make sure there is no overlap between query and ref regions

if(!all(!overlapsAny(newRemapped.gr, newMapped.gr))){
  stop("not sufficint room for remapped gaps in stretched genome")
}


# sort negative strand info
newRemapped.gr[newRemapped.gr$sData == "-"] <- sort(newRemapped.gr[newRemapped.gr$sData == "-"],
                                                    by = ~ seqnames(queRanges) + start(queRanges), decreasing = TRUE)


# identify ties
ol <- findOverlaps(newRemapped.gr)
ol <- ol[!(isSelfHit(ol) | isRedundantHit(ol))]

# shift ties so they become untied
untieShift <- abs(start(newRemapped.gr)[subjectHits(ol)] - end(newRemapped.gr)[queryHits(ol)]) + 1
aggShift <- dplyr::summarise(
  dplyr::group_by(
    dplyr::data_frame(range = subjectHits(ol), untieShift = untieShift),
    range),
  untieShift = sum(untieShift)
)
newRemapped.gr[aggShift$range] <- shift(newRemapped.gr[aggShift$range], aggShift$untieShift)

# make sure there is no overlaps
if(!all(!overlapsAny(newRemapped.gr, newMapped.gr))){
  stop("untied remmaped gaps overlap mapped gaps")
}

# save files 

stretchedRef.gr <- c(newRemapped.gr, newMapped.gr)
stretchedRef.gr <- sort(sortSeqlevels(stretchedRef.gr))
genome(stretchedRef.gr) <- paste("stretched.",genomes["ref"], sep = "")

save(stretchedRef.gr, file = opt$output)



# how many pices per gap

#res <- data.frame(mcols(stretchedRef.gr)) %>% group_by(elementID) %>% summarise(n())


# is there a way to keep the unstreched range
# yes, we add the ref range to the end 
#




