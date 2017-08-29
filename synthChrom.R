
library(dplyr)
library(GenomicRanges)

rm(list = ls())


devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/comparativeGenomics/netScripts/netDataFunctions.R")




specRef = "mm10"
specQue = "hg19"

load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/formattedNetData/",specRef,".",specQue,".netData.RData",sep = ""))
load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/stretchedGapAnnotation/",specRef,".stretch.RData", sep = ""))
load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/shiftData/",specRef,".expand.breaks.RData", sep = ""))


# annotate extra files

refMissingGaps.gr <- GenomicRanges::setdiff(refFillGaps.gr, refGap.gr, ignore.strand = TRUE)
all(!overlapsAny(refMissingGaps.gr, refGap.gr))

mcols(refMissingGaps.gr)$queRanges <- GRanges(seqnames = "empty", ranges = IRanges(start = 1, end = 1))
mcols(refMissingGaps.gr)$sData = "*"
mcols(refMissingGaps.gr)$chainID = NA
mcols(refMissingGaps.gr)$elementID = 1:length(refMissingGaps.gr)
mcols(refMissingGaps.gr)$type = "missingGap"
mcols(refMissingGaps.gr)$featureID = 1:length(refMissingGaps.gr)
mcols(refMissingGaps.gr)$refRanges = granges(refMissingGaps.gr, use.mcols = FALSE)


refSeqGaps.gr <- reduce(refSeqGaps.gr)
mcols(refSeqGaps.gr)$queRanges <- GRanges(seqnames = "empty", ranges = IRanges(start = 1, end = 1))
mcols(refSeqGaps.gr)$sData = "*"
mcols(refSeqGaps.gr)$chainID = NA
mcols(refSeqGaps.gr)$elementID = 1:length(refSeqGaps.gr)
mcols(refSeqGaps.gr)$type = "seqGap"
mcols(refSeqGaps.gr)$featureID = 1:length(refSeqGaps.gr)
mcols(refSeqGaps.gr)$refRanges = granges(refSeqGaps.gr, use.mcols = FALSE)


missingGenome.gr <- sort(c(refMissingGaps.gr, refSeqGaps.gr))
missingGenome.gr <- genoExpandBreak(missingGenome.gr, newSynthRefShift, seqlengths(stretchedRef.gr))



refSynth <- sort(c(stretchedRef.gr,missingGenome.gr))


# check for overlappign ranges
cov <- coverage(refSynth)
sl <- IRanges::slice(cov, lower = 2)
all(unlist(lapply(sl, length)) == 0)


# create synthetic genome
synthGenome <- GRanges(seqnames = seqlevels(stretchedRef.gr), 
                       ranges = IRanges(width = seqlengths(stretchedRef.gr), end = seqlengths(stretchedRef.gr)))
seqlengths(synthGenome) <- seqlengths(stretchedRef.gr)
genome(synthGenome) <- genome(stretchedRef.gr)


# select bin size
binSize <- 2e5

# bin synthetic genome
synthBin.gr <- unlist(slidingWindows(synthGenome, width = binSize, step = binSize))

# get overlapping ranges
ol <- findOverlaps(refSynth, synthBin.gr)
pInt <- pintersect(refSynth[queryHits(ol)], synthBin.gr[subjectHits(ol)])


# sort into data frame
df <- data.frame(width = width(pInt), type = pInt$type, binNo = subjectHits(ol))
sumDF <- summarise(group_by(df, type, binNo), total = sum(width))

# setup bin columns
mcols(synthBin.gr)$refIns <- 0
mcols(synthBin.gr)$refDel <- 0
mcols(synthBin.gr)$queIns <- 0
mcols(synthBin.gr)$queDel <- 0
mcols(synthBin.gr)$missingGap <- 0
mcols(synthBin.gr)$seqGap <- 0

# sort into bin columns
mcols(synthBin.gr[sumDF$binNo[sumDF$type == "refIns"]])$refIns = sumDF$total[sumDF$type == "refIns"]
mcols(synthBin.gr[sumDF$binNo[sumDF$type == "refDel"]])$refDel = sumDF$total[sumDF$type == "refDel"]
mcols(synthBin.gr[sumDF$binNo[sumDF$type == "queIns"]])$queIns = sumDF$total[sumDF$type == "queIns"]
mcols(synthBin.gr[sumDF$binNo[sumDF$type == "queDel"]])$queDel = sumDF$total[sumDF$type == "queDel"]
mcols(synthBin.gr[sumDF$binNo[sumDF$type == "missingGap"]])$missingGap = sumDF$total[sumDF$type == "missingGap"]
mcols(synthBin.gr[sumDF$binNo[sumDF$type == "seqGap"]])$seqGap = sumDF$total[sumDF$type == "seqGap"]



# at this point we get our syntheticBinned genome with the data
save(synthBin.gr, file = paste("~/Desktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/", genomes["ref"], ".synthBin.RData", sep = ""))






