# collect statistics for our gap finding process and ancestral genome analysis



# need to sort out coordinate system tommorow



rm(list = ls())

options(stringsAsFactors = FALSE)

library(RMySQL)
library(GenomicRanges)


## read functions



devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/comparativeGenomics/netScripts/netDataFunctions.R")

## specify files

# genomes
genomes <- c(ref = "hg19",que = "mm10")

#gap file
gapFiles <- c(ref = paste("~/Desktop/RTN_domains/data/comparativeGenomics/netAlignment/processed/", 
                          genomes["ref"], ".", genomes["que"], ".net.gaps", sep = ""),
              que = paste("~/Desktop/RTN_domains/data/comparativeGenomics/netAlignment/processed/", 
                          genomes["que"], ".", genomes["ref"], ".net.gaps", sep = ""))


# fill file
fillFiles <- c(ref = paste("~/Desktop/RTN_domains/data/comparativeGenomics/netAlignment/processed/", 
                           genomes["ref"], ".", genomes["que"], ".net.fills", sep = ""),
               que = paste("~/Desktop/RTN_domains/data/comparativeGenomics/netAlignment/processed/", 
                           genomes["que"], ".", genomes["ref"], ".net.fills", sep = ""))

# ancestral DNA file
ancDNAfiles <- c(ref = paste("~/Desktop/RTN_domains/data/comparativeGenomics/netAlignment/ancestralGenome/", 
                             genomes["ref"], ".merge.ancestral.pass.bed", sep = ""),
                 que =  paste("~/Desktop/RTN_domains/data/comparativeGenomics/netAlignment/ancestralGenome/", 
                              genomes["que"], ".merge.ancestral.pass.bed", sep = ""))


# read in info
for(i in 1:length(genomes)){
  # get chr info
  mychannel <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", db = genomes[i])
  chrInfo <- dbGetQuery(mychannel, "SELECT * FROM chromInfo;")
  #chrInfo <- chrInfo[-(grep(pattern = "_", x = chrInfo$chrom)),]
  
  # get seq gaps, 0 based half open
  seqGaps <- dbGetQuery(mychannel, "SELECT * FROM gap;")
  #seqGaps <- seqGaps[-(grep(pattern = "_", x = seqGaps$chrom)),]
  seqGaps.gr <- GRanges(seqnames = Rle(seqGaps$chrom), 
                        ranges = IRanges(start = seqGaps$chromStart + 1, 
                                         end = seqGaps$chromEnd)
  )
  seqlevels(seqGaps.gr) <- chrInfo$chrom
  seqlengths(seqGaps.gr) <- chrInfo$size
  genome(seqGaps.gr) <- genomes[i]
  
  seqGaps.gr <- sort(sortSeqlevels(seqGaps.gr))
  
  mychannel <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", db = genomes[genomes != genomes[i]])
  altChrInfo <- dbGetQuery(mychannel, "SELECT * FROM chromInfo;")
  
  
  # get gap and fill data for ref and que
  
  for(j in 1:2){
    netOutputFile <- get(paste(c("gap", "fill")[j],"Files", sep = ""))
    
    netOutput <- read.table(file = netOutputFile[i], header = FALSE,
                            col.names = c("seqnames", "start", "end","que.seqnames","que.start", "que.end","strand","chainID"),
                            colClasses = c("character", "integer", "integer","character", "integer", "integer","character","integer"))
    # convert to GRanges, move from 0 based to 1 based
    if(j == 2){
      netOutput[netOutput$strand == "-", c("que.start", "que.end")] <- netOutput[netOutput$strand == "-", c("que.end", "que.start")]
    }
    netOutput.gr <- GRanges(seqnames = netOutput$seqnames,
                            ranges = IRanges(start = netOutput$start,
                                             end  = netOutput$end),
                            queRanges = GRanges(seqnames = netOutput$que.seqnames,
                                                ranges = IRanges(start = netOutput$que.start,
                                                                 end = netOutput$que.end)
                            ),
                            strand = "*",
                            sData = netOutput$strand,
                            chainID = netOutput$chainID
                            
    )
    seqlevels(netOutput.gr) <- chrInfo$chrom
    seqlengths(netOutput.gr) <- chrInfo$size
    genome(netOutput.gr) <- genomes[i]
    
    seqlevels(netOutput.gr$queRanges) <- altChrInfo$chrom
    seqlengths(netOutput.gr$queRanges) <- altChrInfo$size
    genome(netOutput.gr$queRanges) <- genomes[genomes != genomes[i]]
    
    netOutput.gr <- rmSeqGapsFromNetOutput(netOutput = netOutput.gr, seqGaps = seqGaps.gr)
    
    netOutput.gr$queRanges <- sortSeqlevels(netOutput.gr$queRanges)
    netOutput.gr <- sort(sortSeqlevels(netOutput.gr))
    
    netOutput.gr$elementID = 1:length(netOutput.gr)
    
    
    assign(x = paste(names(genomes)[i],c("Gap","Fill")[j],".gr", sep = ""), value = netOutput.gr)
    
  }
  
  
  # get ancestral
  ancDna <- read.table(file = ancDNAfiles[i], header = FALSE,
                       col.names = c("seqnames", "start", "end","coverage"),
                       colClasses = c("character", "integer", "integer","integer"))
  # convert to GRanges, move from 0 based to 1 based
  ancDna.gr <- GRanges(seqnames = ancDna$seqnames,
                       ranges = IRanges(start = ancDna$start,
                                        end = ancDna$end)
  )
  ancDna.gr <- reduce(ancDna.gr)
  seqlevels(ancDna.gr) <- chrInfo$chrom
  seqlengths(ancDna.gr) <- chrInfo$size
  genome(ancDna.gr) <- genomes[i]
  
  ancDna.gr <- sort(sortSeqlevels(ancDna.gr))
  
  assign(x = paste(names(genomes)[i],"AncDna.gr", sep = ""), value = ancDna.gr)
  assign(x = paste(names(genomes)[i],"ChrInfo", sep = ""), value = chrInfo)
  assign(x = paste(names(genomes)[i],"SeqGaps.gr", sep = ""), value = seqGaps.gr)
  
  rm(chrInfo, seqGaps.gr, seqGaps, netOutput, netOutput.gr, ancDna, ancDna.gr)
  
}

all_cons <- dbListConnections(MySQL())
for(con in all_cons) {
  dbDisconnect(con)
}

# have read in all the data of both query and species


refFillGaps.gr <- gaps(refFill.gr)
refFillGaps.gr <- refFillGaps.gr[strand(refFillGaps.gr) == "*"]
refFillGaps.gr <- rmSeqGapsFromNetOutput(netOutput = refFillGaps.gr, seqGaps = refSeqGaps.gr)


queFillGaps.gr <- gaps(queFill.gr)
queFillGaps.gr <- queFillGaps.gr[strand(queFillGaps.gr) == "*"]
queFillGaps.gr <- rmSeqGapsFromNetOutput(netOutput = queFillGaps.gr, seqGaps = queSeqGaps.gr)


# filter recipirical best hits.

ol1 <- findOverlaps(queFill.gr$queRanges, refFill.gr)
ol2 <- findOverlaps(queFill.gr, refFill.gr$queRanges)

queBest.gr <- unique(queFill.gr[queryHits(intersect(ol1,ol2))])
refBest.gr <- unique(refFill.gr[subjectHits(intersect(ol1,ol2))])

queWorst.gr <- unique(queFill.gr[-queryHits(intersect(ol1,ol2))])
refWorst.gr <- unique(refFill.gr[-subjectHits(intersect(ol1,ol2))])

queGapBest.gr <- queGap.gr[overlapsAny(queGap.gr, setdiff(gaps(queBest.gr),queSeqGaps.gr), type = "equal")]
queGapWorst.gr <- queGap.gr[!overlapsAny(queGap.gr, setdiff(gaps(queBest.gr),queSeqGaps.gr), type = "equal")]

refGapBest.gr <- refGap.gr[overlapsAny(refGap.gr, setdiff(gaps(refBest.gr),refSeqGaps.gr), type = "equal")]
refGapWorst.gr <- refGap.gr[!overlapsAny(refGap.gr, setdiff(gaps(refBest.gr),refSeqGaps.gr), type = "equal")]


queGap.gr <- queGapBest.gr
queNonRBHGap.gr <- queGapWorst.gr
queFill.gr <- sort(queBest.gr)
queNonRBHFill.gr <- queWorst.gr

refGap.gr <- refGapBest.gr
refNonRBHGap.gr <- refGapWorst.gr
refFill.gr <- sort(refBest.gr)
refNonRBHFill.gr <- refWorst.gr






# save RObjects

save("genomes","queChrInfo", "queAncDna.gr", "queFill.gr", "queFillGaps.gr", "queGap.gr",
     "queNonRBHGap.gr", "queNonRBHFill.gr","queSeqGaps.gr", 
     "refChrInfo" ,"refAncDna.gr", "refFill.gr", "refFillGaps.gr", "refGap.gr", 
     "refNonRBHGap.gr", "refNonRBHFill.gr", "refSeqGaps.gr", 
     file = paste("~/Desktop/RTN_domains/R_objects/netsAnalysis/formattedNetData/", genomes["ref"],"." ,genomes["que"], ".netData.RData", sep = ""))







