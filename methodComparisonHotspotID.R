#### Now we can compate both datasets on a bin by bin basis 
#### See how the distributions agree with each other

library(igraph)
library(spdep)
library(GenomicRanges)

rm(list = ls())




ZtoP <- function(Zs){
  2*pnorm(-abs(Zs))
}

ZsigAdj <- function(Zs, cutoff){
  score = p.adjust(p = ZtoP(Zs), method = "BH")
  return(score <= cutoff)
}




specRef = "mm10"
specQue = "hg19"

load(paste("~/Desktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/",specRef,".synthBin.RData", sep = ""))
noRepSynthBin.gr <- synthBin.gr
noRepSynthBin.df <- mcols(noRepSynthBin.gr)
noRepSynthBin.df <- noRepSynthBin.df[noRepSynthBin.df$missingGap + noRepSynthBin.df$seqGap < 20000,]

load(paste("~/Desktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/",specRef,".synthBin.rep.RData", sep = ""))
repSynthBin.gr <- synthBin.gr
repSynthBin.df <- mcols(repSynthBin.gr)
repSynthBin.df <- repSynthBin.df[repSynthBin.df$missingGap + repSynthBin.df$seqGap < 20000,]



pdf(file = paste("~/Desktop/RTN_domains/RTN_domain_plots/netGainLoss/methodCompare/",specRef,"scatterCompare.pdf", sep = ""))
layout(matrix(1:4, nrow = 2))
par(mar = c(2,2,2,2), oma = c(5,4,4,2))
smoothScatter(rank(repSynthBin.df$refIns), rank(noRepSynthBin.df$refIns),
              xlab = "", ylab = "", nrpoints = 0)
mtext(paste(specRef ,"gain"), side = 3)

smoothScatter(rank(repSynthBin.df$refDel), rank(noRepSynthBin.df$refDel),
              xlab = "", ylab = "", nrpoints = 0)
mtext(paste(specRef ,"loss"), side = 3)

smoothScatter(rank(repSynthBin.df$queIns), rank(noRepSynthBin.df$queIns),
              xlab = "", ylab = "", nrpoints = 0)
mtext(paste(specQue ,"gain"), side = 3)

smoothScatter(rank(repSynthBin.df$queDel), rank(noRepSynthBin.df$queDel),
              xlab = "", ylab = "", nrpoints = 0)
mtext(paste(specQue ,"loss"), side = 3)

title(main = specRef, outer = TRUE)
mtext(text = "Recent repeats based method (rank)", side = 1, outer = TRUE, line = 1)
mtext(text = "Ancestral bases based method (rank)", side = 2, outer = TRUE, line = 1)

dev.off()
##### so how do we compare our enriched regions

# should probably remove gap dense regions



gapChoice <- c("refIns", "refDel", "queIns", "queDel")
mehtods = c("rep", "noRep")

for(m in mehtods){
  
  df <- data.frame(get(paste(m,"SynthBin.gr", sep = "")))
  df$remainingBases <- (df$end - df$start + 1) - (df$missingGap + df$seqGap)
  
  df <- df[df$remainingBases > 150001,]
  df <- df[complete.cases(df),]
  
  df[,c("refIns", "refDel", "queIns", "queDel")] <- (df[,c("refIns", "refDel", "queIns", "queDel")] * 200000 )/ df$remainingBases 
  
  synthBinNorm.gr <- GRanges(df)
  
  seqinfo(synthBinNorm.gr) <- seqinfo(synthBin.gr)
  
  
  
  
  ol<-findOverlaps(synthBinNorm.gr, maxgap = 3*width(synthBinNorm.gr)[1])
  ol <- ol[!(isRedundantHit(ol))]
  # remove NA hits
  ol <- ol[!is.na(df$refIns[queryHits(ol)])]
  ol <- ol[!is.na(df$refIns[subjectHits(ol)])]
  olMat <- data.frame(ol)
  G <- graph.data.frame(d = olMat,directed=FALSE)
  #weight <- olMat$subjectHits - olMat$queryHits
  #weight <- -(weight - (max(weight) + 1)) / max(weight)
  weight = rep(1, length(ol))
  weight[isSelfHit(ol)] <- 0
  E(G)$weight <- weight
  #A <- as_adjacency_matrix(G,type="both",names=FALSE,sparse=TRUE,edges = FALSE)
  A <- as_adjacency_matrix(G,type="both",names=FALSE,sparse=TRUE,edges = FALSE, attr = "weight")
  wMat <- mat2listw(A)
  wMat = nb2listw(include.self(wMat$neighbours))
  
  gapChoice <- c("refIns", "refDel", "queIns", "queDel")
  dfGscore <- data.frame(synthBinNorm.gr)
  
  for(i in 1:length(gapChoice)){
    score <- df[,gapChoice[i]]
    G <- localG(x = score,wMat)
    dfGscore[,gapChoice[i]] <- G
  }
  
  
  sigRanges = NULL
  for(i in 1:length(gapChoice)){
    z <- ZsigAdj(dfGscore[,gapChoice[i]],cutoff = .05) & dfGscore[,gapChoice[i]] > 0
    sigRanges <- c(sigRanges , list((synthBinNorm.gr[z])))
  }
  names(sigRanges) <- gapChoice
  assign(x = paste(m, "SigRanges", sep = ""), value = sigRanges)
  
}




sigOlMat <- data.frame(t(matrix(NA, nrow = 3, ncol = 4,
                                dimnames = list(c("repUniq","intersect","ancUniq"),
                                                c("refIns", "refDel", "queIns", "queDel")
                                )
)
)
)


for(gChoice in gapChoice){
  
  sigOlMat[gChoice, "intersect"] <- sum(overlapsAny(repSigRanges[[gChoice]], noRepSigRanges[[gChoice]]))
  
  sigOlMat[gChoice, "ancUniq"] <-  length(noRepSigRanges[[gChoice]]) - sigOlMat[gChoice, "intersect"]
  
  sigOlMat[gChoice, "repUniq"] <- length(repSigRanges[[gChoice]]) - sigOlMat[gChoice, "intersect"]
  
}

sigOlMat * 200000/1e6


# select hotspots based on varified sites.
# use things where we can measure the error in the estimates
# gain hotspots are identified by transposon enrichent. 
# loss hotspots are identified by ancestral elements.


save(repSigRanges, noRepSigRanges, 
     file = paste("~/Desktop/RTN_domains/R_objects/netsAnalysis/hotspots/", specRef, "repNoRep.RData", sep = ""))

save(synthBinNorm.gr, 
     file = paste("~/Desktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/", specRef, "synthBinNorm.RData", sep = ""))

