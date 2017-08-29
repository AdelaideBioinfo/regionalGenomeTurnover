
rm(list = ls())

specRef = "mm10"
queRef = "hg19"

# only works for mouse
keepX = FALSE


load(paste("~/Desktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/",specRef,".synthBin.genome.variable.RData", sep = ""))

load(file = paste("~/Desktop/RTN_domains/R_objects/netsAnalysis/hotspots/",specRef,"repNoRep.RData", sep = ""))


# maybe just rewrite the code to get the norm data
# Also normalize the other sections

gapType = c("refIns", "refDel", "queIns", "queDel")


intSigRanges = NULL
for(i in 1:length(gapType)){
  
  if((gapType[i] == "refDel" & specRef == "hg19") | (gapType[i] == "queDel" & specRef == "mm10")){
    intSigRange <-  repSigRanges[[gapType[i]]]
  }else{
    intSigRange <- intersect(noRepSigRanges[[gapType[i]]], repSigRanges[[gapType[i]]])
  }
  
  intSigRanges <- c(intSigRanges, list(intSigRange))
}
names(intSigRanges) <- gapType






## Bin rates for vairables
df<- data.frame(synthBin.gr)
colChoice <- c("dnasePeaks","ladCov","exon", "intron","ctcfMotif", 
               "L1Motif", "prdm9Motif", "ancient", "new_L1", "new_SINE", 
               "old_L1", "recHotSpot","cpgCov")

remainRefBases <- width(synthBin.gr)[1] - (df$queIns + df$refDel + df$seqGap)
df[remainRefBases > 0,colChoice] <- (df[remainRefBases > 0,colChoice]/remainRefBases[remainRefBases > 0]) * width(synthBin.gr)[1]


remainingBases <- (df$end - df$start + 1) - (df$missingGap + df$seqGap)
df <- df[remainingBases > 150001,]
df <- df[complete.cases(df),]
synthBinNormVar.gr <- GRanges(df)
seqinfo(synthBinNormVar.gr) <- seqinfo(synthBin.gr)


chrAll <- seqlevels(synthBinNormVar.gr)
chrAll <- chrAll[-grep("_", chrAll)]
chrAll <- chrAll[!(chrAll == "chrM" | chrAll == "chrY")]

if(specRef == "mm10" & !keepX){
  chrAll <- chrAll[chrAll != "chrX"]
}

synthBinNormVar.gr <- synthBinNormVar.gr[seqnames(synthBinNormVar.gr) %in% chrAll]

dfMcol <- data.frame(mcols(synthBinNormVar.gr))
dfMcol <- dfMcol[,c(colChoice,"gcContent", "dnaseActivity", "recombRate")]

meanDiffCalc <- function(x, df){
  return(colMeans(df[x,]) - colMeans(df[!x,]))
}


zAll = NULL

# do this for each gap type

for(gType in gapType){
  ols <- overlapsAny(synthBinNormVar.gr,intSigRanges[[gType]], minoverlap = 1)
  
  print(sum(ols))
  print(sum(width(intSigRanges[[gType]]))/2e5)
  
  meanDiff <- colMeans(dfMcol[ols,]) - colMeans(dfMcol[!ols,])
  
  a <- replicate(n = 10000, expr = meanDiffCalc(x = sample(ols, size = length(ols), replace = FALSE) , df = dfMcol))
  a <- data.frame(t(a))
  
  
  varType = "intron"
  hist(a[,varType], breaks = 50)
  
  abline(v = meanDiff[varType])
  
  colSd <- apply(X = a, MARGIN = 2, FUN = sd)
  colMean <- colMeans(a)
  
  z <- (meanDiff - colMean) / colSd
  zAll <- rbind(zAll,z)
  print(z)
}

rownames(zAll) <- gapType

zAll

if(specRef == "hg19"){
  save(zAll, file = paste("~/Desktop/RTN_domains/R_objects/netsAnalysis/hotspots/",specRef,".permZ.RData", sep = ""))
}

if(specRef == "mm10"){
  if(keepX){
    save(zAll, file = paste("~/Desktop/RTN_domains/R_objects/netsAnalysis/hotspots/",specRef,".permZkeepX.RData", sep = ""))
  }else{
    save(zAll, file = paste("~/Desktop/RTN_domains/R_objects/netsAnalysis/hotspots/",specRef,".permZnoX.RData", sep = ""))
  }
}


# so now we can just calculate the Z score for enrichemnt levels
# anything within 3 sd we can remove
# or we can change the size of the square based on the P value

