

library(dplyr)
library(topGO)
library(GO.db)
library(wordcloud)
library(GenomicRanges)

rm(list = ls())


specRefs = c("hg19", "mm10")
specQues = c("mm10","hg19")

gapType = "queIns"

gapNames = c("gain","loss","gain", "loss")
names(gapNames) <- c("refIns", "refDel", "queIns", "queDel")

for(spec in 1:2){
  
  specRef = specRefs[spec]
  specQue = specQues[spec]
  
  for(gapType in names(gapNames)){
    
    
    
    
    if(specRef == "hg19"){
      library(TxDb.Hsapiens.UCSC.hg19.knownGene)
      txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
      library(org.Hs.eg.db)
      org.db <- org.Hs.eg.db
    }
    if(specRef == "mm10"){
      library(TxDb.Mmusculus.UCSC.mm10.knownGene)
      txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
      library(org.Mm.eg.db)
      org.db <- org.Mm.eg.db
    }
    
    
    
    
    load(paste("~/Desktop/RTN_domains/R_objects/netsAnalysis/hotspots/", specRef, "repNoRep.RData", sep = ""))
    load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/", specRef, ".synthBin.RData",sep = ""))
    load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/shiftData/", specRef, ".expand.breaks.RData", sep = ""))
    load(paste("~/Desktop/RTN_domains/R_objects/netsAnalysis/syntheticBinnedGenome/", specRef, "synthBinNorm.RData", sep = ""))
    # both of these things work through entrez gene IDs
    
    # now we can do the overlaps and get the necesary gene IDs
    
    devtools::source_url("http://raw.githubusercontent.com/ReubenBuck/RTN_domains_scripts/master/comparativeGenomics/netScripts/netDataFunctions.R")
    
    
    
    
    
    geneFeatures <- sort(sortSeqlevels(GenomicFeatures::genes(txdb)))
    gene <- genoExpandBreak(geneFeatures, synthGenome = newSynthRefShift, seqlengths(synthBinNorm.gr))
    
    if((gapType == "refDel" & specRef == "hg19") | (gapType == "queDel" & specRef == "mm10")){
       sigRange <- repSigRanges[[gapType]]
      }else{
    sigRange <- GenomicRanges::intersect(repSigRanges[[gapType]], noRepSigRanges[[gapType]])
     }
    
    ol <- findOverlaps(gene, sigRange)
    pInt <- pintersect(gene[queryHits(ol)], sigRange[subjectHits(ol)])
    
    olWidth <- data.frame(pInt) %>% 
      group_by(gene_id) %>%
      summarise(width = sum(width)) 
    
    geneWidth <-  width(GenomicFeatures::genes(txdb))
    names(geneWidth) <- GenomicFeatures::genes(txdb)$gene_id
    
    
    
    olWidth <- olWidth[olWidth$width/geneWidth[olWidth$gene_id] == 1,]
    keys <- olWidth$gene_id
    geneSymbol <- OrganismDbi::select(org.db, keys=keys, columns = c("SYMBOL","GO"))
    # save the gene symbol info elsewhere
    write.table(geneSymbol,quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE,
                file = paste("~/Desktop/RTN_domains/data/comparativeGenomics/hotspotsGenes/", 
                             specRef,"_", gapType,"_","hotspotGenes.txt", sep = ""))
    
    
    
    
    ol <- findOverlaps(gene, reduce(synthBinNorm.gr), type = "within")
    
    
    allKeysNames <- gene$gene_id[queryHits(ol)]
    allKeys <- rep(0, length(allKeysNames))
    names(allKeys) <- allKeysNames
    allKeys[keys] = 1
    
    sigGene <- function(allGene){
      return(allGene == 1)
    }
    
    if(specRef == "hg19"){
      sampleGOdata <- new("topGOdata",
                          description = "Simple session", ontology = "BP",
                          allGenes = allKeys,
                          geneSel = sigGene,
                          nodeSize = 10,
                          mapping = "org.Hs.eg.db",
                          annotationFun = annFUN.org
      )
    }
    if(specRef == "mm10"){
      sampleGOdata <- new("topGOdata",
                          description = "Simple session", ontology = "BP",
                          allGenes = allKeys,
                          geneSel = sigGene,
                          nodeSize = 10,
                          mapping = "org.Mm.eg.db",
                          annotationFun = annFUN.org
      )
    }
    
    resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")
    resultElim <- runTest(sampleGOdata, algorithm = "elim", statistic = "fisher")
    resultWeight <- runTest(sampleGOdata, algorithm = "weight", statistic = "fisher")
    resultParentChild <- runTest(sampleGOdata, algorithm = "parentchild", statistic = "fisher")
    
    resMethod <- c("classicFisher", "elimFisher","weightFisher" ,"parentChildFisher")
    allRes <- NULL
    for(i in resMethod){
      res <- GenTable(sampleGOdata, classicFisher = resultFisher, elimFisher = resultElim, 
                      weightFisher = resultWeight, parentChildFisher = resultParentChild,
                      topNodes = length(score(resultFisher)),orderBy = i)
      allRes <- c(allRes, list(res))
    }
    names(allRes) <-  resMethod
    
    cFisher <- allRes$classicFisher
    save(cFisher, file = paste("~/Desktop/RTN_domains/R_objects/netsAnalysis/goLists/",specRef, gapType,"GoTermLists.RData",sep = "") )
    
    algo <- c("Elim", "Weight", "ParentChild")
    algoName <- c("elim", "weight", "parent child")
    names(algoName) <- algo
    
    
    algoName <- c("classic", algoName)
    names(algoName) <- resMethod
    
    printResAll <- NULL
    for(i in resMethod){
      printRes <- data.frame(Algorithm = c(algoName[i], rep("",9)),
                             allRes[[i]][1:10,c("GO.ID", "Term", "Significant", "Expected", i)])
      rownames(printRes) <- NULL
      colnames(printRes)[ncol(printRes)] <- "p-value"
      printResAll <- rbind(printResAll,printRes)
    }
    
    
    # now we can give these name 
    # labels and captions
    xTab <- xtable::xtable(printResAll, 
                           caption = paste("Top 10 biological process GO terms for genes located in", specRef, gapNames[gapType], "hotspots.",
                                           "P-values for each GO term were calculated using the fisher statistic combined with one of four separate algorithms that each take the GO hierarchy into account (described in methods)"), 
                           label = paste("tab:",specRef,gapType,"goTerm", sep = ""))
    termTab <- print(xTab,include.rownames = FALSE, hline.after = c(0,0,10,20,30))
    
    fileConn<-file(paste("~/Desktop/RTN_domains/RTN_domain_plots/netGainLoss/goTables/",specRef, gapType,"GoTermTab.tex",sep = ""))
    writeLines(termTab, fileConn)
    close(fileConn)
    
    
  
    
  }
  
}


 
 
  
 
 
 
 
 