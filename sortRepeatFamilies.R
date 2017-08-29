#### alignments and repeats
library(data.table)
library(MASS)
library(reshape)
library(dplyr)
library(GenomicRanges)





ldaLine <- function(fit){
  gmean <- fit$prior%*%fit$means
  const <- drop(gmean%*%fit$scaling)
  
  slope <- -fit$scaling[1]/fit$scaling[2]
  intercept <- const/fit$scaling[2]
  
  res <- c(intercept, slope)
  names(res) <- c("intercept", "slope")
  return(res)
}


specRef = "hg19"
specQue = "mm10"

load(paste("Desktop/RTN_domains/R_objects/netsAnalysis/formattedNetData/",specRef,".",specQue,".netData.RData",sep = ""))


load(paste("~/Desktop/RTN_domains/R_objects/rmskTables/",genomes[1],"/",genomes[1],".RData",sep = ""))

refRep <- rep

refRepAll.gr <- GRanges(seqnames = refRep$genoChr, 
                        ranges = IRanges(start = refRep$genoStart + 1, end = refRep$genoEnd),
                        repName = refRep$repName,
                        repClass = refRep$repClass,
                        perDiv = refRep$perDiv)


# query genome

load(paste("~/Desktop/RTN_domains/R_objects/rmskTables/",genomes[2],"/",genomes[2],".RData",sep = ""))

queRep <- rep

queRepAll.gr <- GRanges(seqnames = queRep$genoChr, 
                        ranges = IRanges(start = queRep$genoStart + 1, end = queRep$genoEnd),
                        repName = queRep$repName,
                        repClass = queRep$repClass,
                        perDiv = queRep$perDiv)





queRepAll.gr <- queRepAll.gr[!(queRepAll.gr$repClass == "Simple_repeat" | 
                                 queRepAll.gr$repClass == "Satellite" |  
                                 queRepAll.gr$repClass == "Low_complexity" |
                                 queRepAll.gr$repClass == "Satellite/centr")]
queRepAll.gr <- queRepAll.gr[width(queRepAll.gr) < 20000]


queRepGap.gr <- queRepAll.gr[overlapsAny(queRepAll.gr, GenomicRanges::setdiff(queGap.gr, queAncDna.gr), type = "within")]
queRepFill.gr <- queRepAll.gr[overlapsAny(queRepAll.gr, GenomicRanges::union(queAncDna.gr,queFill.gr), type = "within")]


dfQueGap <- data.frame(queRepGap.gr)
dfQueGap <- dfQueGap[ complete.cases(dfQueGap), ]
dfQueGap$type = "gap"
dfQueFill <- data.frame(queRepFill.gr)
dfQueFill <- dfQueFill[ complete.cases(dfQueFill), ]
dfQueFill$type = "fill"


dfQueAll <- rbind(dfQueGap, dfQueFill)

dfRepQue <- dfQueAll %>% 
  group_by(repName, type) %>% 
  summarise(divMean = mean(perDiv), number = n(), width = sum(width), 
            bot25 = quantile(perDiv, probs = .25),top25 = quantile(perDiv, probs = .75),
            mid50 = quantile(perDiv, probs = .5))

dfRepQue <- melt(setDT(dfRepQue), id=1:2)
dfRepQue <- cast(dfRepQue, repName  ~ type + variable)
dfRepQue[is.na(dfRepQue)] <- 0
dfRepQue$totalNo <- dfRepQue$fill_number + dfRepQue$gap_number
dfRepQue$totalWidth <- dfRepQue$fill_width + dfRepQue$gap_width



refRepAll.gr <- refRepAll.gr[!(refRepAll.gr$repClass == "Simple_repeat" | 
                                 refRepAll.gr$repClass == "Satellite" |  
                                 refRepAll.gr$repClass == "Low_complexity"|
                                 refRepAll.gr$repClass == "Satellite/centr" )]
refRepAll.gr <- refRepAll.gr[width(refRepAll.gr) < 20000]

refRepGap.gr <- refRepAll.gr[overlapsAny(refRepAll.gr, GenomicRanges::setdiff(refGap.gr, refAncDna.gr), type = "within")]
refRepFill.gr <- refRepAll.gr[overlapsAny(refRepAll.gr, GenomicRanges::union(refAncDna.gr,refFill.gr), type = "within")]


dfRefGap <- data.frame(refRepGap.gr)
dfRefGap <- dfRefGap[ complete.cases(dfRefGap), ]
dfRefGap$type = "gap"
dfRefFill <- data.frame(refRepFill.gr)
dfRefFill <- dfRefFill[ complete.cases(dfRefFill), ]
dfRefFill$type = "fill"

dfRefAll <- rbind(dfRefGap, dfRefFill)


dfRepRef <- dfRefAll %>% 
  group_by(repName, type) %>% 
  summarise(divMean = mean(perDiv), number = n(), width = sum(width), 
            bot25 = quantile(perDiv, probs = .25),top25 = quantile(perDiv, probs = .75), 
            mid50 = quantile(perDiv, probs = .5))

dfRepRef <- melt(setDT(dfRepRef), id=1:2)
dfRepRef <- cast(dfRepRef, repName  ~ type + variable)
dfRepRef[is.na(dfRepRef)] <- 0
dfRepRef$totalNo <- dfRepRef$fill_number + dfRepRef$gap_number
dfRepRef$totalWidth <- dfRepRef$fill_width + dfRepRef$gap_width




dfQueAllLda <- dfQueAll
familyGapOL <- dfRepQue$gap_width / dfRepQue$totalWidth
names(familyGapOL) <- dfRepQue$repName
dfQueAllLda$gapOL <-  familyGapOL[dfQueAllLda$repName]

dfRefAllLda <- dfRefAll
familyGapOL <- dfRepRef$gap_width / dfRepRef$totalWidth
names(familyGapOL) <- dfRepRef$repName
dfRefAllLda$gapOL <-  familyGapOL[dfRefAllLda$repName]


#fit <- lda(data = rbind(dfQueAllLda, dfRefAllLda), type ~ perDiv + gapOL)

fitRef <- lda(data = dfRefAllLda, type ~ perDiv + gapOL)
fitQue <- lda(data = dfQueAllLda, type ~ perDiv + gapOL)








pdf(file = "~/Desktop/RTN_domains/RTN_domain_plots/netGainLoss/sortRepeats/sortRepeats.pdf", width = 12, height = 6)

# we want to see if there is a difference between the divergence levels for each repeat
layout(matrix(1:2, nrow=1))
par(mar=c(5,5,2,2), oma = c(1,1,1,1))
for(i in c("Ref","Que")){
repSum <- get(paste("dfRep",i, sep = ""))

plot(repSum$gap_mid50, repSum$gap_width/repSum$totalWidth, pch = 16, cex = .3,
     ylab = "Overlap with non-ancestral sequence (%)",
     xlab = "Divergence from consensus (%)", type = "n", xlim = c(0,40),
     main = genomes[tolower(i)], yaxt = "n")
axis(side = 2,at = seq(0,1,.2), labels = seq(0,1,.2) *100)

rect(xleft = repSum$gap_bot25,
     xright = repSum$gap_top25,
     ytop = (repSum$gap_width/repSum$totalWidth) + (repSum$gap_width/5e8),
     ybottom = (repSum$gap_width/repSum$totalWidth) - (repSum$gap_width/5e8),
     col = scales::alpha("black", .2), border = NA)

rect(xleft = repSum$fill_bot25,
     xright = repSum$fill_top25,
     ytop = ( (repSum$gap_width/repSum$totalWidth)) + (repSum$fill_width/5e8),
     ybottom = ( (repSum$gap_width/repSum$totalWidth)) - (repSum$fill_width/5e8),
     col = scales::alpha("red", .2), border = NA)

#res <- ldaLine(fit)
#abline( res[c("intercept", "slope")], lty = 2, col = 3)
res <- ldaLine(get(paste("fit",i, sep = "")))
abline( res[c("intercept", "slope")], lty = 2, lwd = 2)

legend("bottomleft", legend = c("Non-ancestral", "Ancestral"), 
       col = c(scales::alpha("black", .2),scales::alpha("red", .2)),
       bty = "n", pch = 15)


}

dev.off()


# 
dfRepRefBoth <- dfRefAll %>% 
  group_by(repName) %>% 
  summarise(perDiv = mean(perDiv), width = sum(width))

familyGapOL <- dfRepRef$gap_width / dfRepRef$totalWidth
names(familyGapOL) <- dfRepRef$repName
dfRepRefBoth$gapOL <-  familyGapOL[dfRepRefBoth$repName]

dfRepQueBoth <- dfQueAll %>% 
  group_by(repName) %>% 
  summarise(perDiv = mean(perDiv), width = sum(width))

familyGapOL <- dfRepQue$gap_width / dfRepQue$totalWidth
names(familyGapOL) <- dfRepQue$repName
dfRepQueBoth$gapOL <-  familyGapOL[dfRepQueBoth$repName]


# now we can classify families based on lda

# the total sequnce assigned to each group is pretty constant in the query
# maybe not so much in the ref

refClass <- data.frame(class = predict(fitRef, dfRepRefBoth)$class, width = dfRepRefBoth$width)
rownames(refClass) <- dfRepRefBoth$repName

queClass <- data.frame(class = predict(fitQue, dfRepQueBoth)$class, width = dfRepQueBoth$width)
rownames(queClass) <- dfRepQueBoth$repName



queSharedName <- queClass[intersect(rownames(queClass), rownames(refClass)), ]
refSharedName <- refClass[intersect(rownames(queClass), rownames(refClass)), ]

sum(queSharedName$class == refSharedName$class)/nrow(queSharedName)

sum(queSharedName$width[queSharedName$class == refSharedName$class])
sum(queSharedName$width[queSharedName$class != refSharedName$class])



# four groups 
# with different names
# is this important
# does it take too long to explain


queShareFam <- queClass[intersect(rownames(queClass), rownames(refClass)), ] %>%
  group_by(class) %>%
  summarise( width = sum(width), n = n())


queNewFam <- queClass[setdiff(rownames(queClass), rownames(refClass)), ] %>%
  group_by(class) %>%
  summarise( width = sum(width), n = n())

queNameTable <- data.frame(matrix(NA, nrow = 2, ncol = 2,
                       dimnames = list(c("nonSharedFamily","sharedFamily"), 
                                       c("nonAncestral", "ancestral"))))
queNameTable <- list(width = queNameTable, n = queNameTable)

queNameTable$width["nonSharedFamily", "nonAncestral"] <- queNewFam[2,2]
queNameTable$width["sharedFamily", "nonAncestral"] <- queShareFam[2,2]
queNameTable$width["nonSharedFamily", "ancestral"] <- queNewFam[1,2]
queNameTable$width["sharedFamily", "ancestral"] <- queShareFam[1,2]

queNameTable$n["nonSharedFamily", "nonAncestral"] <- queNewFam[2,3]
queNameTable$n["sharedFamily", "nonAncestral"] <- queShareFam[2,3]
queNameTable$n["nonSharedFamily", "ancestral"] <- queNewFam[1,3]
queNameTable$n["sharedFamily", "ancestral"] <- queShareFam[1,3]





refShareFam <- refClass[intersect(rownames(queClass), rownames(refClass)), ] %>%
  group_by(class) %>%
  summarise( width = sum(width), n = n())

refNewFam <- refClass[setdiff(rownames(refClass), rownames(queClass)), ] %>%
  group_by(class) %>%
  summarise( width = sum(width), n = n())

refNameTable <- data.frame(matrix(NA, nrow = 2, ncol = 2,
                                  dimnames = list(c("nonSharedFamily","sharedFamily"), 
                                                  c("nonAncestral", "ancestral"))))

refNameTable <- list(width = refNameTable, n = refNameTable)

refNameTable$width["nonSharedFamily", "nonAncestral"] <- refNewFam[2,2]
refNameTable$width["sharedFamily", "nonAncestral"] <- refShareFam[2,2]
refNameTable$width["nonSharedFamily", "ancestral"] <- refNewFam[1,2]
refNameTable$width["sharedFamily", "ancestral"] <- refShareFam[1,2]

refNameTable$n["nonSharedFamily", "nonAncestral"] <- refNewFam[2,3]
refNameTable$n["sharedFamily", "nonAncestral"] <- refShareFam[2,3]
refNameTable$n["nonSharedFamily", "ancestral"] <- refNewFam[1,3]
refNameTable$n["sharedFamily", "ancestral"] <- refShareFam[1,3]


refNameTable$width/1e6
queNameTable$width/1e6


refNameTable$n
queNameTable$n

# sequeces potentially derived from repeats active during divergence.





refShare <- refClass[intersect(rownames(queClass), rownames(refClass)), ]
refNonShare <- refClass[setdiff(rownames(refClass), rownames(queClass)), ]

queShare <- queClass[intersect(rownames(queClass), rownames(refClass)), ]
queNonShare <- queClass[setdiff(rownames(queClass), rownames(refClass)), ]


pdf(file = "~/Desktop/RTN_domains/RTN_domain_plots/netGainLoss/sortRepeats/repeatDivergence.pdf", width = 12, height = 6 )
layout(matrix(c(1,2),nrow =1))
par(mar=c(5,5,2,2), oma = c(1,1,1,1))

plot(density(dfRefAll$perDiv[
  dfRefAll$repName %in% rownames(refShare[refShare$class == "gap",])]),
  main = "hg19", ylim = c(0,.13), xlim = c(0,50),
  col = 2, xlab = "Divergence from consensus (%)", lwd = 2)
lines(density(dfRefAll$perDiv[
  dfRefAll$repName %in% rownames(refShare[refShare$class == "fill",])]), col = 4, lwd = 2)
lines(density(dfRefAll$perDiv[
  dfRefAll$repName %in% rownames(refNonShare[refNonShare$class == "gap",])]), col= 1, lwd = 2)
#lines(density(dfRefAll$perDiv[
#  dfRefAll$repName %in% rownames(refNonShare[refNonShare$class == "fill",])]), col= 4)
legend("topright", bty = "n",
       legend = c("recent",
                  "divergence",
                  "ancestor"), 
       title = "Peroid of activity",
       col = c(1,2,4), lty = 1, lwd = 2)




plot(density(dfQueAll$perDiv[
  dfQueAll$repName %in% rownames(queShare[queShare$class == "gap",])]),
  main = "mm10", ylim = c(0,.13), xlim = c(0,50),
  col = 2, xlab = "Divergence from consensus (%)", lwd = 2)
lines(density(dfQueAll$perDiv[
  dfQueAll$repName %in% rownames(queShare[queShare$class == "fill",])]), col = 4, lwd = 2)
lines(density(dfQueAll$perDiv[
  dfQueAll$repName %in% rownames(queNonShare[queNonShare$class == "gap",])]), col = 1, lwd = 2)
#lines(density(dfQueAll$perDiv[
#  dfQueAll$repName %in% rownames(queNonShare[queNonShare$class == "fill",])]), col= 4)

legend("topright", bty = "n",
       legend = c("recent",
                  "divergence",
                  "ancestor"), 
       title = "Peroid of activity",
       col = c(1,2,4), lty = 1, lwd = 2)

dev.off()



# method overlap

newRep.gr <- refRepAll.gr[refRepAll.gr$repName %in% rownames(refClass[refClass$class == "gap",])]
repIns <- GenomicRanges::intersect(refGap.gr, newRep.gr)
repDel <- GenomicRanges::setdiff(refGap.gr, newRep.gr)

ancIns <- GenomicRanges::setdiff(refGap.gr, refAncDna.gr)
ancDel <- GenomicRanges::intersect(refGap.gr, refAncDna.gr)


gapSum <- sum(width(refGap.gr))

conTable <- matrix(data = NA, nrow = 2, ncol = 2, 
                   dimnames = list(c("repIns", "repDel"),
                                   c("ancIns", "ancDel")))

conTable["repIns", "ancIns"] <- sum(width(GenomicRanges::intersect(repIns, ancIns)))
conTable["repDel", "ancIns"] <- sum(width(GenomicRanges::intersect(repDel, ancIns)))
conTable["repIns", "ancDel"] <- sum(width(GenomicRanges::intersect(repIns, ancDel)))
conTable["repDel", "ancDel"] <- sum(width(GenomicRanges::intersect(repDel, ancDel)))


refConTable <- conTable / 1e6




newRep.gr <- queRepAll.gr[queRepAll.gr$repName %in% rownames(queClass[queClass$class == "gap",])]
repIns <- GenomicRanges::intersect(queGap.gr, newRep.gr)
repDel <- GenomicRanges::setdiff(queGap.gr, newRep.gr)

ancIns <- GenomicRanges::setdiff(queGap.gr, queAncDna.gr)
ancDel <- GenomicRanges::intersect(queGap.gr, queAncDna.gr)


gapSum <- sum(width(queGap.gr))

conTable <- matrix(data = NA, nrow = 2, ncol = 2, 
                   dimnames = list(c("repIns", "repDel"),
                                   c("ancIns", "ancDel")))

conTable["repIns", "ancIns"] <- sum(width(GenomicRanges::intersect(repIns, ancIns)))
conTable["repDel", "ancIns"] <- sum(width(GenomicRanges::intersect(repDel, ancIns)))
conTable["repIns", "ancDel"] <- sum(width(GenomicRanges::intersect(repIns, ancDel)))
conTable["repDel", "ancDel"] <- sum(width(GenomicRanges::intersect(repDel, ancDel)))


queConTable <- conTable / 1e6


queConTable
refConTable


# maybe classify as SINE, LINE, LTR, DNA


# classify repeats here and save the opjects, new, old, div

refRepAll.gr$time <- as.character(NA)

refRepAll.gr[refRepAll.gr$repName %in% rownames(refShare[refShare$class == "gap",])]$time <- "share;new"
refRepAll.gr[refRepAll.gr$repName %in% rownames(refShare[refShare$class == "fill",])]$time <- "share;old"
refRepAll.gr[refRepAll.gr$repName %in% rownames(refNonShare[refNonShare$class == "gap",])]$time <- "uniq;new"
refRepAll.gr[refRepAll.gr$repName %in% rownames(refNonShare[refNonShare$class == "fill",])]$time <- 'uniq;old'

table(refRepAll.gr[is.na(refRepAll.gr$time)]$repName)

refRepAll.gr <- refRepAll.gr[!is.na(refRepAll.gr$time)]


seqlevels(refRepAll.gr) <- refChrInfo$chrom
seqlengths(refRepAll.gr) <- refChrInfo$size
refRepAll.gr <- sort(sortSeqlevels(refRepAll.gr))
genome(refRepAll.gr) <- genomes["ref"]



queRepAll.gr$time <- as.character(NA)

queRepAll.gr[queRepAll.gr$repName %in% rownames(queShare[queShare$class == "gap",])]$time <- "share;new"
queRepAll.gr[queRepAll.gr$repName %in% rownames(queShare[queShare$class == "fill",])]$time <- "share;old"
queRepAll.gr[queRepAll.gr$repName %in% rownames(queNonShare[queNonShare$class == "gap",])]$time <- "uniq;new"
queRepAll.gr[queRepAll.gr$repName %in% rownames(queNonShare[queNonShare$class == "fill",])]$time <- 'uniq;old'

table(queRepAll.gr[is.na(queRepAll.gr$time)]$repName)

queRepAll.gr <- queRepAll.gr[!is.na(queRepAll.gr$time)]


seqlevels(queRepAll.gr) <- queChrInfo$chrom
seqlengths(queRepAll.gr) <- queChrInfo$size
queRepAll.gr <- sort(sortSeqlevels(queRepAll.gr))
genome(queRepAll.gr) <- genomes["que"]





allRep.gr <- refRepAll.gr
save(allRep.gr, file = paste("~/Desktop/RTN_domains/R_objects/rmskTables/",genomes["ref"],
                      "/",genomes["ref"],"NewOldFamily.RData", sep = "")
)

newRep.gr <- refRepAll.gr[grep("new",refRepAll.gr$time)]
save(newRep.gr, file = paste("~/Desktop/RTN_domains/R_objects/rmskTables/",genomes["ref"],
                             "/",genomes["ref"],"New.RData", sep = "")
)


allRep.gr <- queRepAll.gr
save(allRep.gr, file = paste("~/Desktop/RTN_domains/R_objects/rmskTables/",genomes["que"],
                             "/",genomes["que"],"NewOldFamily.RData", sep = "")
)

newRep.gr <- queRepAll.gr[grep("new",queRepAll.gr$time)]
save(newRep.gr, file = paste("~/Desktop/RTN_domains/R_objects/rmskTables/",genomes["que"],
                             "/",genomes["que"],"New.RData", sep = "")
)

