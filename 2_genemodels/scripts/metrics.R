args <- commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages({
  library("Biostrings")
  library(stringr)
})

pred <- readAAStringSet(args[1])
elegansprot <- readAAStringSet(args[2])
blast = read.table(args[3], header = FALSE, sep = "", dec = ".")

colnames(blast) <- c("query","subject","ident","length","mismatch","gapop","qstart","qend","sstart","send","eval","bitscore")

prednames <- names(pred)
predseq <- paste(pred)
preddf <- data.frame(prednames,predseq)
eprotnames <- names(elegansprot)
eprotseq <- paste(elegansprot)
eprotdf <- data.frame(eprotnames,eprotseq)

#print(preddf[1:5,2])
#print(eprotdf[1:5,])
#print(blast[1:5,])

hits80 <- blast[blast$ident > 80,]
orderedhits <- hits80[order(hits80$query, hits80$eval),]
filteredhits <- orderedhits[!duplicated(orderedhits$query),]
filteredhits <- filteredhits[order(filteredhits$query),]


predhits <- preddf[preddf$prednames %in% filteredhits$query,]
predhits <- predhits[order(predhits$prednames),]
prothits <- eprotdf[eprotdf$eprotnames %in% filteredhits$subject,]
predhits$len <- str_count(predhits$predseq)
prothits$len <- str_count(prothits$eprotseq)

predhits$subject <- filteredhits$subject
predhits$slen <- prothits$len[match(predhits$subject, prothits$eprotnames)]
predhits$delta <- (predhits$len / predhits$slen)
ones <- predhits[predhits$len == predhits$slen,]
fives <- predhits[abs((predhits$len - predhits$slen)/predhits$slen) < 0.05 ,]

print(nrow(predhits))
print(nrow(ones))
print(nrow(fives))

mindelta = min(predhits$delta)
maxdelta = max(predhits$delta)
totalcount <- paste("Total Hits:",nrow(predhits),sep=" ")
onescount <- paste("1:1 Hits:",nrow(ones),sep=" ")
fivescount <- paste("5% Off Hits:",nrow(fives),sep=" ")

filenameRatio <- paste("ratiopl",args[1], sep="-")
filenameReport <- paste("report",args[1], sep="-")

write.table(predhits$delta , sep=",",  col.names=FALSE, row.names = FALSE,file = filenameRatio)
write.table(totalcount, append=TRUE, col.names=FALSE, row.names = FALSE,file = filenameReport)
write.table(onescount, append=TRUE, col.names=FALSE, row.names = FALSE,file = filenameReport)
write.table(fivescount, append=TRUE, col.names=FALSE, row.names = FALSE,file = filenameReport)
