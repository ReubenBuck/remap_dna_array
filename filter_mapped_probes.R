rm(list = ls())
options(stringsAsFactors = FALSE)

args = commandArgs(trailingOnly=TRUE)
# datProbe datA datB out



datProbe <- read.csv(args[1], skip = 7)
datProbe <- datProbe[datProbe$AlleleA_ProbeSeq != "",]
dubAddressSnp <- datProbe[datProbe$AlleleB_ProbeSeq != "", "Name"]

datA <- read.table(args[2], 
                   col.names = c("qaccver", "saccver", "pident", "length", 
                                 "mismatch", "gapopen", "qstart", "qend", 
                                 "sstart", "send", "evalue", "bitscore", "sstrand", "btop")
)
datA[datA$sstrand == "minus" ,c("sstart","send")] <- datA[datA$sstrand == "minus" ,c("send","sstart")]
datA$pos[datA$sstrand == "plus" &
           !(datA$qaccver %in% dubAddressSnp & datA$qend == 50) ] <- datA$send[datA$sstrand == "plus"&
                                                                                 !(datA$qaccver %in% dubAddressSnp & datA$qend == 50)] + 1

datA$pos[datA$sstrand == "plus" & is.na(datA$pos)] <- datA$send[datA$sstrand == "plus" & is.na(datA$pos)]

datA$pos[datA$sstrand == "minus" &
           !(datA$qaccver %in% dubAddressSnp & datA$qend == 50) ] <- datA$sstart[datA$sstrand == "minus"&
                                                                                   !(datA$qaccver %in% dubAddressSnp & datA$qend == 50)] - 1
datA$pos[datA$sstrand == "minus" & is.na(datA$pos)] <- datA$sstart[datA$sstrand == "minus" & is.na(datA$pos)]

datB <- read.table(args[3], 
                   col.names = c("qaccver", "saccver", "pident", "length", 
                                 "mismatch", "gapopen", "qstart", "qend", 
                                 "sstart", "send", "evalue", "bitscore", "sstrand", "btop")
)
datB[datB$sstrand == "minus" ,c("sstart","send")] <- datB[datB$sstrand == "minus" ,c("send","sstart")]
datB$pos[datB$sstrand == "plus" &
           !(datB$qaccver %in% dubAddressSnp & datB$qend == 50)] <- datB$send[datB$sstrand == "plus"&
                                                                                !(datB$qaccver %in% dubAddressSnp & datB$qend == 50)] + 1
datB$pos[datB$sstrand == "plus" & is.na(datB$pos)] <- datB$send[datB$sstrand == "plus" & is.na(datB$pos)]

datB$pos[datB$sstrand == "minus" &
           !(datB$qaccver %in% dubAddressSnp & datB$qend == 50)] <- datB$sstart[datB$sstrand == "minus"
                                                                                & !(datB$qaccver %in% dubAddressSnp & datB$qend == 50)] - 1
datB$pos[datB$sstrand == "minus" & is.na(datB$pos)] <- datB$sstart[datB$sstrand == "minus" & is.na(datB$pos)]

if(any(!(datA$qaccver %in% datProbe$Name))){
  stop("snp names in A do not match manifest snp names")
}

if(any(!(datB$qaccver %in% datProbe$Name))){
  stop("snp names in B do not match manifest snp names")
}

library(schoolmath)
errNames <- c("3unmapped", "5unmapped",
         "mismatch", "gap", "lowerQual", "equalQual",
         "onlyA", "mismatchA", "gapA", "lowerQualA", "equalQualA",
         "onlyB", "mismatchB", "gapB", "lowerQualB", "equalQualB",
         "diffAB", "A>B", "A<B", "noMatch")
err <- primes()[2:(length(errNames) + 1)]
names(err) <- errNames

# we just multiply these together to get the warning status

datA$address <- "S"
datA$address[datA$qaccver %in% dubAddressSnp] <- "A"
datB$address <- "B"
datAll <- rbind(datA, datB)
datAll <- datAll[order(datAll$qaccver),]
datAll$qual <- datAll$qend - datAll$qstart + 1 - datAll$mismatch - datAll$gapopen
#datAll$qual[datAll$qaccver %in% dubAddressSnp & datAll$qstart == 1 & datAll$btop == "49"] <- datAll$qual[datAll$qaccver %in% dubAddressSnp & datAll$qstart == 1 & datAll$btop == "49"] + 1
# remove a point if non overlapping A and B
datMerge <- merge(datAll[datAll$address == "A",], datAll[datAll$address == "B",], "qaccver")
datMerge <- datMerge[datMerge$saccver.x == datMerge$saccver.y & datMerge$pos.x == datMerge$pos.y,]
noOLAB <- (1:nrow(datAll))[!(apply(datAll[,c("qaccver","saccver","pos")], MARGIN = 1, FUN = paste, collapse = "|") %in% 
         apply(datMerge[,c("qaccver","saccver.x","pos.x")], MARGIN = 1, FUN = paste, collapse = "|")) &
         datAll$address != "S"]
datAll$qual[noOLAB] <- datAll$qual[noOLAB] - 1

datAll <- datAll[order(datAll$qual, decreasing = TRUE),]

# now we get to matching the first element

datFirst <- datAll[!duplicated(datAll$qaccver),]
datNext <- datAll[duplicated(datAll$qaccver),]

newMap <- data.frame(name = datProbe$Name, chr = NA, pos = NA, row.names = datProbe$Name)
newMap[datFirst$qaccver,c("chr","pos")] <- datFirst[,c("saccver","pos")]

statusList <- list()
length(statusList) <- nrow(newMap)
names(statusList) <- rownames(newMap)

# now we can run our tests to add the numbers to the array

# second matchs for single mapped stuff
merS <- merge(datFirst[datFirst$address == "S",], datNext[datNext$address == "S",], by = "qaccver")
merS <- merS[order(merS$qual.y, decreasing = TRUE),]
merS <- merS[!duplicated(merS$qaccver),]

opt <- merS[merS$qual.x == merS$qual.y,"qaccver"]
opt <- opt[opt %in% datFirst$qaccver]
statusList[opt] <- lapply(statusList[opt], c, err["equalQual"])
rm(opt)

opt <- merS[merS$qual.x > merS$qual.y,"qaccver"]
opt <- opt[opt %in% datFirst$qaccver]
statusList[opt] <- lapply(statusList[opt], c, err["lowerQual"])
rm(opt)

# count mismatches
# count indels last

# might need a dat all maybe
# seccond matched in A probes
merA <- merge(datAll[datAll$address == "A",], datAll[datAll$address == "A",], by = "qaccver")
merA[merA$qual.x < merA$qual.y, c(grep("\\.x", colnames(merA)), grep("\\.y", colnames(merA)))] <- merA[merA$qual.x < merA$qual.y, c(grep("\\.y", colnames(merA)), grep("\\.x", colnames(merA)))]
merA <- merA[order(merA$qual.x, decreasing = TRUE),]
merA <- merA[apply(merA[,grep("\\.x", colnames(merA))],MARGIN = 1,FUN = paste, collapse = "|") != 
               apply(merA[,grep("\\.y", colnames(merA))],MARGIN = 1,FUN = paste, collapse = "|"),]
merA <- merA[apply(merA[,c("qaccver","saccver.x","pos.x")], MARGIN = 1, FUN = paste, collapse = "|") %in% 
                     apply(datFirst[datFirst$address != "S",c("qaccver","saccver","pos")], MARGIN = 1, paste, collapse = "|"),]
merA <- merA[order(merA$qual.y, decreasing = TRUE),]
merA <- merA[!duplicated(merA$qaccver),]

opt <- merA[merA$qual.x == merA$qual.y,"qaccver"]
opt <- opt[opt %in% datFirst$qaccver]
statusList[opt] <- lapply(statusList[opt], c, err["equalQualA"])
rm(opt)


opt <- merA[merA$qual.x > merA$qual.y,"qaccver"]
opt <- opt[opt %in% datFirst$qaccver]
statusList[opt] <- lapply(statusList[opt], c, err["lowerQualA"])
rm(opt)

rm(merA)

# seccond match in B probe
merB <- merge(datAll[datAll$address == "B",], datAll[datAll$address == "B",], by = "qaccver")
merB[merB$qual.x < merB$qual.y, c(grep("\\.x", colnames(merB)), grep("\\.y", colnames(merB)))] <- merB[merB$qual.x < merB$qual.y, c(grep("\\.y", colnames(merB)), grep("\\.x", colnames(merB)))]
merB <- merB[order(merB$qual.x, decreasing = TRUE),]
merB <- merB[apply(merB[,grep("\\.x", colnames(merB))],MARGIN = 1,FUN = paste, collapse = "|") != 
               apply(merB[,grep("\\.y", colnames(merB))],MARGIN = 1,FUN = paste, collapse = "|"),]
merB <- merB[apply(merB[,c("qaccver","saccver.x","pos.x")], MARGIN = 1, FUN = paste, collapse = "|") %in% 
               apply(datFirst[datFirst$address != "S",c("qaccver","saccver","pos")], MARGIN = 1, paste, collapse = "|"),]
merB <- merB[order(merB$qual.y, decreasing = TRUE),]
merB <- merB[!duplicated(merB$qaccver),]

opt <- merB[merB$qual.x == merB$qual.y,"qaccver"]
opt <- opt[opt %in% datFirst$qaccver]
statusList[opt] <- lapply(statusList[opt], c, err["equalQualB"])
rm(opt)

opt <- merB[merB$qual.x > merB$qual.y,"qaccver"]
opt <- opt[opt %in% datFirst$qaccver]
statusList[opt] <- lapply(statusList[opt], c, err["lowerQualB"])
rm(opt)
rm(merB)

# AB comparisons
merAB <- merge(datAll[datAll$address == "A",], datAll[datAll$address == "B",], by = "qaccver", all = TRUE)

opt <- merAB[is.na(merAB$pos.y), "qaccver"]
opt <- opt[opt %in% datFirst$qaccver]
statusList[opt] <- lapply(statusList[opt], c, err["onlyA"])
rm(opt)

opt <- merAB[is.na(merAB$pos.x), "qaccver"]
opt <- opt[opt %in% datFirst$qaccver]
statusList[opt] <- lapply(statusList[opt], c, err["onlyB"])
rm(opt)

merAB <- merAB[complete.cases(merAB),]

a <- merAB[apply(merAB[,c("qaccver","saccver.x","pos.x")], MARGIN = 1, FUN = paste, collapse = "|") %in% 
             apply(datFirst[datFirst$address == "A",c("qaccver","saccver","pos")], MARGIN = 1, paste, collapse = "|"),]
a <- a[!(a$qaccver %in% a[a$saccver.x == a$saccver.x & a$pos.x == a$pos.y,"qaccver"]),]
a <- a[order(a$qual.x, decreasing = TRUE),]
a <- a[order(a$qual.y, decreasing = TRUE),]
a <- a[!duplicated(a$qaccver),]

b <- merAB[apply(merAB[,c("qaccver","saccver.y","pos.y")], MARGIN = 1, FUN = paste, collapse = "|") %in% 
             apply(datFirst[datFirst$address == "B",c("qaccver","saccver","pos")], MARGIN = 1, paste, collapse = "|"),]
b <- b[!(b$qaccver %in% b[b$saccver.x == b$saccver.x & b$pos.x == b$pos.y,"qaccver"]),]
b <- b[order(b$qual.x, decreasing = TRUE),]
b <- b[order(b$qual.y, decreasing = TRUE),]
b <- b[!duplicated(b$qaccver),]

merAB <- rbind(a, b)


# 47 A == B report different position
# 53 A > B, using A
# 59 B > A, using B

opt <- merAB[merAB$qual.x == merAB$qual.y, "qaccver"]
opt <- opt[opt %in% datFirst$qaccver]
statusList[opt] <- lapply(statusList[opt], c, err["diffAB"])
rm(opt)

opt <- merAB[merAB$qual.x > merAB$qual.y, "qaccver"]
opt <- opt[opt %in% datFirst$qaccver]
statusList[opt] <- lapply(statusList[opt], c, err["A>B"])
rm(opt)

opt <- merAB[merAB$qual.x < merAB$qual.y, "qaccver"]
opt <- opt[opt %in% datFirst$qaccver]
statusList[opt] <- lapply(statusList[opt], c, err["A<B"])
rm(opt)

# non Map probes
opt <- newMap$name[is.na(newMap$chr)]
statusList[opt] <- lapply(statusList[opt], c, err["noMatch"])
rm(opt)

# mismatches

opt <- datFirst$qaccver[datFirst$mismatch > 0]
statusList[opt] <- lapply(statusList[opt], c, err["mismatch"])
rm(opt)

# indels
opt <- datFirst$qaccver[datFirst$gapopen > 0]
statusList[opt] <- lapply(statusList[opt], c, err["gap"])
rm(opt)


# for S address
# 61
datFirst$btop <- paste(datFirst$qstart - 1, datFirst$btop, sep = "|")
datShort <- datFirst[datFirst$qstart > 1,]
opt <- datShort$qaccver[datShort$qstart > 1]
statusList[opt] <- lapply(statusList[opt], c, err["5unmapped"])
rm(opt)

datFirst$btop <- paste(datFirst$btop, 50 - datFirst$qend, sep = "|")
datShort <- datFirst[datFirst$qend < 50,]
opt <- datShort$qaccver[50 - datShort$qend > 0]
statusList[opt] <- lapply(statusList[opt], c, err["3unmapped"])
rm(opt)

sList <- sapply(statusList, FUN = prod)

newMap$status <- NA
newMap[names(sList),"status"] <- sList

newMap$BTOP <- NA
newMap[datFirst$qaccver,"BTOP"] <- datFirst$btop


write.table(newMap, file = args[4], 
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


