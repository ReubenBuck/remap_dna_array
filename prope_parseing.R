rm(list = ls())
options(stringsAsFactors = FALSE)

dat <- read.table("~/Desktop/hills_remap/out/250k.blast_out.txt", sep = "\t",
                  col.names = c("qaccver", "saccver", "pident", "length", "mismatch",
                                "gapopen", "qstart", "qend", "sstart",
                                "send", "evalue", "bitscore", "btop"))

dat$strand = "+"
dat$strand[dat$send < dat$sstart] = "-"
dat[dat$strand == "-", c("sstart","send")] <- dat[dat$strand == "-", c("send","sstart")]


front <- dat[grep("front", dat$qaccver),]
front$qaccver <- gsub("_front", "", front$qaccver)
front$side = "front"

back <- dat[grep("back", dat$qaccver),]
back$qaccver <- gsub("_back", "", back$qaccver)
back$side = "back"


dupProbesFront <- unique(front$qaccver[queryHits(ol[duplicated(queryHits(ol))])])

,subjectHits(ol[duplicated(subjectHits(ol))])))


# all sequences where chr is the same
# all sequences where there is a gap of 1 between them

library(GenomicRanges)

front.gr <- GRanges(seqnames = front$saccver, 
                    ranges = IRanges(start = front$send + 1, width = 1))
front.gr$start[front$strand == "-"] <- front[front$strand == "-", "sstart"] - 1
front.gr$end[front$strand == "-"] <- front[front$strand == "-", "sstart"] - 1

back.gr <- GRanges(seqnames = back$saccver, 
                    ranges = IRanges(start = back$sstart - 1, width = 1))
back.gr$start[back$strand == "-"] <- back[back$strand == "-", "send"] + 1
back.gr$end[back$strand == "-"] <- back[back$strand == "-", "send"] + 1

ol <- findOverlaps(front.gr, back.gr)

ol[duplicated(queryHits(ol))]

# A1_117677838-75_T_F_2301561504, both flanking regions have 2 perfect matches in the genome




head(front)

mer <- merge(front,back, by = 1)

pmax(mer[,c("sstart.x", "sstart.y", "send.x", "send.y")]) - pmin

pMatch <- dat[dat$btop == 75,]
npMatch <- dat[dat$btop < 75,]


sum(duplicated(pMatch$qaccver))


# work our way down the list!


# perferct matches
# both sides
# duplicates
# orientation






