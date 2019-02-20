rm(list = ls())
options(stringsAsFactors = FALSE)


args = commandArgs(trailingOnly=TRUE)

dat <- read.csv(args[1], fill = TRUE, skip = 7)


df <- data.frame(paste(">", dat$Name, sep = ""), dat$AlleleA_ProbeSeq)
df <- df[df[,2] != "",]

write.table(df, file = paste(args[2],"Aseq.fa", sep = "_"), 
            sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)

df <- data.frame(paste(">", dat$Name, sep = ""), dat$AlleleB_ProbeSeq)
df <- df[df[,2] != "",]

write.table(df, file = paste(args[2],"Bseq.fa", sep = "_"), 
            sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)
