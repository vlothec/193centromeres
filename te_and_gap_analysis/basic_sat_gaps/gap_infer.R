#!/bin/Rscript

suppressMessages(library(dplyr))
args <- commandArgs(TRUE)

file.i <- args[1]

print(file.i)

data.i <- read.delim(paste0(getwd(),"/",file.i), header=FALSE, comment.char = "")
colnames(data.i) <- c("chr","start","end","TRASH_class","strand")

data.i <- data.i[with(data.i, order(chr,start)),]
data.i <- data.i %>%
  group_by(chr) %>%
  mutate(gap.start = lag(end),
         gap.end = start) %>%
  mutate(gap.length = gap.end-gap.start)

data.i <- data.i[which(data.i$gap.length > 0),]
data.i <- data.i[,c(1,6,7,4,5,8)]

write.table(data.i,paste0(getwd(),"/",file.i,".gaps"), row.names = F, col.names = T, sep = "\t", quote=F)
