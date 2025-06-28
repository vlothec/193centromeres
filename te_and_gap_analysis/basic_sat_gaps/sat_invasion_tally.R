#!/bin/Rscript

suppressMessages(library(dplyr))
args <- commandArgs(TRUE)

file.i <- args[1]
print(file.i)
file.i.name <- file.i %>%
  {gsub(".*/","",.)}

file.ii <- args[2]
print(file.ii)
file.ii.name <- file.ii %>%
  {gsub(".*/","",.)}

species.i <- gsub("\\..*","",file.ii.name) 
print(species.i)


data.i <- read.table(file.i, header = T)
cent.occ <- mapply(function(cent_col) ifelse(length(grep("IN",cent_col)) > 0, "IN","OUT"),data.i$TRUE_FALSE)
data.i$cent.occ <- cent.occ
table(data.i$cent.occ)

nrow(data.i)
data.i <- data.i[!grepl("long_terminal_repeat|target_site_duplication",data.i$feature),]
nrow(data.i)
data.i <- data.i[!grepl("ID=repeat_region_",data.i$atr),]
nrow(data.i)


data.ii <- read.table(file.ii, header = T, comment.char = "",sep = "\t")

data.ii$feature <- mapply(function(x) gsub(".*#","",x),data.ii$X.TE)
data.ii$chr <- mapply(function(x) gsub("(.*):.*","\\1",x),data.ii$X.TE)
data.ii$start <- mapply(function(x) gsub(".*:([0-9]+)\\.\\..*#.*","\\1",x),data.ii$X.TE) %>% as.numeric()
data.ii$end <- mapply(function(x) gsub(".*:.*\\.\\.([0-9]+)-.*#.*","\\1",x),data.ii$X.TE) %>% as.numeric()
data.ii$method <- mapply(function(x) gsub(".*:.*\\.\\.[0-9]+-(.*)#.*","\\1",x),data.ii$X.TE)
data.ii$n.hmm <- mapply(function(x) length(unlist(gregexpr("\\|",x))),data.ii$Domains)

# clean the data little bit
nrow(data.ii)

test.consistency <- mapply(function(x,y) ifelse(length(grep(x,y)) > 0,T,F),data.ii$Order,data.ii$feature)
data.ii <- data.ii[test.consistency,]
nrow(data.ii)

data.ii <- data.ii[which(data.ii$Domains != "none"),]
nrow(data.ii)

data.ii <- data.ii[which(data.ii$Clade != "unknown"),]
nrow(data.ii)


data.iv <- left_join(data.i,data.ii,by = c("chr","start","end"))
nrow(data.iv)

length(which(is.na(data.iv$Clade)==F))
data.iv <- data.iv[which(is.na(data.iv$Clade)==F),]

data.iv$ft.length <- mapply(function(x,y) (y-x)+1, data.iv$start, data.iv$end)

write.table(data.iv,paste0(getwd(),"/",file.ii.name,".anno"),col.names = T,row.names = F,quote = F,sep = "\t")

data.iv.sum <- data.iv %>%
  group_by(Clade,cent.occ) %>%
  summarise(n.count = n(),
            nt.sum = sum(ft.length),
            species = species.i) %>%
  arrange(cent.occ,nt.sum) %>%
  as.data.frame()

write.table(data.iv.sum,paste0(getwd(),"/",file.ii.name,".tally"),col.names = T,row.names = F,quote = F,sep = "\t")
