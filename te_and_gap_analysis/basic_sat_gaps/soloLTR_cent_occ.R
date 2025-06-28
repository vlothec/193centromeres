#!/bin/Rscript

suppressMessages(library(dplyr))
options(scipen = 999)
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


data.i <- read.table(file.i,comment.char="",sep="\t")
colnames(data.i) <- c("chr","method","feature","start","end","score","strand","phase","attributes")

feature.length <- mapply(function (x,y) (y-x)+1,data.i$start,data.i$end)
data.i$feature.length <- feature.length

# load centromere boundaries
{
  cent.files <- file.ii
  for (cent in cent.files){
    cent.data.i <-read.table(cent, header=F, sep="\t")
    colnames(cent.data.i) <- c("chr","start","end","feature","starnd","filt.1","filt.2")
    
    if(nrow(cent.data.i[which(cent.data.i$filt.1 == TRUE & cent.data.i$filt.2 == FALSE),c(1,2,3,4)]) > 0){
      cent.data.i <- cent.data.i[which(cent.data.i$filt.1 == TRUE & cent.data.i$filt.2 == FALSE),c(1,2,3,4)]
      cent.data.i$feature <- "array"
      
    }
    
  }
  
  data.ii <- cent.data.i
  colnames(data.ii) <- c("chr","start","end")
  colnames(data.ii)[c(2,3)] <- c("start.cent","end.cent")
  
  data.ii <- data.ii[which(data.ii$end.cent > 0),]
  
  data.ii$length.cent <- mapply(function(x,y) (y-x)+1, data.ii$start.cent, data.ii$end.cent)
  data.i <- data.i[which(data.i$chr %in% unique(data.ii$chr)),]
  
  data.i$cent.occ <- "OUT"
  for(chr in 1:length(unique(data.ii$chr))){
    chr.i <- unique(data.ii$chr)[chr]
    cent.data.ii <- data.ii[which(cent.data.i$chr %in% chr.i),]
    cent.data.ii <- cent.data.ii[with(cent.data.ii,order(start.cent)),]
    for(arr in 1:nrow(cent.data.ii)){
      data.i$cent.occ[which(data.i$chr == chr.i & data.i$cent.occ == "OUT" & data.i$start >= cent.data.ii$start[arr] & data.i$end <= cent.data.ii$end[arr])] <- "IN"
    }
  }

  table(data.i$cent.occ)
  
  data.i <- data.i[which(data.i$chr %in% unique(data.ii$chr)),]
  
}

data.i %>%
  group_by(cent.occ) %>%
  summarise(n.counts = n(),
            sum.length = sum(feature.length))

fam.cls <- mapply(function (attributes) gsub(".*;Lineage=([A-Za-z_]+.*);.*;.*","\\1",attributes),data.i$attributes)
data.i$fam.cls <- fam.cls

data.i.sum <- data.i %>%
  group_by(fam.cls,cent.occ) %>%
  summarise(n.counts = n(),
            sum.length = sum(feature.length)) %>%
  arrange(desc(sum.length)) %>%
  print(n = 22)

cent.space <- data.ii$length.cent %>% sum
genome.space <- data.i %>%
  group_by(chr) %>%
  summarise(max.end = max(end)) %>%
  select(max.end) %>%
  sum

data.i.sum$space.prop <- mapply(function (cent.occ,sum.length,genome.space,cent.space) ifelse(cent.occ == FALSE,sum.length/(genome.space-cent.space),sum.length/cent.space),data.i.sum$cent.occ,data.i.sum$sum.length,genome.space,cent.space)

sink(paste0(file.i.name,".cent.sum"))
print(data.i.sum,n = nrow(data.i.sum))
sink()

write.table(data.i,paste0(file.i.name,"_wga_soloLTR.gff3.cent_occurrence"),col.names = T,row.names = F, quote = F,sep = "\t")
