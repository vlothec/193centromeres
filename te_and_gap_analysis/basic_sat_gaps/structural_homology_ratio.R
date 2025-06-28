#!/bin/Rscript

suppressMessages(library(dplyr))
args <- commandArgs(TRUE)

file.i <- args[1]
print(file.i)

species.i <- gsub(".*/","",file.i) %>%
  {gsub("_.*","",.)}
print(species.i)

data.i <- read.delim(paste0(file.i), header=T, sep = "\t")
# get rid of non centromeric arrays
data.i <- data.i[which(data.i$feature != "long_terminal_repeat" & data.i$feature != "target_site_duplication"),]
# get feature length
data.i$feature.length <- mapply(function(x,y) abs(y-x), data.i$start, data.i$end)
# rename variable atr_6 to make its name meaningful
colnames(data.i)[which(data.i == "homology", arr.ind=T) %>% {.[,2]} %>% unique] <- "ft.method"

# parse centromere occurence info
if(length(grep("TRUE_NA",colnames(data.i))) == 0){
  data.i$TRUE_FALSE.red <- mapply(function(x) ifelse(length(grep("IN",x)) > 0, "IN", "OUT"), data.i$TRUE_FALSE)
}
if(length(grep("TRUE_FALSE.red",colnames(data.i))) > 0 & length(which(is.na(data.i$TRUE_TRUE == F))) > 0){
  data.i$TRUE_TRUE.red <- mapply(function(x) ifelse(length(grep("IN",x)) > 0, "IN", "OUT"), data.i$TRUE_TRUE)
}
if(length(grep("TRUE_FALSE.red",colnames(data.i))) > 0 & length(which(is.na(data.i$TRUE_TRUE == F))) == 0){
  data.i$TRUE_TRUE.red <- NA
}

if(length(grep("TRUE_NA",colnames(data.i))) == 0){
  data.i$filt.collapsed <- mapply(function(x,y) ifelse(x == "IN" | y == "IN", "IN", "OUT"), data.i$TRUE_FALSE.red, data.i$TRUE_TRUE.red)
}

# following if conditions are intended to deal with holocentric
if(length(grep("TRUE_NA",colnames(data.i))) > 0){
  data.i <- data.i[,!grepl("TRUE_FALSE",colnames(data.i))]
}

if(length(which(is.na(data.i$TRUE_NA == F))) > 0 & length(grep("TRUE_NA",colnames(data.i))) > 0){
  data.i$TRUE_NA.red <- mapply(function(x) ifelse(length(grep("IN",x)) > 0, "IN", "OUT"), data.i$TRUE_NA)
}
if(length(which(is.na(data.i$TRUE_NA == F))) == 0 & length(grep("TRUE_NA",colnames(data.i))) > 0){
  data.i$TRUE_NA.red <- NA
}

if(length(grep("TRUE_NA",colnames(data.i))) > 0){
  data.i$filt.collapsed <- data.i$TRUE_NA.red
}

# proper chr list
file.ii <- "proper_chr.txt"
data.ii <- read.delim(paste0(file.ii), header=T, sep = "\t")
colnames(data.ii) <- c("num.id","species","chr","size")
# keep target species
data.ii <- data.ii[which(data.ii$species %in% species.i),]
# keep proper chr
data.i <- data.i[which(data.i$chr %in% unique(data.ii$chr)),]
# infer genome size
gs.i <- data.ii$size %>% sum


if(length(grep("TRUE_NA",colnames(data.i))) == 0){
  tally.data.i.TRUE_FALSE <- data.i %>%
    group_by(TRUE_FALSE.red,ft.method) %>%
    mutate(n.counts = n(),
           nt.sum = sum(feature.length)) %>%
    ungroup() %>%
    select(ft.method, TRUE_FALSE.red, n.counts, nt.sum) %>%
    distinct() %>%
    mutate(prop.n.counts = n.counts/sum(n.counts),
           prop.nt.sum = nt.sum/sum(nt.sum), .by = "TRUE_FALSE.red") %>%
    arrange(TRUE_FALSE.red)
  
  tally.data.i.TRUE_FALSE$gs.i <- gs.i
  tally.data.i.TRUE_FALSE$prop.nt.sum.genome_wide <- tally.data.i.TRUE_FALSE$nt.sum/gs.i
  
  write.table(tally.data.i.TRUE_FALSE,paste0(species.i,"_tally_filt_TRUE_FALSE_structural_homology_ratio"),row.names = F, col.names = T, sep = "\t", quote = F)
  
  
  if(length(which(is.na(data.i$TRUE_TRUE == F))) > 0){
    tally.data.i.TRUE_TRUE <- data.i %>%
      group_by(TRUE_TRUE.red,ft.method) %>%
      mutate(n.counts = n(),
             nt.sum = sum(feature.length)) %>%
      ungroup() %>%
      select(ft.method, TRUE_TRUE.red, n.counts, nt.sum) %>%
      distinct() %>%
      mutate(prop.n.counts = n.counts/sum(n.counts),
             prop.nt.sum = nt.sum/sum(nt.sum), .by = "TRUE_TRUE.red") %>%
      arrange(TRUE_TRUE.red)
    tally.data.i.TRUE_TRUE$gs.i <- gs.i
    tally.data.i.TRUE_TRUE$prop.nt.sum.genome_wide <- tally.data.i.TRUE_TRUE$nt.sum/gs.i
    
    write.table(tally.data.i.TRUE_TRUE,paste0(species.i,"_tally_filt_TRUE_TRUE_structural_homology_ratio"),row.names = F, col.names = T, sep = "\t", quote = F)
  }else{
    print(paste0(species.i,": no TRUE_TRUE obs found"))
  }
  
  tally.data.i.filt.collapsed <- data.i %>%
    group_by(filt.collapsed,ft.method) %>%
    mutate(n.counts = n(),
           nt.sum = sum(feature.length)) %>%
    ungroup() %>%
    select(ft.method, filt.collapsed, n.counts, nt.sum) %>%
    distinct() %>%
    mutate(prop.n.counts = n.counts/sum(n.counts),
           prop.nt.sum = nt.sum/sum(nt.sum), .by = "filt.collapsed") %>%
    arrange(filt.collapsed)
  tally.data.i.filt.collapsed$gs.i <- gs.i
  tally.data.i.filt.collapsed$prop.nt.sum.genome_wide <- tally.data.i.filt.collapsed$nt.sum/gs.i
  
  write.table(tally.data.i.filt.collapsed,paste0(species.i,"_tally_filt_collapsed_structural_homology_ratio"),row.names = F, col.names = T, sep = "\t", quote = F)
}


if(length(grep("TRUE_NA",colnames(data.i))) > 0){
  tally.data.i.TRUE_NA <- data.i %>%
    group_by(TRUE_NA.red,ft.method) %>%
    mutate(n.counts = n(),
           nt.sum = sum(feature.length)) %>%
    ungroup() %>%
    select(ft.method, TRUE_NA.red, n.counts, nt.sum) %>%
    distinct() %>%
    mutate(prop.n.counts = n.counts/sum(n.counts),
           prop.nt.sum = nt.sum/sum(nt.sum), .by = "TRUE_NA.red") %>%
    arrange(TRUE_NA.red)
  
  tally.data.i.TRUE_NA$gs.i <- gs.i
  tally.data.i.TRUE_NA$prop.nt.sum.genome_wide <- tally.data.i.TRUE_NA$nt.sum/gs.i
  
  write.table(tally.data.i.TRUE_NA,paste0(species.i,"_tally_filt_TRUE_NA_structural_homology_ratio"),row.names = F, col.names = T, sep = "\t", quote = F)

}
