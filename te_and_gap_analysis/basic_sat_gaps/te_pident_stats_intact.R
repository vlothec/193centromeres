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
data.i <- data.i[which(data.i$feature != "long_terminal_repeat" & data.i$feature != "target_site_duplication" & data.i$feature != "repeat_region"),]
# get feature length
data.i$feature.length <- mapply(function(x,y) abs(y-x)+1, data.i$start, data.i$end)

# get pident col
data.i$pident <- mapply(function(x) gsub(".*dentity=([0-9\\.]+|NA);.*","\\1",x),data.i$atr) %>%
  as.numeric()
data.i$ft.method <- mapply(function(x) gsub(".*Method=(homology|structural).*","\\1",x),data.i$atr)


# TEs with no pident data
print(data.i$feature[which(data.i$ft.method == "structural" & is.na(data.i$pident) == T)] %>% table)

data.i <- data.i[which(is.na(data.i$pident) == F),]
print(data.i$ft.method %>% table)
print(data.i$feature %>% table)

data.i <- data.i %>%
  mutate(feature.0 = case_when(grepl("LTR",feature) & !grepl("non",feature) ~ "LTRRT",
                               grepl("TIR_transposon",feature) ~ "TIR",
                               .default = "unk"))
print(data.i$feature.0 %>% table)

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

# data.i <- data.i[which(data.i$ft.method == "structural"),]
# data.i %>% nrow()

counts.f <- function(x){
  hist(x, breaks=seq(0,1,0.001), plot=F)$counts %>%
    {paste0(.,collapse = ";")}
}
mids.f <- function(x){
  hist(x, breaks=seq(0,1,0.001), plot=F)$mids %>%
    {paste0(.,collapse = ";")}
}

if(length(grep("TRUE_NA",colnames(data.i))) == 0){
  tally.data.i.TRUE_FALSE <- data.i %>%
    group_by(TRUE_FALSE.red,ft.method,feature.0) %>%
    summarise(obs.values = paste0(pident, collapse = ";"),
              counts.hist = counts.f(pident),
              mids.hist = mids.f(pident),
              lower.whisker.bp = boxplot.stats(pident)$stats[1],
              lower.hinge.bp = boxplot.stats(pident)$stats[2],
              median.bp = boxplot.stats(pident)$stats[3],
              upper.hinge.bp = boxplot.stats(pident)$stats[4],
              upper.whisker.bp = boxplot.stats(pident)$stats[5],
              n.obs = boxplot.stats(pident)$n)
  tally.data.i.TRUE_FALSE$gs.i <- gs.i
  
  write.table(tally.data.i.TRUE_FALSE,paste0(species.i,"_tally_filt_TRUE_FALSE_pident"),row.names = F, col.names = T, sep = "\t", quote = F)
  
  # all
  tally.data.i.TRUE_FALSE <- data.i %>%
    group_by(TRUE_FALSE.red,ft.method) %>%
    summarise(obs.values = paste0(pident, collapse = ";"),
              counts.hist = counts.f(pident),
              mids.hist = mids.f(pident),
              lower.whisker.bp = boxplot.stats(pident)$stats[1],
              lower.hinge.bp = boxplot.stats(pident)$stats[2],
              median.bp = boxplot.stats(pident)$stats[3],
              upper.hinge.bp = boxplot.stats(pident)$stats[4],
              upper.whisker.bp = boxplot.stats(pident)$stats[5],
              n.obs = boxplot.stats(pident)$n)
  tally.data.i.TRUE_FALSE$gs.i <- gs.i
  
  write.table(tally.data.i.TRUE_FALSE,paste0(species.i,"_tally_filt_TRUE_FALSE_all_pident"),row.names = F, col.names = T, sep = "\t", quote = F)
  
  if(length(which(is.na(data.i$TRUE_TRUE == F))) > 0){
    tally.data.i.TRUE_TRUE <- data.i %>%
      group_by(TRUE_TRUE.red,ft.method,feature.0) %>%
      summarise(obs.values = paste0(pident, collapse = ";"),
                counts.hist = counts.f(pident),
                mids.hist = mids.f(pident),
                lower.whisker.bp = boxplot.stats(pident)$stats[1],
                lower.hinge.bp = boxplot.stats(pident)$stats[2],
                median.bp = boxplot.stats(pident)$stats[3],
                upper.hinge.bp = boxplot.stats(pident)$stats[4],
                upper.whisker.bp = boxplot.stats(pident)$stats[5],
                n.obs = boxplot.stats(pident)$n)
    tally.data.i.TRUE_TRUE$gs.i <- gs.i
    
    write.table(tally.data.i.TRUE_TRUE,paste0(species.i,"_tally_filt_TRUE_TRUE_pident"),row.names = F, col.names = T, sep = "\t", quote = F)
    # all
    tally.data.i.TRUE_TRUE <- data.i %>%
      group_by(TRUE_TRUE.red,ft.method) %>%
      summarise(obs.values = paste0(pident, collapse = ";"),
                counts.hist = counts.f(pident),
                mids.hist = mids.f(pident),
                lower.whisker.bp = boxplot.stats(pident)$stats[1],
                lower.hinge.bp = boxplot.stats(pident)$stats[2],
                median.bp = boxplot.stats(pident)$stats[3],
                upper.hinge.bp = boxplot.stats(pident)$stats[4],
                upper.whisker.bp = boxplot.stats(pident)$stats[5],
                n.obs = boxplot.stats(pident)$n)
    tally.data.i.TRUE_TRUE$gs.i <- gs.i
    
    write.table(tally.data.i.TRUE_TRUE,paste0(species.i,"_tally_filt_TRUE_TRUE_all_pident"),row.names = F, col.names = T, sep = "\t", quote = F)
    
  }else{
    print(paste0(species.i,": no TRUE_TRUE obs found"))
  }
  
  tally.data.i.filt.collapsed <- data.i %>%
    group_by(filt.collapsed,ft.method,feature.0) %>%
    summarise(obs.values = paste0(pident, collapse = ";"),
              counts.hist = counts.f(pident),
              mids.hist = mids.f(pident),
              lower.whisker.bp = boxplot.stats(pident)$stats[1],
              lower.hinge.bp = boxplot.stats(pident)$stats[2],
              median.bp = boxplot.stats(pident)$stats[3],
              upper.hinge.bp = boxplot.stats(pident)$stats[4],
              upper.whisker.bp = boxplot.stats(pident)$stats[5],
              n.obs = boxplot.stats(pident)$n)
  tally.data.i.filt.collapsed$gs.i <- gs.i
  
  write.table(tally.data.i.filt.collapsed,paste0(species.i,"_tally_filt_collapsed_pident"),row.names = F, col.names = T, sep = "\t", quote = F)
  
  # all
  tally.data.i.filt.collapsed <- data.i %>%
    group_by(filt.collapsed,ft.method) %>%
    summarise(obs.values = paste0(pident, collapse = ";"),
              counts.hist = counts.f(pident),
              mids.hist = mids.f(pident),
              lower.whisker.bp = boxplot.stats(pident)$stats[1],
              lower.hinge.bp = boxplot.stats(pident)$stats[2],
              median.bp = boxplot.stats(pident)$stats[3],
              upper.hinge.bp = boxplot.stats(pident)$stats[4],
              upper.whisker.bp = boxplot.stats(pident)$stats[5],
              n.obs = boxplot.stats(pident)$n)
  tally.data.i.filt.collapsed$gs.i <- gs.i
  
  write.table(tally.data.i.filt.collapsed,paste0(species.i,"_tally_filt_collapsed_all_pident"),row.names = F, col.names = T, sep = "\t", quote = F)
  
}

if(length(grep("TRUE_NA",colnames(data.i))) > 0){
  tally.data.i.TRUE_NA <- data.i %>%
    group_by(TRUE_NA.red,ft.method,feature.0) %>%
    summarise(obs.values = paste0(pident, collapse = ";"),
              counts.hist = counts.f(pident),
              mids.hist = mids.f(pident),
              lower.whisker.bp = boxplot.stats(pident)$stats[1],
              lower.hinge.bp = boxplot.stats(pident)$stats[2],
              median.bp = boxplot.stats(pident)$stats[3],
              upper.hinge.bp = boxplot.stats(pident)$stats[4],
              upper.whisker.bp = boxplot.stats(pident)$stats[5],
              n.obs = boxplot.stats(pident)$n)
  
  tally.data.i.TRUE_NA$gs.i <- gs.i
  
  write.table(tally.data.i.TRUE_NA,paste0(species.i,"_tally_filt_TRUE_NA_pident"),row.names = F, col.names = T, sep = "\t", quote = F)
  
  # all
  tally.data.i.TRUE_NA <- data.i %>%
    group_by(TRUE_NA.red,ft.method) %>%
    summarise(obs.values = paste0(pident, collapse = ";"),
              counts.hist = counts.f(pident),
              mids.hist = mids.f(pident),
              lower.whisker.bp = boxplot.stats(pident)$stats[1],
              lower.hinge.bp = boxplot.stats(pident)$stats[2],
              median.bp = boxplot.stats(pident)$stats[3],
              upper.hinge.bp = boxplot.stats(pident)$stats[4],
              upper.whisker.bp = boxplot.stats(pident)$stats[5],
              n.obs = boxplot.stats(pident)$n)
  
  tally.data.i.TRUE_NA$gs.i <- gs.i
  
  write.table(tally.data.i.TRUE_NA,paste0(species.i,"_tally_filt_TRUE_NA_all_pident"),row.names = F, col.names = T, sep = "\t", quote = F)
  
}
