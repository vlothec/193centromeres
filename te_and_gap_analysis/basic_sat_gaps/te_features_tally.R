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
data.i$TRUE_FALSE.red <- mapply(function(x) ifelse(length(grep("IN",x)) > 0, "IN", "OUT"), data.i$TRUE_FALSE)
if(length(which(is.na(data.i$TRUE_TRUE == F))) > 0){
  data.i$TRUE_TRUE.red <- mapply(function(x) ifelse(length(grep("IN",x)) > 0, "IN", "OUT"), data.i$TRUE_TRUE)
}else{
  data.i$TRUE_TRUE.red <- NA
}
data.i$filt.collapsed <- mapply(function(x,y) ifelse(x == "IN" | y == "IN", "IN", "OUT"), data.i$TRUE_FALSE.red, data.i$TRUE_TRUE.red)


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

# load TE_general_cls
file.iii <-"TE_general_cls"
data.iii <- read.delim(paste0(file.iii), header=F, sep = " ")
colnames(data.iii) <- c("cls.0","feature")

prev.nrow <- data.i %>% nrow
data.i <- left_join(data.i,data.iii,by = "feature")
post.nrow <- data.i %>% nrow
tally.per.cls.0 <- data.i %>%
  group_by(cls.0) %>%
  summarise(n.counts = n()) %>%
  ungroup() %>%
  mutate(sum.n.counts = sum(n.counts))
tally.per.cls.0$prev.nrow <- prev.nrow
tally.per.cls.0$post.nrow <- post.nrow

tally.per.cls.0$match.test <- mapply(function(x,y,z) ifelse(x == y & x == z, T, F), tally.per.cls.0$sum.n.counts, tally.per.cls.0$prev.nrow, tally.per.cls.0$post.nrow)
print(tally.per.cls.0)


tally.data.i.total <- data.i %>%
  group_by(feature) %>%
  mutate(n.counts = n(),
         nt.sum = sum(feature.length)) %>%
  ungroup() %>%
  select(cls.0, feature, n.counts, n.counts, nt.sum) %>%
  distinct() %>%
  mutate(prop.n.counts = n.counts/sum(n.counts),
         prop.nt.sum = nt.sum/sum(nt.sum))
tally.data.i.total$gs.i <- gs.i
tally.data.i.total$prop.nt.sum.genome_wide <- tally.data.i.total$nt.sum/gs.i

write.table(tally.data.i.total,paste0(species.i,"_tally_total"),row.names = F, col.names = T, sep = "\t", quote = F)

tally.data.i.TRUE_FALSE <- data.i %>%
  group_by(feature,TRUE_FALSE.red) %>%
  mutate(n.counts = n(),
         nt.sum = sum(feature.length)) %>%
  ungroup() %>%
  select(cls.0, feature, TRUE_FALSE.red, n.counts, nt.sum) %>%
  distinct() %>%
  mutate(prop.n.counts = n.counts/sum(n.counts),
         prop.nt.sum = nt.sum/sum(nt.sum))
tally.data.i.TRUE_FALSE$gs.i <- gs.i
tally.data.i.TRUE_FALSE$prop.nt.sum.genome_wide <- tally.data.i.TRUE_FALSE$nt.sum/gs.i

write.table(tally.data.i.TRUE_FALSE,paste0(species.i,"_tally_filt_TRUE_FALSE"),row.names = F, col.names = T, sep = "\t", quote = F)


if(length(which(is.na(data.i$TRUE_TRUE == F))) > 0){
  tally.data.i.TRUE_TRUE <- data.i %>%
    group_by(feature,TRUE_TRUE.red) %>%
    mutate(n.counts = n(),
           nt.sum = sum(feature.length)) %>%
    ungroup() %>%
    select(cls.0, feature, TRUE_TRUE.red,n.counts, nt.sum) %>%
    distinct() %>%
    mutate(prop.n.counts = n.counts/sum(n.counts),
           prop.nt.sum = nt.sum/sum(nt.sum))
  tally.data.i.TRUE_TRUE$gs.i <- gs.i
  tally.data.i.TRUE_TRUE$prop.nt.sum.genome_wide <- tally.data.i.TRUE_TRUE$nt.sum/gs.i
  
  write.table(tally.data.i.TRUE_TRUE,paste0(species.i,"_tally_filt_TRUE_TRUE"),row.names = F, col.names = T, sep = "\t", quote = F)
}else{
  print(paste0(species.i,": no TRUE_TRUE obs found"))
}

tally.data.i.filt.collapsed <- data.i %>%
  group_by(feature,filt.collapsed) %>%
  mutate(n.counts = n(),
         nt.sum = sum(feature.length)) %>%
  ungroup() %>%
  select(cls.0, feature, filt.collapsed,n.counts, nt.sum) %>%
  distinct() %>%
  mutate(prop.n.counts = n.counts/sum(n.counts),
         prop.nt.sum = nt.sum/sum(nt.sum))
tally.data.i.filt.collapsed$gs.i <- gs.i
tally.data.i.filt.collapsed$prop.nt.sum.genome_wide <- tally.data.i.filt.collapsed$nt.sum/gs.i

write.table(tally.data.i.filt.collapsed,paste0(species.i,"_tally_filt_collapsed"),row.names = F, col.names = T, sep = "\t", quote = F)
