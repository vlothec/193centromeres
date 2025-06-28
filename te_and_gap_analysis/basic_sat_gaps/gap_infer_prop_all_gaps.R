#!/bin/Rscript

suppressMessages(library(dplyr))
args <- commandArgs(TRUE)

file.i <- args[1]

print(file.i)

data.i <- read.delim(paste0(getwd(),"/",file.i), header=FALSE, comment.char = "")
colnames(data.i) <- c("TRASH.chr","TRASH.start","TRASH.end","TRASH.class","TRASH.strand","GAP.length.id","EDTA.chr","EDTA.start","EDTA.end","EDTA.feature","EDTA.score","EDTA.strand","overlap")

data.i$GAP.length <- data.i$GAP.length.id %>%
  {gsub(".*_","",.)} %>%
  {as.numeric(.)}

n.rejected.gaps <- which(data.i$GAP.length > 100000 & data.i$EDTA.feature %in% c("target_site_duplication","long_terminal_repeat")) %>%
  length()

print(length(which(data.i$EDTA.feature == "repeat_region")))

data.i <- data.i[which(data.i$GAP.length < 100000),]

data.i <- data.i[!grepl("target_site_duplication|long_terminal_repeat",data.i$EDTA.feature),]

# additional vars
data.i <- data.i %>%
  mutate(GAP.type = gsub(".*repeats_([A-Z_A-Z]+).bed.gaps.TEanno.*","\\1",file.i),
         EDTA.length = (EDTA.end-EDTA.start),
         prop.GAP.covered.raw = overlap/GAP.length,
         prop.EDTA.covered.raw = overlap/((EDTA.end-EDTA.start))) %>%
  group_by(GAP.length.id) %>%
  mutate(prop.GAP.covered = sum(overlap)/unique(GAP.length),
         prop.EDTA.covered = sum(overlap)/sum(EDTA.length),
         prop.GAP.unannotated = ifelse((GAP.length - sum(overlap)) > 0, (GAP.length - sum(overlap))/GAP.length, 0))

write.table(data.i,paste0(getwd(),"/",file.i,".parsed"), row.names = F, col.names = T, sep = "\t", quote=F)
