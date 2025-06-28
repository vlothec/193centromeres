#!/bin/Rscript

suppressMessages(library(dplyr))
options(scipen = 999)
args <- commandArgs(TRUE)

species.i <- args[1]


df.complete <- data.frame()

data.i <- try(read.table(paste0(species.i,"_all_in_gaps_solo_list.txt"),comment.char="",sep=" "))

if(class(data.i) != "try-error"){
  colnames(data.i) <- c("gap.id","fam","ft.coord","gap.length","feature.length","overlap")
  
  data.i$fam <- gsub("soloLTR#","",data.i$fam)
  data.i.sum <- data.i %>%
    group_by(fam) %>%
    summarise(species = species.i,
              ft = "soloLTR",
              th = "all_in_gaps",
              n.counts = n(),
              nt.sum = sum(feature.length))
  
  df.complete <- rbind(df.complete,data.i.sum)
  rm(data.i.sum, data.i)
}

# soloLTR all

data.i <- try(read.table(paste0(species.i,"_all_solo_list.txt"),comment.char="",sep=" "))

if(class(data.i) != "try-error"){
  colnames(data.i) <- c("gap.id","fam","ft.coord","gap.length","feature.length","overlap")
  
  data.i$fam <- gsub("soloLTR#","",data.i$fam)
  data.i.sum <- data.i %>%
    group_by(fam) %>%
    summarise(species = species.i,
              ft = "soloLTR",
              th = "all_in_overlap",
              n.counts = n(),
              nt.sum = sum(feature.length))
  
  df.complete <- rbind(df.complete,data.i.sum)
  rm(data.i.sum, data.i)
}

# soloLTR 100 bp

data.i <- try(read.table(paste0(species.i,"_solo_list_100_bp.txt"),comment.char="",sep=" "))

if(class(data.i) != "try-error"){
  colnames(data.i) <- c("gap.id","fam","ft.coord","gap.length","feature.length","overlap")

  data.i$fam <- gsub("soloLTR#","",data.i$fam)
  data.i.sum <- data.i %>%
    group_by(fam) %>%
    summarise(species = species.i,
              ft = "soloLTR",
              th = "all_in_overlap_100",
              n.counts = n(),
              nt.sum = sum(feature.length))
  
  df.complete <- rbind(df.complete,data.i.sum)
  rm(data.i.sum, data.i)
}

# intact all in gaps

data.i <- try(read.table(paste0(species.i,"_all_in_gaps_intact_list.txt"),comment.char="",sep=" "))

if(class(data.i) != "try-error"){
  colnames(data.i) <- c("gap.id","fam","ft.coord","gap.length","feature.length","overlap")
  
  data.i.sum <- data.i %>%
    group_by(fam) %>%
    summarise(species = species.i,
              ft = "intact",
              th = "all_in_gaps",
              n.counts = n(),
              nt.sum = sum(feature.length))
  
  df.complete <- rbind(df.complete,data.i.sum)
  rm(data.i.sum, data.i)
}


# intact all

data.i <- try(read.table(paste0(species.i,"_all_intact_list.txt"),comment.char="",sep=" "))

if(class(data.i) != "try-error"){
  colnames(data.i) <- c("gap.id","fam","ft.coord","gap.length","feature.length","overlap")
  
  data.i.sum <- data.i %>%
    group_by(fam) %>%
    summarise(species = species.i,
              ft = "intact",
              th = "all_in_overlap",
              n.counts = n(),
              nt.sum = sum(feature.length))
  
  df.complete <- rbind(df.complete,data.i.sum)
  rm(data.i.sum, data.i)
}

# intact 200 bp

data.i <- try(read.table(paste0(species.i,"_intact_list_200_bp.txt"),comment.char="",sep=" "))

if(class(data.i) != "try-error"){
  colnames(data.i) <- c("gap.id","fam","ft.coord","gap.length","feature.length","overlap")
  
  data.i.sum <- data.i %>%
    group_by(fam) %>%
    summarise(species = species.i,
              ft = "intact",
              th = "all_in_overlap_200",
              n.counts = n(),
              nt.sum = sum(feature.length))
  
  df.complete <- rbind(df.complete,data.i.sum)
  rm(data.i.sum, data.i)
}

write.table(df.complete[,c(2,1,3,4,5,6)],paste0(species.i,"_fam_tally.txt"),sep = "\t", quote = F, col.names = T, row.names = F)
