#!/bin/Rscript

suppressMessages(library(dplyr))
args <- commandArgs(TRUE)

file.i <- args[1]
print(file.i)

species.i <- gsub(".*/","",file.i) %>%
  {gsub("_.*","",.)}

data.i <- read.delim(paste0(file.i), header=FALSE, sep = "\t")
colnames(data.i) <- c("chr","start","end","TRASH_class","strand","filt_1","filt_2")
# get rid of non centromeric arrays
data.i <- data.i[which(data.i$filt_1 != "FALSE"),]
data.i$filt_comb <- mapply(function(x,y) paste0(x,"/",y), data.i$filt_1, data.i$filt_2)
data.i$filt_comb_row <- mapply(function(x,y,z) paste0(x,"/",y,"/",z), data.i$filt_1, data.i$filt_2, rownames(data.i))
cent.type <- mapply(function(x) ifelse(length(grep("NA",x)) > 0, "holocentrict","monocentric"), data.i$filt_comb_row) %>%
  unique()


file.ii <- args[2]
print(file.ii)

if(length(grep("intact.gff3",file.ii)) > 0){
  data.ii <- read.delim(paste0(file.ii), header=FALSE, sep = "\t")
  head(data.ii)
  colnames(data.ii) <- c("chr","method","feature","start","end","score","strand","phase","atr")
  for(i in 1:length(seq(10,24,1))){
    name.i <- seq(10,24,1)[i]
    data.ii <- cbind(data.ii,NA)
    colnames(data.ii)[ncol(data.ii)] <- paste0("atr_", name.i)
  }
}

if(length(grep("_edta_filtered.csv.reassigned",file.ii)) > 0){
  data.ii <- read.delim(paste0(file.ii), header=FALSE, sep = ",")
  head(data.ii)
  colnames(data.ii) <- c("row","chr","method","feature","start","end","score","strand","phase",paste0("atr_",seq(1,ncol(data.ii)-9,1)))
}


cent.occurence.f <- function(te.start,array.start,array.end){
  ifelse(te.start >= array.start & te.start <= array.end, "IN","OUT")
}

for(chr in 1:length(unique(data.ii$chr))){
  chr.i <- unique(data.ii$chr)[chr]
  # subset te annotation
  data.ii.chr.i <- data.ii[which(data.ii$chr %in% chr.i),]
  
  if(length(which(data.i$chr %in% chr.i)) > 0){
    # subset array annotation
    data.i.chr.i <- data.i[which(data.i$chr %in% chr.i),]
    
    df <- data.frame()
    df <- df[seq(1,nrow(data.ii.chr.i),1),FALSE]
    filt.v <- vector()
    for(row in 1:nrow(data.i.chr.i)){
      filt.v <- c(filt.v,data.i.chr.i$filt_comb_row[row])
      
      df.i <- data.frame("A" = mapply(cent.occurence.f, data.ii.chr.i$start, data.i.chr.i$start[row], data.i.chr.i$end[row]))
      colnames(df.i) <- data.i.chr.i$filt_comb_row[row]
      df <- cbind(df,df.i)
    }
    
    filt.v.red <- gsub("([A-Z]{4,5}/[A-Z]{2,5}).*","\\1",filt.v) %>%
      unique() %>%
      sort
    for(filt in 1:length(unique(filt.v.red))){
      filt.i <- filt.v.red[filt]
      if(ncol(data.frame(df[,grep(filt.i,filt.v)])) > 1){
        data.ii.chr.i$filt.collapsed <- apply(df[,grep(filt.i,filt.v)],1,function(x) paste(x,collapse = "/"))
        colnames(data.ii.chr.i)[ncol(data.ii.chr.i)] <- filt.i %>%
          {gsub("/","_",.)}
        assign(paste0("data.chr.",chr),data.ii.chr.i)
      }else{
        data.ii.chr.i$filt.collapsed <- df[,grep(filt.i,filt.v)]
        colnames(data.ii.chr.i)[ncol(data.ii.chr.i)] <- filt.i %>%
          {gsub("/","_",.)}
        assign(paste0("data.chr.",chr),data.ii.chr.i)
      }
    }
  }else{
    assign(paste0("data.chr.",chr),data.ii.chr.i)
  }
}


# mapply(ncol,mget(ls(pattern="data.chr."))) %>% sort

complete.cols.f <- function(df,cent.type){
  filt.v.red <- c("TRUE_FALSE","TRUE_TRUE")
  cent.type.i <- cent.type
  if(ncol(df) == 24 & cent.type.i == "monocentric"){
    print(unique(df[,grepl("chr$",colnames(df))]))
    df$TRUE_FALSE <- NA_character_
    df$TRUE_TRUE <- NA_character_
    return(df)
  }
  if(ncol(df) == 24 & cent.type.i != "monocentric"){
    print(unique(df[,grepl("chr$",colnames(df))]))
    df$TRUE_FALSE <- NA_character_
    df$TRUE_NA <- NA_character_
    return(df)
  }
  if(ncol(df) == 25){
    print(unique(df[,grepl("chr$",colnames(df))]))
    df <- cbind(df, "def" = replicate(nrow(df),NA))
    colnames(df)[26] <- filt.v.red[which(filt.v.red != gsub("data.chr.[a-zA-Z0-9]+.","",colnames(df)[ncol(df)-1]))]
    colnames(df)[26] <- paste0(gsub("(data.chr.[a-zA-Z0-9]+.).*","\\1",colnames(df)[ncol(df)-1]),colnames(df)[26])
    df <- df[,c(1:24,c(25,26)[order(colnames(df)[25:26])])]
    return(df)
  }
  if(ncol(df) == 26){
    print(unique(df[,grepl("chr$",colnames(df))]))
    return(df)
  }
}

df.complete <- data.frame()
for(obj in 1:length(ls(pattern="data.chr."))){
  df.i <- mget(ls(pattern="data.chr.")[obj]) %>%
    data.frame()
  complete.cols.f(df.i,cent.type) -> df.i
  colnames(df.i) <- colnames(df.i) %>%
    {gsub("data.chr.[a-zA-Z0-9]+.","",.)}
  df.complete <- rbind(df.complete,df.i)
}

if(length(grep("intact.gff3",file.ii)) > 0){
  write.table(df.complete, paste0(species.i,"_mod.EDTA.intact.gff3.array_boundaries"),row.names = F, col.names = T, sep = "\t", quote = F)
}

if(length(grep("_edta_filtered.csv.reassigned",file.ii)) > 0){
  write.table(df.complete, paste0(species.i,"_edta_filtered.csv.reassigned.array_boundaries"),row.names = F, col.names = T, sep = "\t", quote = F)
}


# df.complete %>% group_by(chr,TRUE_FALSE) %>% summarise(cent.occ = n()) %>% arrange(desc(cent.occ)) %>% print(n=60)
print(paste0(" \n",species.i,"\n "))
df.complete %>% group_by(chr,TRUE_FALSE) %>% summarise(cent.occ = n()) %>% arrange(desc(cent.occ)) %>% filter(grepl("IN",TRUE_FALSE)) %>% print(n=60)
