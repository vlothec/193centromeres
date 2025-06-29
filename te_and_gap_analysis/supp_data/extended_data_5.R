
library(cowplot); library(ggnewscale); library(ggpubr); library(grid); library(gridExtra); library(RColorBrewer); library(tidyverse);

file.path <- "data/"
# data from cent_prop_combined.sh
list.files.v <- list.files(file.path, pattern = "_centromere_nt_prop_combined")

list.files.v <- list.files.v[!grepl("TRUE_NA",list.files.v)]
species.list <- list.files.v %>%
  {gsub("\\..*","",.)} %>%
  {unique(.)}


TE_general_cls <- read.table("TE_general_cls.txt", quote="\"", comment.char="")
colnames(TE_general_cls) <- c("cls.0","EDTA.feature")
TE_general_cls <- TE_general_cls[!grepl("repeat_region",TE_general_cls$cls.0),]

TE_general_cls_extended <- read.delim("TE_general_cls_extended", header=FALSE)
colnames(TE_general_cls_extended)[c(3)] <- colnames(TE_general_cls)[c(2)]
colnames(TE_general_cls_extended)[c(2)] <- "TE.cls"
colnames(TE_general_cls_extended)[c(1)] <- "cls.0"
TE_general_cls_extended <- TE_general_cls_extended[c(1,2)] %>%
  distinct()


gap <- 2
gap.i <- c("TRUE_TRUE","TRUE_FALSE","FALSE_TRUE","FALSE_FALSE")[gap]

species_order_list <- read.table("df_order.txt", quote="\"", comment.char="")
species_order_list$V1 <- as.character(species_order_list$V1)
colnames(species_order_list)[c(1)] <- "species"

# taxonomic order
{
data.attributes <- read.delim("dtol_plant_attributes.txt", header = T, sep = "\t")
data.attributes <- data.attributes %>%
  mutate(species = gsub("\\..*","",fasta))

data.attributes %>% str

species_order_list <- left_join(species_order_list,data.attributes[,c(9,2,5)],by = "species")

species_order_list$group.0 <- mapply(function(x) case_when(x == "Dicots" ~ "plant",
                                                        x == "Monocots" ~ "plant",
                                                        x == "Aves" ~ "chordate",
                                                        x == "Mammalia" ~ "chordate",
                                                        x == "Actinopterygii" ~ "chordate",
                                                        x == "Hymenoptera" ~ "invertebrate",
                                                        x == "Annelida" ~ "invertebrate",
                                                        x == "Coleoptera" ~ "invertebrate",
                                                        x == "Diptera" ~ "invertebrate",
                                                        x == "Hemiptera" ~ "invertebrate",
                                                        x == "Cnidaria" ~ "invertebrate",
                                                        x == "Lepidoptera" ~ "invertebrate",
                                                        x == "Neuroptera" ~ "invertebrate",
                                                        x == "Tunicata" ~ "chordate",
                                                        x == "Reptilia" ~ "chordate",
                                                        x == "Bryozoa" ~ "invertebrate",
                                                        x == "Mollusca" ~ "invertebrate",
                                                        x == "Plecoptera" ~ "invertebrate",
                                                        x == "Echinodermata" ~ "invertebrate",
                                                        .default = NA_character_),species_order_list$Group)
}

species_order_list <- species_order_list[with(species_order_list,order(group.0)),]


pdf(paste0("figure_s7.pdf"),onefile = T, width = 30.5,height = 16.5)
par(mfrow = c(8, 10),font = 1, cex.main = 0.75, family = "mono", mai = c(.5,1.835,.5,0.2))
for(species in 1:length(species_order_list$species)){
  species.i <- species_order_list$species[species]
  list.files.species <- list.files.v[grepl(species.i,list.files.v)]
  
  color.l <- c("plant" = "#98C648","chordate" = "#E34084", "invertebrate" = "#3D3CC1")
  col.main.i <- color.l[which(names(color.l) == species_order_list$group.0[species])]
  
  print(species.i)
  
  cent.prop.raw <- read.delim(paste0(file.path,"/",list.files.species), header=FALSE, comment.char = "")
  colnames(cent.prop.raw) <-c("nt.count","filter_edge","min.gap","side","feature","TE.cls","step")
  
  cent.prop.raw <- cent.prop.raw[which(cent.prop.raw$step != 4),]
  
  if(length(grep("TRUE_NA",cent.prop.raw$filter_edge)) == 0){
    
    # first condition: 0 && upper
    min.gap.length.v <- c(0,50,100,150,200,250)
    min.gap <- 1
    min.gap.length <- min.gap.length.v[min.gap]
    
    side <- 1
    side.i <- c("upper","lower")[side]
    
    cent.prop.i <- cent.prop.raw[which(cent.prop.raw$filter_edge == gap.i & cent.prop.raw$side == side.i & cent.prop.raw$min.gap == min.gap.length),]
    
    # panel 3 combined
    data.i.sum <- cent.prop.i[c(which(cent.prop.i$step == 2 & cent.prop.i$feature == "TE_nt"),which(cent.prop.i$step == 3  & cent.prop.i$feature == "TE_nt"),which(cent.prop.i$step == 1)),]
    data.i.sum$TE.cls[which(data.i.sum$step == 2 & data.i.sum$feature == "TE_nt")] <- "all_TEs"
    data.i.sum$TE.cls[which(data.i.sum$step == 1 & data.i.sum$feature == "cent_nt")] <- "centromere"
    data.i.sum$TE.cls[which(data.i.sum$step == 1 & data.i.sum$feature == "gap_nt")] <- "gaps"
    
    # infer TRASH contribution
    data.i.sum.TRASH <- data.frame(nt.count = data.i.sum$nt.count[which(data.i.sum$step == 1 & data.i.sum$feature == "cent_nt")] - data.i.sum$nt.count[which(data.i.sum$step == 1 & data.i.sum$feature == "gap_nt")],
                                   filter_edge = data.i.sum$filter_edge[which(data.i.sum$step == 1 & data.i.sum$feature == "cent_nt")],
                                   min.gap = data.i.sum$min.gap[which(data.i.sum$step == 1 & data.i.sum$feature == "cent_nt")],
                                   side = data.i.sum$side[which(data.i.sum$step == 1 & data.i.sum$feature == "cent_nt")],
                                   feature = "TR_nt",
                                   TE.cls = "TRASH",
                                   step = 1)
    
    data.i.sum <- rbind(data.i.sum,data.i.sum.TRASH)
    
    # infer unaccounted nt
    data.i.sum.unk <- data.frame(nt.count = data.i.sum$nt.count[which(data.i.sum$step == 1 & data.i.sum$feature == "gap_nt")] - data.i.sum$nt.count[which(data.i.sum$step == 2 & data.i.sum$feature == "TE_nt")],
                                 filter_edge = data.i.sum$filter_edge[which(data.i.sum$step == 2 & data.i.sum$feature == "TE_nt")],
                                 min.gap = data.i.sum$min.gap[which(data.i.sum$step == 2 & data.i.sum$feature == "TE_nt")],
                                 side = data.i.sum$side[which(data.i.sum$step == 2 & data.i.sum$feature == "TE_nt")],
                                 feature = "unk_nt",
                                 TE.cls = "unknown",
                                 step = 2)
    
    data.i.sum <- rbind(data.i.sum,data.i.sum.unk)
    
    data.i.sum <- left_join(data.i.sum,TE_general_cls_extended,by= "TE.cls")
    
    TE.cls <- mapply(function(a,b,c) ifelse(is.na(b) == F & c == 3, b, a), data.i.sum$TE.cls,data.i.sum$cls.0,data.i.sum$step)
    data.i.sum$TE.cls <- TE.cls
    
    data.i.sum <- data.i.sum[,-ncol(data.i.sum)]
    
    data.i.sum <- data.i.sum %>%
      mutate(prop = round((nt.count/max(nt.count))*100,digits = 5),
             axis.label = ifelse(TE.cls != "centromere", paste0(TE.cls, ": ", str_pad(format(round((nt.count/max(nt.count))*100,digits = 2),nsmall = 2),width = 5,side = "left"),"%"), TE.cls),
             nt.count = nt.count/1000000)
    
    data.i.sum$axis.label <- data.i.sum$axis.label %>%
      {factor(.,levels = rev(unique(data.i.sum$axis.label)[c(9,11,10,12,1:8)]))}
    data.i.sum <- data.i.sum[order(data.i.sum$axis.label),]
    
    scf_name <- paste0(substr(data.attributes$Genus[which(data.attributes$species == species.i)],1,1),". ",data.attributes$Species[which(data.attributes$species == species.i)])
    if(nrow(data.i.sum) > 0){
      barplot(replicate(length(data.i.sum$nt.count),max(data.i.sum$nt.count)),
              names.arg = data.i.sum$axis.label,
              col="#fbfbfb",
              border = "white",
              horiz = T,
              main = bquote(italic(.(scf_name))),
              col.main = col.main.i,
              las = 1,
              cex.axis = 0.75,
              cex.names = 0.75,
              axisnames = T)
      
      barplot(data.i.sum$nt.count,
              names.arg = data.i.sum$axis.label,
              col=c(replicate(length(unique(data.i.sum$TE.cls[which(data.i.sum$step == 3)])),"orange"), "brown","brown","blue","blue","black"),
              border = "white",
              horiz = T,
              las = 1,
              cex.axis = 0.75,
              cex.names = 0.75,
              axisnames = T,
              add = T)
    }else{
      plot.new()
    }
    
  }
  rm(list = setdiff(ls(),c("species.list","list.files.v","file.path","TE_general_cls","gap","gap.i","species_order_list", "TE_general_cls_extended","data.attributes")))
}
dev.off()
