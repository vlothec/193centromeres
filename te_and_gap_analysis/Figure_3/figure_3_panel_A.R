
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

gap <- 2
gap.i <- c("TRUE_TRUE","TRUE_FALSE","FALSE_TRUE","FALSE_FALSE")[gap]

non.filter.data <- data.frame()
upper.filter.data <- data.frame()
lower.filter.data <- data.frame()
for(species in 1:length(species.list)){
  species.i <- species.list[species]
  list.files.species <- list.files.v[grepl(species.i,list.files.v)]
  
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
    
    data.i.sum <- data.i.sum %>%
      mutate(prop = round((nt.count/max(nt.count))*100,digits = 5),
             axis.label = ifelse(TE.cls != "centromere", paste0(TE.cls, ": ", str_pad(format(round((nt.count/max(nt.count))*100,digits = 2),nsmall = 2),width = 5,side = "left"),"%"), TE.cls),
             nt.count = nt.count/1000000)
    
    data.i.sum$axis.label <- data.i.sum$axis.label %>%
      {factor(.,levels = rev(unique(data.i.sum$axis.label)[c(9,11,10,12,1:8)]))}
    data.i.sum <- data.i.sum[order(data.i.sum$axis.label),]
    
    data.i.sum <- data.i.sum %>%
      mutate(species = species.i)
    
    non.filter.data <- rbind(non.filter.data,data.i.sum)
    
    # second condition: 0 && upper
    min.gap.length.v <- c(0,50,100,150,200,250)
    min.gap <- 6
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
    
    data.i.sum <- data.i.sum %>%
      mutate(prop = round((nt.count/max(nt.count))*100,digits = 5),
             axis.label = ifelse(TE.cls != "centromere", paste0(TE.cls, ": ", str_pad(format(round((nt.count/max(nt.count))*100,digits = 2),nsmall = 2),width = 5,side = "left"),"%"), TE.cls),
             nt.count = nt.count/1000000)
    
    data.i.sum$axis.label <- data.i.sum$axis.label %>%
      {factor(.,levels = rev(unique(data.i.sum$axis.label)[c(9,11,10,12,1:8)]))}
    data.i.sum <- data.i.sum[order(data.i.sum$axis.label),]
    
    data.i.sum <- data.i.sum %>%
      mutate(species = species.i)
    
    upper.filter.data <- rbind(upper.filter.data,data.i.sum)
    
    # second condition: 0 && upper
    min.gap.length.v <- c(0,50,100,150,200,250)
    min.gap <- 6
    min.gap.length <- min.gap.length.v[min.gap]
    
    side <- 2
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
    
    data.i.sum <- data.i.sum %>%
      mutate(prop = round((nt.count/max(nt.count))*100,digits = 5),
             axis.label = ifelse(TE.cls != "centromere", paste0(TE.cls, ": ", str_pad(format(round((nt.count/max(nt.count))*100,digits = 2),nsmall = 2),width = 5,side = "left"),"%"), TE.cls),
             nt.count = nt.count/1000000)
    
    data.i.sum$axis.label <- data.i.sum$axis.label %>%
      {factor(.,levels = rev(unique(data.i.sum$axis.label)[c(9,11,10,12,1:8)]))}
    data.i.sum <- data.i.sum[order(data.i.sum$axis.label),]
    
    data.i.sum <- data.i.sum %>%
      mutate(species = species.i)
    
    lower.filter.data <- rbind(lower.filter.data,data.i.sum)
    
  }
  rm(list = setdiff(ls(),c("species.list","list.files.v","file.path","TE_general_cls","gap","gap.i","non.filter.data","upper.filter.data","lower.filter.data")))
}



# combined version
non.filter.data <- non.filter.data[with(non.filter.data,order(prop)),]
non.filter.data$species <- non.filter.data$species %>%
  {factor(.,levels = unique(non.filter.data$species[which(non.filter.data$step == 1 & non.filter.data$TE.cls == "TRASH")]))}
non.filter.data <- non.filter.data[with(non.filter.data,order(species)),]

non.filter.data$TE.cls <- non.filter.data$TE.cls %>%
  {factor(.,levels = unique(non.filter.data$TE.cls)[c(10,8,9,1:7,11,12)])}

non.filter.data <- non.filter.data[order(non.filter.data$TE.cls),]


non.filter.data$species <- non.filter.data$species %>%
  {gsub("\\..*","",.)}


# next line can be uncommented in case we want to sort yaxis by taxonomic order
{
  # upload phylo data
  data.attributes <- read.delim("dtol_plant_attributes.txt", header = T, sep = "\t")
  data.attributes <- data.attributes %>%
    mutate(species = gsub("\\..*","",fasta))
  
  data.attributes %>% str

  non.filter.data <- left_join(non.filter.data,data.attributes[,c(9,2,5)],by = "species")
  
  non.filter.data$group.0 <- mapply(function(x) case_when(x == "Dicots" ~ "plant",
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
                                                     .default = NA_character_),non.filter.data$Group)
}


non.filter.data$group.0 <- non.filter.data$group.0 %>%
  {factor(.,levels = c("invertebrate","chordate","plant"))}


# tally
distinct(non.filter.data[,c(10,13)]) %>% {table(.[,2])}

non.filter.data$nt.count <- non.filter.data$nt.count*1000000

# saving y axis labels to apply in the siding plots
# real.y.lables <- distinct(non.filter.data[,c(10,13)]) %>% arrange(group.0)
# write.table(real.y.lables,"real_y_lables.txt",sep = "\t",quote = F,col.names = T,row.names = F)
df.order <- unique(non.filter.data$species)
write.table(df.order,"df_order.txt",sep = "\t",quote = F,col.names = F,row.names = F)

# text.size <- 8
text.size <- 6.5
bar.width <- 0.425
p.1 <- ggplot(data = non.filter.data[c(which(non.filter.data$step == 1 & non.filter.data$TE.cls == "TRASH"),which(non.filter.data$step == 2)),]) +
  geom_col(data = non.filter.data[which(non.filter.data$step == 1 & non.filter.data$TE.cls == "centromere"),],aes(x = prop, y = fct_inorder(species)), color = NA, fill = "#DDDDDD", width = bar.width) +
  geom_col(aes(x = prop, y = fct_inorder(species), fill = TE.cls), color = NA,  width = bar.width) +
  scale_fill_manual(values = c("TRASH" = "#D83067","all_TEs" = "#487CE3", "unknown" = "#DDDDDD")) +
  # scale_x_continuous(labels = ~ ifelse(. < 0, " ",.), expand = expansion(0,0)) +
  scale_y_discrete(expand = expansion(add = .85)) +
  facet_grid(group.0 ~ ., scales = "free_y", space = "free_y") +
  guides(color = "none") +
  coord_cartesian(expand = T) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.key.size = unit(.35,"cm"),
        legend.spacing = unit(2,"cm"),
        legend.text = element_text(size = text.size),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "grey95"),
        plot.margin = unit(c(t = 5,r = 5,b = 5,l = 5),"mm"),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.title = element_blank(),
        axis.text.y = element_text(),
        strip.text.y = element_blank(),
        text = element_text(size = text.size))


pdf(paste0("feature_prop_", gap.i,".all_in_species.non_filter.phylotaxa.pdf"),onefile = T, width = 3.35,height = 7.15)
print(p.1)
dev.off()


scf_name <- mapply(function (Genus,Species) paste0(gsub("(^[A-Za-z]{1}).*","\\1",Genus),". ", Species), data.attributes$Genus, data.attributes$Species)
data.attributes$scf_name <- scf_name

non.filter.data <- left_join(non.filter.data,data.attributes[,c(9,10)],by = "species")

p.1 <- ggplot(data = non.filter.data[c(which(non.filter.data$step == 1 & non.filter.data$TE.cls == "TRASH"),which(non.filter.data$step == 2)),]) +
  geom_col(data = non.filter.data[which(non.filter.data$step == 1 & non.filter.data$TE.cls == "centromere"),],aes(x = prop, y = fct_inorder(scf_name)), color = NA, fill = "#DDDDDD", width = bar.width) +
  geom_col(aes(x = prop, y = fct_inorder(scf_name), fill = TE.cls), color = NA,  width = bar.width) +
  scale_fill_manual(values = c("TRASH" = "#D83067","all_TEs" = "#487CE3", "unknown" = "#DDDDDD")) +
  # scale_x_continuous(labels = ~ ifelse(. < 0, " ",.), expand = expansion(0,0)) +
  scale_y_discrete(expand = expansion(add = .85)) +
  facet_grid(group.0 ~ ., scales = "free_y", space = "free_y") +
  guides(color = "none") +
  coord_cartesian(expand = T) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.key.size = unit(.35,"cm"),
        legend.spacing = unit(2,"cm"),
        legend.text = element_text(size = text.size),
        panel.background = element_rect(fill = "grey97", color = NA),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "grey95"),
        plot.margin = unit(c(t = 5,r = 5,b = 5,l = 5),"mm"),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.title = element_blank(),
        axis.text.y = element_text(face = "italic"),
        strip.text.y = element_blank(),
        text = element_text(size = text.size))

pdf(paste0("feature_prop_", gap.i,".all_in_species.non_filter.phylotaxa.scf_names.pdf"),onefile = T, width = 3.35,height = 7.15)
print(p.1)
dev.off()

