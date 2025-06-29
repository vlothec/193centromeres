
options(scipen = 999)

library(cowplot); library(ggnewscale); library(ggpubr); library(ggtree); library(ggtreeExtra); library(grid); library(gridExtra); library(RColorBrewer); library(tidyverse);

tree.file <- "plant.Gypsy_Athila.retree2.max100.mafft.fasttree"

tree <- read.tree(file = paste0(tree.file))

tree$tip.label[grepl("_sce",tree$tip.label)]
tree <- treeio::root(tree,outgroup = tree$tip.label[grepl("_sce",tree$tip.label)])

tree$tip.label <- {tree$tip.label %>% gsub("\\&","_",.)}

q <- read.delim("data_for_figure_4_panel_a.tsv",sep = "\t",header = T)

# tree aesthetics
layout.i <- "circular"
layout.i <- "fan"
offset.v <- -0.05
p <- ggtree(tree, layout = layout.i,size=.075, open.angle=7.5)

# get boot strap data
relevant.nodes <- c(7983,10300,10299,10432,10448,12160,10231,13296,10210,10208,14642,9702,9703,9131,9119)

d <- p$data
d <- d[!d$isTip,]
d$label <- as.numeric(d$label)
d$label <- round(d$label * 100)
d <- d[which(d$node %in% relevant.nodes),]

# define colors
col.v <- c("Satellite" = "#EA346F", "Holocentric_satellite" = "#F5BF44", "Transposon" = "#4E87F7", "Mixed" = "#7943E3", "ref" = "black", "unk" = "black")
light.v <- c("Satellite" = "#fad9e4", "Holocentric_satellite" = "#faedcf", "Transposon" = "#d4e1fa", "Mixed" = "#d2befa", "ref" = "grey80", "unk" = "grey80")

species.l <- unique(q$species)
species.l <- species.l[!grepl("unk|ATHILA_exemplars",species.l)]
pdf("tree_circ_by_species.pdf",height = 5.75, width = 5.75,onefile = T)
for(i in 1:length(species.l)){
  
  genome.i <- species.l[i]
  print(genome.i)
  class.i <- unique(q$Class[which(q$genome == genome.i)])
  
  if(length(which(q$genome == genome.i & q$cent.status == "IN")) > 0){
    cent.v <- q$rownames.tag[which(q$genome == genome.i & q$cent.status == "IN")] %>% {paste(.,collapse = "|")}
    
    n.elements <- q$rownames.tag[which(q$genome == genome.i)] %>% length
    
    layout.i <- "circular"
    layout.i <- "fan"
    offset.v <- -0.05
    p <- ggtree(tree, layout = layout.i,size=.075, open.angle=7.5) + 
      geom_tiplab(align = TRUE, aes(subset=grepl(genome.i,label,fixed=FALSE)==TRUE,label = " "),color = light.v[grepl(class.i,names(light.v))], linetype = "solid", linesize = 0.05,offset = .3) +
      geom_tiplab(align = TRUE, aes(subset=grepl(cent.v,label,fixed=FALSE)==TRUE,label = ""),color = col.v[grepl(class.i,names(light.v))], linetype = "solid", linesize = 0.05,offset = .3) +
      geom_label(data=d, aes(label=label), color = "black",size.unit = "pt",size = 3.25,label.padding = unit(0.125, "lines"), label.r = unit(0.00,"lines"),label.size = 0,fill = "grey95", alpha = .75)
  }else{
    layout.i <- "circular"
    layout.i <- "fan"
    offset.v <- -0.05
    p <- ggtree(tree, layout = layout.i,size=.075, open.angle=7.5) + 
      geom_tiplab(align = TRUE, aes(subset=grepl(genome.i,label,fixed=FALSE)==TRUE,label = " "),color = light.v[grepl(class.i,names(light.v))], linetype = "solid", linesize = 0.05,offset = .3) +
      geom_label(data=d, aes(label=label), color = "black",size.unit = "pt",size = 3.25,label.padding = unit(0.125, "lines"), label.r = unit(0.00,"lines"),label.size = 0,fill = "grey95", alpha = .75)
  }
  
  p.1 <- p +
    geom_fruit(data = q,
               geom = "geom_tile",
               offset = 0,
               mapping = aes(y = ID, x = dom.comb, fill = Class),
               pwidth=0.05,
               color = "transparent",
               axis.params=list(),
               grid.params=list(vline = F,
                                color = "transparent")) +
    scale_fill_manual(values = col.v, guide = "none") +
    ggnewscale::new_scale_fill() +
    ggnewscale::new_scale_color()
  
  # pident
  q$target.species <- "targ_0"
  q$target.species[which(q$genome %in% genome.i)] <- "targ_1"
  
  
  q$target.species <- q$target.species %>%
    {factor(.,levels = c("targ_1","targ_0"))}
  
  q <- q[with(q,order(target.species)),]
  
  col.v.i <- c("targ_0" = "grey90","targ_1"= col.v[grepl(class.i,names(col.v))])
  names(col.v.i)[c(2)] <- "targ_1"
  p.2 <- p.1 +
    geom_fruit(data = rbind(q[which(q$target.species == "targ_0"),],
                            q[which(q$target.species == "targ_1"),]),
               geom = "geom_point",
               mapping = aes(y = ID, x = p.iden, color = target.species),
               offset = 0,
               size = .005,
               pwidth = .1,
               axis.params=list(axis = "x",
                                limits = c(min(q$p.iden,na.rm = T),1)),
               grid.params=list(vline = T,
                                color = "grey97")) +
    scale_color_manual(values = col.v.i) +
    ggnewscale::new_scale_color() +
    labs(title = paste0(genome.i,"\n",class.i,"\nn: ", n.elements)) +
    theme(legend.position = "none",
          panel.background = element_rect(colour = NA, fill = NA),
          plot.title = element_text(hjust = .5)) 
  
  print(p.2)
  
  assign(paste0("plot_",i),p.2)
  
  rm(p)
  rm(p.1)
  rm(p.2)
}
dev.off()

