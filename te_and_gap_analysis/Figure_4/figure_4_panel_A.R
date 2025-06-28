

options(scipen = 999)

library(cowplot); library(ggnewscale); library(ggpubr); library(ggtree); library(ggtreeExtra); library(grid); library(gridExtra); library(RColorBrewer); library(tidyverse);

tree.file <- "plant.Gypsy_Athila.retree2.max100.mafft.fasttree"

tree <- read.tree(file = paste0(tree.file))
tree$tip.label[grepl("_sce",tree$tip.label)]
tree <- treeio::root(tree,outgroup = tree$tip.label[grepl("_sce",tree$tip.label)])

tree$tip.label <- {tree$tip.label %>% gsub("\\&","_",.)}

q <- read.delim("data_for_figure_4_panel_a.tsv",sep = "\t",header = T)

q$genome[which(q$genome == "drMalSylv7")] <- "drMalDome5"
cent.v <- q$rownames.tag[which(q$genome == "drGeuUrba1" & q$cent.status == "IN")] %>% {paste(.,collapse = "|")}
cent.v.mal <- q$rownames.tag[which(q$genome == "drMalDome5" & q$cent.status == "IN")] %>% {paste(.,collapse = "|")}
cent.v.ros <- q$rownames.tag[which(q$genome == "rosCan_S27_v1" & q$cent.status == "IN")] %>% {paste(.,collapse = "|")}
cent.v.ara <- q$rownames.tag[which(q$genome == "ddAraThal4" & q$cent.status == "IN")] %>% {paste(.,collapse = "|")}

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

p <- ggtree(tree, layout = layout.i,size=.075, open.angle=7.5) + 
  geom_tiplab(align = TRUE, aes(label = " "),color = "grey95", linetype = "solid", linesize = 0.05,offset = offset.v) +
  geom_tiplab(align = TRUE, aes(subset=grepl("ddAraThal4",label,fixed=FALSE)==TRUE,label = " "),color = light.v[1], linetype = "solid", linesize = 0.05,offset = .3) +
  geom_tiplab(align = TRUE, aes(subset=grepl("drGeuUrba",label,fixed=FALSE)==TRUE,label = " "),color = "#c4d8ff", linetype = "solid", linesize = 0.05,offset = .4) +
  geom_tiplab(align = TRUE, aes(subset=grepl("drMalDome5",label,fixed=FALSE)==TRUE,label = " "),color = "#e3ebfa", linetype = "solid", linesize = 0.05,offset = .6) +
  geom_tiplab(align = TRUE, aes(subset=grepl("rosCan_S27_v1",label,fixed=FALSE)==TRUE,label = " "),color = light.v[4], linetype = "solid", linesize = 0.05,offset = .5) +
  geom_tiplab(align = TRUE, aes(subset=grepl(cent.v.ara,label,fixed=FALSE)==TRUE,label = ""),color = col.v[1], linetype = "solid", linesize = 0.05,offset = .3) + 
  geom_tiplab(align = TRUE, aes(subset=grepl(cent.v,label,fixed=FALSE)==TRUE,label = ""),color = col.v[3], linetype = "solid", linesize = 0.05,offset = .4) + 
  geom_tiplab(align = TRUE, aes(subset=grepl(cent.v.mal,label,fixed=FALSE)==TRUE,label = ""),color = col.v[3], linetype = "solid", linesize = 0.05,offset = .6)+  
  geom_tiplab(align = TRUE, aes(subset=grepl(cent.v.ros,label,fixed=FALSE)==TRUE,label = ""),color = col.v[4], linetype = "solid", linesize = 0.05,offset = .5)+
  geom_label(data=d, aes(label=label), color = "black",size.unit = "pt",size = 3.25,label.padding = unit(0.125, "lines"), label.r = unit(0.00,"lines"),label.size = 0,fill = "grey95", alpha = .75)


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
q$target.species[which(q$genome %in% c("drGeuUrba1","drMalDome5"))] <- "targ_1"
q$target.species[which(q$genome %in% c("rosCan_S27_v1"))] <- "targ_2"
q$target.species[which(q$genome %in% c("ddAraThal4"))] <- "targ_3"


q$target.species <- q$target.species %>%
  {factor(.,levels = c("targ_3","targ_2","targ_1","targ_0"))}

q <- q[with(q,order(target.species)),]

p.2 <- p.1 +
  geom_fruit(data = rbind(q[which(q$target.species == "targ_0"),],
                          q[which(q$target.species == "targ_1"),],
                          q[which(q$target.species == "targ_2"),],
                          q[which(q$target.species == "targ_3"),]),
             geom = "geom_point",
             mapping = aes(y = ID, x = p.iden, color = target.species),
             offset = 0,
             size = .005,
             pwidth = .1,
             axis.params=list(axis = "x",
                              limits = c(min(q$p.iden,na.rm = T),1),
                              line.color == "transparent",
                              nbreak = 2,
                              text.size = 0),
             grid.params=list(vline = T,
                              color = "grey97")) +
  scale_color_manual(values = c("targ_0" = "grey90","targ_1"="#4E87F7","targ_2"="#7943E3","targ_3"="#EA346F")) +
  ggnewscale::new_scale_color() +
  theme(legend.position = "none",
        panel.background = element_rect(colour = NA, fill = NA))

pdf("figure_4_panel_A.pdf",height = 5.75, width = 5.75)
p.2
dev.off()
