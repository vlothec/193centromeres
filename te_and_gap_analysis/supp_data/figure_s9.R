
library(cowplot); library(ggnewscale); library(ggpubr); library(grid); library(gridExtra); library(RColorBrewer); library(tidyverse);

# data from sat_invasion_tally.R
data.i <- read.delim("all.mod.EDTA.intact.fa.rexdb-plant.cls.tsv.tally",header = F)
colnames(data.i) <- c("Clade","cent.occ","n.count","nt.sum","species")

data.i %>% str
# remove mixture elements
data.i <- data.i[-which(data.i$Clade == "mixture"),]

data.i <- data.i %>%
  group_by(species,cent.occ) %>%
  mutate(nt.sum.prop = nt.sum/sum(nt.sum),
         n.count.prop = n.count/sum(n.count),
         com.label = paste0(species," - (n: ",sum(n.count),")"))

copia_lineages <- read.table("copia_lineages.txt", quote="\"", comment.char="")
copia_lineages$cls <- "copia"
ty3_lineages <- read.table("ty3_lineages.txt", quote="\"", comment.char="")
ty3_lineages$cls <- "gypsy"
lineages.l <- rbind(copia_lineages,ty3_lineages)
colnames(lineages.l)[c(1)] <- "Clade"

data.i <- left_join(data.i,lineages.l,by = "Clade")

data.i$cls <- data.i$cls %>%
  {factor(.,levels = c("copia","gypsy"))}

dtol_taxonomic_order <- data.frame("species" = c("lpJunEffu1","dcPolAvic1","drHedHeli1","daSheArve1","daLinVulg1","daMisOron1","daScuGale1","daBalNigr1","drChaAngu1","ddAraThal4","ddMerAnnu1","ddEupPepu3","drMedArab1","dhQueRobu3","drParJuda1","drPotAnse1"),
                                   "LP" = seq(1,length(unique(data.i$species)),1))

dtol_taxonomic_order <- dtol_taxonomic_order[with(dtol_taxonomic_order,order(LP)),]
dtol_taxonomic_order <- dtol_taxonomic_order[which(dtol_taxonomic_order$species %in% unique(data.i$species)),]

species.list <- dtol_taxonomic_order$species

data.i$Clade <- data.i$Clade %>%
  {factor(.,levels = lineages.l$Clade)}

data.i <- data.i[with(data.i,order(Clade,cls)),]

text.size <- 12
p.1 <- list()
p.2 <- list()
for(sp in 1:length(species.list)){
  sp.i <- species.list[sp]
  linewidth <- 0.5
  
  if(sp == 1){
    font.col <- "black"
  }else{
    font.col <- "transparent"
  }
  

  p.1[[length(p.1) + 1]] <- ggplot(data.i[which(data.i$species == sp.i & data.i$cent.occ == "IN"),],) +
    geom_col(data = data.frame("species" = sp.i, "Clade" = unique(data.i$Clade[which(data.i$cent.occ == "IN")])),aes(x = species, y = 1), fill = "grey95",color = "transparent", width=4.5) +
    geom_col(data = data.frame("species" = sp.i, "Clade" = unique(data.i$Clade[which(data.i$species == sp.i & data.i$cent.occ == "IN")])),aes(x = species, y = 1), fill = "transparent",color = "black", width=4.5) +
    geom_col(aes(x = species, y = nt.sum.prop,fill = cls),color = "black", width=4.5) +
    geom_text(aes(x = species, y = 0, label = paste0("n: ",n.count,"\n",round(nt.sum/1000,digits = 1)," kb")), size.unit = "pt", size = 5.15) +
    scale_color_manual(values = c("transparent","black")) +
    scale_fill_manual(values = c("copia" = "#D83067","gypsy" = "#487CE3")) +
    coord_polar("y", start=0) +
    facet_grid(com.label~Clade,labeller = labeller(.rows = ~ gsub(" - ","\n",.))) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.line.y = element_line(colour = NA, linewidth = linewidth),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.title = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = NA, colour = NA, linewidth = linewidth*2),
          panel.ontop = F,
          plot.margin = margin(t=0,r=0,b=0,l=0),
          plot.title = element_text(hjust = 0),
          strip.text.y = element_text(angle = 0,hjust = 1,size = 9),
          strip.text.x = element_text(angle = 0,hjust = .5,vjust = .5,size = 9,color = font.col),
          text = element_text(size = text.size))
  
  p.2[[length(p.2) + 1]] <- ggplot(data.i[which(data.i$species == sp.i & data.i$cent.occ == "IN"),],) +
    geom_col(data = data.frame("species" = sp.i, "Clade" = unique(data.i$Clade[which(data.i$cent.occ == "IN")])),aes(x = species, y = 1), fill = "grey95",color = "transparent", width=4.5) +
    geom_col(data = data.frame("species" = sp.i, "Clade" = unique(data.i$Clade[which(data.i$species == sp.i & data.i$cent.occ == "IN")])),aes(x = species, y = 1), fill = "transparent",color = "black", width=4.5) +
    geom_col(aes(x = species, y = n.count.prop,fill = cls),color = "black", width=4.5) +
    geom_text(aes(x = species, y = 0, label = paste0("n: ",n.count,"\n",round(nt.sum/1000,digits = 1)," kb")), size.unit = "pt", size = 5.15) +
    scale_color_manual(values = c("transparent","black")) +
    scale_fill_manual(values = c("copia" = "#D83067","gypsy" = "#487CE3")) +
    coord_polar("y", start=0) +
    facet_grid(com.label~Clade,labeller = labeller(.rows = ~ gsub(" - ","\n",.))) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.line.y = element_line(colour = NA, linewidth = linewidth),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.title = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = NA, colour = NA, linewidth = linewidth*2),
          panel.ontop = F,
          plot.margin = margin(t=0,r=0,b=0,l=0),
          plot.title = element_text(hjust = 0),
          strip.text.y = element_text(angle = 0,hjust = 1,size = 9),
          strip.text.x = element_text(angle = 0,hjust = .5,vjust = .5,size = 9,color = font.col),
          text = element_text(size = text.size))
  
}


pdf("figure_s9.pdf",height = 15,width = 15,onefile = T)
print(do.call("grid.arrange", c(p.2, ncol=1, nrow=length(unique(data.i$species)))))
print(do.call("grid.arrange", c(p.1, ncol=1, nrow=length(unique(data.i$species)))))
dev.off()
