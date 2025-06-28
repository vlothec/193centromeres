
library(cowplot); library(ggnewscale); library(ggpubr); library(grid); library(gridExtra); library(RColorBrewer); library(tidyverse);

options(scipen=999)

# proportions

intact.prop <- 0.290919
frag.prop <- 0.584015
prot.prop <- 0.02194767 # PROT found in GAPS
non_athila_LTRRT.prop <- 0.0358482
unclassified_LTRRT.prop <- 0.0163805
nonLTR.prop <- 0.0123491

cent_feature_prop <- data.frame("chr" = "all_CENT",
                                           "ft.comb.red.ii" = c("other TEs","GETHILA"),
                                           "total.cent.feature.space" = c(0.0646,(intact.prop+frag.prop+prot.prop)))

text.size <- 12

pdf("cent_feature_prop_barplot_all_anno.pdf", width = 9, height = 1.15, onefile = T)
ggplot(cent_feature_prop) +
  # grid
  # ATHILA intact
  geom_vline(xintercept = intact.prop, linetype = "solid",color = "white",linewidth = 1.25) +
  
  # all ATHILA
  geom_vline(xintercept = (intact.prop+frag.prop+prot.prop), linetype = "solid",color = "white",linewidth = 1.25) +
  
  # any TE
  geom_vline(xintercept = (intact.prop+frag.prop+non_athila_LTRRT.prop+unclassified_LTRRT.prop+nonLTR.prop+prot.prop), linetype = "solid",color = "white",linewidth = 1.25) +
  
  # bars
  geom_col(data = data.frame("chr" = "all_CENT", "total.cent.feature.space" = 1),aes(y = chr, x = total.cent.feature.space), fill = "grey92", color = "black", width = 0.2) +
  geom_col(aes(y = chr, x = total.cent.feature.space, fill = fct_inorder(ft.comb.red.ii)), color = "black", width = 0.2,show.legend=TRUE) +
  
  # ATHILA intact
  geom_col(data = data.frame("chr" = "all_CENT", "total.cent.feature.space" = intact.prop),aes(y = chr, x = total.cent.feature.space), fill = NA, color = "black", width = 0.2) +

  # labels
  # ATHILA intact label
  geom_text(x = intact.prop/2, y = "all_CENT", label = "INTACT",size.unit = "pt",size =10) +
  geom_text(x = intact.prop/2, y = 1.2, label = "29.1%",size.unit = "pt",size =10) +

  # ATHILA fragment label
  geom_text(x = intact.prop+(frag.prop/2), y = "all_CENT", label = "FRAGMENT",size.unit = "pt",size =10) +
  geom_text(x = intact.prop+(frag.prop/2), y = 1.2, label = "60.7%",size.unit = "pt",size =10) +

  # other TEs label
  geom_text(x = intact.prop+frag.prop+prot.prop+((non_athila_LTRRT.prop+unclassified_LTRRT.prop+nonLTR.prop)/2), y = 1.2, label = "6.5%",size.unit = "pt",size =10) +
  
  # unk label
  geom_text(x = ((intact.prop+frag.prop+prot.prop)+(non_athila_LTRRT.prop+unclassified_LTRRT.prop+nonLTR.prop))+((1-((intact.prop+frag.prop+prot.prop)+(non_athila_LTRRT.prop+unclassified_LTRRT.prop+nonLTR.prop)))/2), y = 1.2, label = "3.7%",size.unit = "pt",size =10) +
  
  scale_x_continuous(expand = expansion(0), labels = c(0,25,50,75,100)) +
  scale_y_discrete(expand = expansion(add = c(.3, .3))) +
  scale_fill_manual(values = rev(c("GETHILA" = "orange",
                                   "other TEs" = "grey50",
                                   "unknown" = "grey92")),drop = FALSE) +
  coord_cartesian(expand = T,clip = "off") +
  theme_minimal() + 
  guides(
    fill = guide_legend(override.aes = list(linewidth = .5,byrow = TRUE))
  ) +
  theme(axis.title = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_blank(),
        legend.position = "right",
        legend.key = element_rect(),
        legend.key.size = unit(.85,"line"),
        legend.text = element_text(),
        legend.title = element_blank(),
        legend.spacing.y = unit(5,"line"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "grey97",colour = NA),
        text = element_text(size = text.size)) 

dev.off()

