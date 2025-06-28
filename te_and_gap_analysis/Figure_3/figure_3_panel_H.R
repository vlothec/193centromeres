

library(cowplot); library(ggnewscale); library(ggpubr); library(grid); library(gridExtra); library(RColorBrewer); library(tidyverse);

file.path <- "data/"
# data from gap_infer_prop_all_gaps.R
list.files.v <- list.files(file.path, pattern = "all.gaps.parsed$")

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

data.attributes <- read.delim("dtol_plant_attributes.txt", header = T, sep = "\t")
data.attributes %>% str
scf_name <- mapply(function (Genus,Species) paste0(gsub("(^[A-Za-z]{1}).*","\\1",Genus),". ", Species), data.attributes$Genus, data.attributes$Species)
data.attributes$scf_name <- scf_name
data.attributes$species <- data.attributes$fasta %>%
  {gsub("\\..*","",.)}
data.attributes <- data.attributes[,c(10,9)]

min.gap.length.v <- c(0,50,100,150,200,250)
min.gap <- 6
min.gap.length <- min.gap.length.v[min.gap]
gap <- 2

species.list <- c("drHedHeli1","ddEupPepu3","daBalNigr1","daSheArve1","daMisOron1","daLinVulg1")

plt_ <- TRUE
gap.v <- vector()
p.1 <- list()
p.2 <- list()
p.3 <- list()
p.4 <- list()
p.5 <- list()

for(species in 1:length(species.list)){
  
  species.i <- species.list[species]
  print(species.i)
  list.files.species <- list.files.v[grepl(species.i,list.files.v)]
  
  gap.i <- c("TRUE_TRUE","TRUE_FALSE","FALSE_TRUE","FALSE_FALSE")[gap]
  
  if(length(grep(gap.i,list.files.species)) > 0){
    
    file.i <- list.files.species[grepl(gap.i,list.files.species)]
    data.i <- read.delim(paste0(file.path,"/",file.i), header=TRUE, comment.char = "")
    
    data.i <- data.i[which(data.i$GAP.length >= 0),]
    data.i <- data.i %>%
      mutate(EDTA.start = ifelse(EDTA.start == -1, abs(EDTA.start*as.numeric(rownames(.))),EDTA.start),
             EDTA.end = ifelse(EDTA.end == -1, abs(EDTA.end*as.numeric(rownames(.))),EDTA.end))
    
    GAP.type.i <- unique(data.i$GAP.type)
    gaps.n <- length(unique(data.i$GAP.length.id))
    
    if(nrow(data.i) > 0){
      data.i <- left_join(data.i,TE_general_cls,by = "EDTA.feature")
    }
    
    data.i <- data.i %>%
      mutate(gap.rel.start = TRASH.start - TRASH.start,
             gap.rel.end = TRASH.end - TRASH.start,
             EDTA.rel.start = EDTA.start - TRASH.start,
             EDTA.rel.end = EDTA.end - TRASH.start) %>%
      mutate(EDTA.rel.start = ifelse(EDTA.rel.start < 0, 0, EDTA.rel.start),
             EDTA.rel.end = ifelse(EDTA.rel.end < 0, 0, EDTA.rel.end),
             EDTA.rel.end = ifelse(EDTA.rel.end > GAP.length, GAP.length, EDTA.rel.end)) %>%
      group_by(GAP.length.id) %>%
      mutate(GAP.architecture = NA)
    
    data.i.relevant <- data.i[which(data.i$EDTA.start != data.i$EDTA.end),]
    data.i.relevant <- data.i.relevant[order(data.i.relevant$GAP.length,decreasing = T),]
    
    if(length(species.i) == 1){
      
      intact_list <- try(read.table(paste0(species.i,"_intact_list_200_bp.txt"), quote="\"", comment.char="",sep = " "))
      
      if(length(intact_list) > 0){
        colnames(intact_list)[c(1)] <- "GAP.length.id"
        intact_list$target.status <- "intact"
      }
      
      solo_list <- try(read.table(paste0(species.i,"_solo_list_100_bp.txt"), quote="\"", comment.char="",sep = " "))
      
      if(class(solo_list) != "try-error"){
        colnames(solo_list)[c(1)] <- "GAP.length.id"
        solo_list$target.status <- "solo"
        solo_list <- solo_list[which(solo_list$V4 > 150),] 
      }
      
      intact_solo_list <- rbind(intact_list,solo_list)
      
      data.i.relevant <- left_join(data.i.relevant,distinct(intact_solo_list),by = "GAP.length.id")
    }else{
      data.i.relevant$target.status <- NA
    }
    
    
    data.i.relevant$GAP.length.id <- data.i.relevant$GAP.length.id %>%
      {factor(.,levels = rev(unique(data.i.relevant$GAP.length.id)))} %>%
      as.numeric()
    
    data.i.relevant <- data.i.relevant[with(data.i.relevant,order(GAP.length.id)),]
    

    # axis limits
    if (length(unique(data.i.relevant$GAP.length.id)) > 636) {
      ylim <- c(0, length(unique(data.i.relevant$GAP.length.id)))
    } else {
      ylim <- c(-(636-length(unique(data.i.relevant$GAP.length.id))), length(unique(data.i.relevant$GAP.length.id)))
    }
    xlim <- c(-1,max(data.i.relevant$GAP.length)/1000)
    
    # all in
    GAP.architecture.v <- c("mosaic","mult_TE","single_TE")
    data.i.relevant$cls.0 <- data.i.relevant$cls.0 %>%
      {factor(.,levels = unique(TE_general_cls$cls.0))}
    
    data.i.relevant$species <- species.i
    data.i.relevant <- left_join(data.i.relevant,data.attributes,by = "species")
    
    # data.i.relevant %>% str
    
    gap.n.i <- data.i.relevant$GAP.length.id %>% unique() %>% length()
    gap.v <- c(gap.v,gap.n.i)
    
    
    legend.label <- TE_general_cls_extended$cls.0
    names(legend.label) <- TE_general_cls_extended$TE.cls
    
    if(plt_ == T){
    
    p.1[[length(p.1)+1]] <- ggplot(data.i.relevant) +
      geom_text(data = data.i.relevant[which(data.i.relevant$target.status == "intact"),],aes(x = -0.85, y = GAP.length.id,label = "-"), color = "black", size = 3, vjust = 0.34) +
      geom_text(data = data.i.relevant[which(data.i.relevant$target.status == "solo"),],aes(x = -0.85, y = GAP.length.id,label = "-"), color = "grey80", size = 3, vjust = 0.34) +
      geom_segment(aes(x = 0, xend = GAP.length/1000, y = GAP.length.id, yend = GAP.length.id), color = "grey90", lineend = "round") +
      geom_segment(aes(x = EDTA.rel.start/1000, xend = EDTA.rel.end/1000, y = GAP.length.id, yend = GAP.length.id, color = cls.0)) +
      scale_color_manual(values = c("LTR_retrotransposon" = "darkblue", "Class_I_Retrotransposon" = "orange", "nonLTR_retrotransposon" = "brown", "Class_II_DNA_Transposon" = "#F8D36B", "repeat_region" = "red", "rRNA_gene" = "cyan", "TE_unclass" = "black"), limits = c("LTR_retrotransposon", "Class_I_Retrotransposon", "nonLTR_retrotransposon", "Class_II_DNA_Transposon", "repeat_region", "rRNA_gene", "TE_unclass"),labels = legend.label) +
      scale_y_continuous(expand = expansion(mult = c(0.015))) +
      scale_x_continuous(position = "top") +
      labs(title = paste0(" \n",unique(data.i.relevant$scf_name)," (n: ",gap.n.i,")\n ")) +
      guides(color = guide_legend(byrow = TRUE)) +
      coord_cartesian(xlim = c(-1,100), ylim = ylim) +
      theme_transparent() +
      theme(legend.position = "none",
            legend.direction = "vertical",
            legend.title = element_blank(),
            legend.key.size = unit(0.15, "cm"),
            legend.spacing.y = unit(0.15,"cm"),
            legend.text = element_text(size = 4),
            plot.title = element_text(hjust = 0.5, vjust = 1,face = "italic", margin = margin(t = 0, b = 0, l = 0, r = 0, unit = "pt")),
            panel.background = element_rect(fill = NA, colour = NA),
            panel.grid.major.x = element_line(colour = "grey99",linetype = "dotted"),
            axis.line.x.top = element_line(),
            axis.ticks.x.top = element_line(),
            axis.text.x.top = element_text(),
            strip.text.x = element_text(),
            strip.background = element_rect(colour = NA, fill = NA),
            strip.placement = "outside",
            text = element_text(size = 12))
    
    p.2[[length(p.2)+1]] <- ggplot(data.i.relevant) +
      geom_text(data = data.i.relevant[which(data.i.relevant$target.status == "intact"),],aes(x = -0.85, y = GAP.length.id,label = "-"), color = "black", size = 3, vjust = 0.34) +
      geom_text(data = data.i.relevant[which(data.i.relevant$target.status == "solo"),],aes(x = -0.85, y = GAP.length.id,label = "-"), color = "grey80", size = 3, vjust = 0.34) +
      geom_segment(aes(x = 0, xend = GAP.length/1000, y = GAP.length.id, yend = GAP.length.id), color = "grey90", lineend = "round") +
      geom_segment(aes(x = EDTA.rel.start/1000, xend = EDTA.rel.end/1000, y = GAP.length.id, yend = GAP.length.id, color = cls.0)) +
      scale_colour_manual(values = c("LTR_retrotransposon" = "darkblue", "Class_I_Retrotransposon" = "orange", "nonLTR_retrotransposon" = "brown", "Class_II_DNA_Transposon" = "#F8D36B", "repeat_region" = "red", "rRNA_gene" = "cyan", "TE_unclass" = "black"), limits = c("LTR_retrotransposon", "Class_I_Retrotransposon", "nonLTR_retrotransposon", "Class_II_DNA_Transposon", "repeat_region", "rRNA_gene", "TE_unclass"),labels = legend.label) +
      scale_y_continuous(expand = expansion(mult = c(0.015))) +
      scale_x_continuous(position = "top") +
      labs(title = paste0(" \n",unique(data.i.relevant$scf_name)," (n: ",gap.n.i,")\n ")) +
      guides(color = guide_legend(byrow = TRUE)) +
      coord_cartesian(xlim = c(-0.1,20), ylim = ylim) +
      theme_transparent() +
      theme(legend.position = "none",
            legend.direction = "vertical",
            legend.title = element_blank(),
            legend.key.size = unit(0.15, "cm"),
            legend.spacing.y = unit(0.15,"cm"),
            legend.text = element_text(size = 4),
            plot.title = element_text(hjust = 0.5, vjust = 1,face = "italic", margin = margin(t = 0, b = 0, l = 0, r = 0, unit = "pt")),
            panel.background = element_rect(fill = NA, colour = NA),
            panel.grid.major.x = element_line(colour = "grey99",linetype = "dotted"),
            axis.line.x.top = element_line(),
            axis.ticks.x.top = element_line(),
            axis.text.x.top = element_text(),
            strip.text.x = element_text(),
            strip.background = element_rect(colour = NA, fill = NA),
            strip.placement = "outside",
            text = element_text(size = 12))
    
    p.5[[length(p.5)+1]] <- ggplot(data.i.relevant) +
      geom_text(data = data.i.relevant[which(data.i.relevant$target.status == "intact"),],aes(x = -0.85, y = GAP.length.id,label = "-"), color = "black", size = 3, vjust = 0.34) +
      geom_text(data = data.i.relevant[which(data.i.relevant$target.status == "solo"),],aes(x = -0.85, y = GAP.length.id,label = "-"), color = "grey80", size = 3, vjust = 0.34) +
      geom_segment(aes(x = 0, xend = GAP.length/1000, y = GAP.length.id, yend = GAP.length.id), color = "grey90", lineend = "round") +
      geom_segment(aes(x = EDTA.rel.start/1000, xend = EDTA.rel.end/1000, y = GAP.length.id, yend = GAP.length.id, color = cls.0)) +
      scale_colour_manual(values = c("LTR_retrotransposon" = "darkblue", "Class_I_Retrotransposon" = "orange", "nonLTR_retrotransposon" = "brown", "Class_II_DNA_Transposon" = "#F8D36B", "repeat_region" = "red", "rRNA_gene" = "cyan", "TE_unclass" = "black"), limits = c("LTR_retrotransposon", "Class_I_Retrotransposon", "nonLTR_retrotransposon", "Class_II_DNA_Transposon", "repeat_region", "rRNA_gene", "TE_unclass"),labels = legend.label) +
      scale_y_continuous(expand = expansion(mult = c(0.015))) +
      scale_x_continuous(position = "top") +
      labs(title = paste0(" \n",unique(data.i.relevant$scf_name)," (n: ",gap.n.i,")\n ")) +
      guides(color = guide_legend(byrow = TRUE)) +
      coord_cartesian(xlim = c(-0.1,25), ylim = ylim) +
      theme_transparent() +
      theme(legend.position = "none",
            legend.direction = "vertical",
            legend.title = element_blank(),
            legend.key.size = unit(0.15, "cm"),
            legend.spacing.y = unit(0.15,"cm"),
            legend.text = element_text(size = 4),
            plot.title = element_text(hjust = 0.5, vjust = 1,face = "italic", margin = margin(t = 0, b = 0, l = 0, r = 0, unit = "pt")),
            panel.background = element_rect(fill = NA, colour = NA),
            panel.grid.major.x = element_line(colour = "grey99",linetype = "dotted"),
            axis.line.x.top = element_line(),
            axis.ticks.x.top = element_line(),
            axis.text.x.top = element_text(),
            strip.text.x = element_text(),
            strip.background = element_rect(colour = NA, fill = NA),
            strip.placement = "outside",
            text = element_text(size = 12))
    
    }

  }
  
  
  rm(list = setdiff(ls(),c("file.path","list.files.v","list.files.v2","species.list","TE_general_cls","min.gap.length.v","min.gap.length","gap","gap.length.v","gap.v","plt_","p.1","p.2","data.attributes","TE_general_cls_extended","p.5")))

  }


pdf(paste0(file.path,"gap_coverage_EDTA_cls_gap_color_target_species_red_",length(species.list),".pdf"),onefile = T, width = (length(species.list)*18)/7,height = 9)
cowplot::plot_grid(plotlist=p.1,ncol = length(species.list),align = "hv", axis = "l")
cowplot::plot_grid(plotlist=p.2,ncol = length(species.list),align = "hv", axis = "l")
cowplot::plot_grid(plotlist=p.5,ncol = length(species.list),align = "hv", axis = "l")
dev.off()

