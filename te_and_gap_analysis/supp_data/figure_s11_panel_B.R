

species.v <- c("daLinVulg1.1.fa","daMisOron1.1.fa","drHedHeli1.1.fa","ddEupPepu3.1.fa","daBalNigr1.1.fa","daSheArve1.1.fa")


for(species in 1:length(species.v)){
  species.i <- species.v[species]
  
  # arrange data
  {
    data.ii <- read.delim(paste0(species.i,".TRUE_FALSE.all_cent_gaps.strand.all_gaps.fa.fasta_to_one"), header=FALSE)
    colnames(data.ii) <- c("header","seq","length")
    
    data.ii$gap.coord <- mapply(function(x) gsub(">(.*)#.*","\\1",x), data.ii$header)
    data.ii$GAP.length.id.name <- mapply(function(x) gsub(">.*#(.*)","\\1",x), data.ii$header)
    
    data.ii <- data.ii[with(data.ii,order(length)),]
    
    data.ii$chr <- mapply(function(x) gsub("(.*):[0-9]+\\.\\..*","\\1",x),data.ii$gap.coord)
    data.ii$start <- mapply(function(x) gsub(".*:([0-9]+)\\.\\..*","\\1",x),data.ii$gap.coord) %>% as.numeric()
    data.ii$end <- mapply(function(x) gsub(".*\\.\\.([0-9]+)","\\1",x),data.ii$gap.coord) %>% as.numeric()
    
    
    # upload array boundaries
    data.i <- read.delim(paste0(species.i ,"_all_centromeric_arrays.tsv"), header=FALSE)
    colnames(data.i) <- c("chr","start","end","TRASH.cls","strand","filt_1","filt_2")
    
    data.i <- data.i[which(data.i$filt_1 == "TRUE" & data.i$filt_2 == "FALSE"),]
    data.i$filt_comb <- mapply(function(x,y) paste0(x,"/",y), data.i$filt_1, data.i$filt_2)
    data.i$filt_comb_row <- mapply(function(x,y,z) paste0(x,"/",y,"/",z), data.i$filt_1, data.i$filt_2, rownames(data.i))
    cent.type <- mapply(function(x) ifelse(length(grep("NA",x)) > 0, "holocentrict","monocentric"), data.i$filt_comb_row) %>%
      unique()
    
    data.ii$cent.start <- NA_real_
    data.ii$cent.end <- NA_real_
    data.ii$cent.n <- NA_real_
    for(chr in 1:length(unique(data.i$chr))){
      chr.i <- unique(data.i$chr)[chr]
      cent.start.i <- data.i$start[which(data.i$chr == chr.i)]
      cent.end.i <- data.i$end[which(data.i$chr == chr.i)]
      
      if(length(cent.start.i) >= 1){
        for(arr in 1:length(cent.start.i)){
          data.ii.idx <- which(data.ii$chr == chr.i & data.ii$start >= cent.start.i[arr] & data.ii$end <= cent.end.i[arr])
          data.ii$cent.start[data.ii.idx] <- cent.start.i[arr]
          data.ii$cent.end[data.ii.idx] <- cent.end.i[arr]
          data.ii$cent.n[data.ii.idx] <- arr
        }
      }
    }
    
    scale.range.f <- function(x,y,z){
      (x-y)/(z-y)
    }
    
    data.ii <- data.ii[with(data.ii,order(chr,start)),]
    data.ii <- data.ii %>%
      group_by(chr,cent.n) %>%
      mutate(rescaled.start = scale.range.f(start,cent.start,cent.end),
             rescaled.end = scale.range.f(end,cent.start,cent.end))
    
    data.ii$length <- mapply(function(x,y) y-x, data.ii$start, data.ii$end)
    data.ii$rescaled.length <- mapply(function(x,y) y-x, data.ii$rescaled.start, data.ii$rescaled.end)
    data.ii$chr_cent.n <- mapply(function(x,y) paste0(x,"/",y),data.ii$chr,data.ii$cent.n)
    
    data.ii$cent.n.length <- mapply(function(x,y) x-y,data.ii$cent.end,data.ii$cent.start)
    data.ii$rescaled.bp <- mapply(function(x) 1/x,data.ii$cent.n.length)
    
    scale.df <- data.ii %>%
      select(chr_cent.n,rescaled.bp) %>%
      distinct()
    
  }
  
  
  # bin the data
  bin.v <- c(0.01,0.02)[2]
  for(bin in 1:1){
    binwidth <- bin.v[bin]
    data.ii.sum.all <- data.ii %>%
      mutate(dist.start = findInterval(rescaled.start, seq(0,max(rescaled.end),length.out = (round(max(rescaled.end) / binwidth))))) %>%
      ungroup() %>%
      group_by(dist.start) %>%
      summarise(n.counts = n())
    assign(paste0("all_bin_",bin),data.ii.sum.all)
    
    data.ii.sum.chr <- data.ii %>%
      mutate(dist.start = findInterval(rescaled.start, seq(0,max(rescaled.end),length.out = (round(max(rescaled.end) / binwidth))))) %>%
      ungroup() %>%
      group_by(chr_cent.n,dist.start) %>%
      summarise(n.counts = n())
    assign(paste0("chr_bin_",bin),data.ii.sum.chr)
    
  }
  
  total.gaps <- data.ii$GAP.length.id.name %>% unique() %>% length()
  prop.gaps <- round(total.gaps*0.01) 
  
  print(total.gaps)
  
  # plot
  # simplest version
  bin <- 1
  p.dens <- ggplot(all_bin_1) +
    geom_segment(x = 0, xend = 1, y = prop.gaps, yend = prop.gaps, color = "grey50", linetype = "dotted") +
    geom_col(aes(x = dist.start*bin.v[bin], y = n.counts),color = "black",fill="black") +
    scale_y_continuous(expand = expansion(mult = c(0,.1))) +
    geom_text(data = data.frame("prop.gaps" = prop.gaps,"n.counts" = -0.05),aes(y = prop.gaps, x = n.counts), label = "1%", color = "grey50", size = 10*0.35, hjust=0) +
    labs(title = paste0(species.i,"\n",total.gaps,"\n")) +
    # 0.02 
    coord_cartesian(clip = "on",expand = F, xlim = c(-0.05,1)) +
    
    # 0.01
    # coord_cartesian(clip = "on",expand = F, xlim = c(-0.1,1)) +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 1),
          plot.margin = margin(t = 1, r = 0, b = 0, l = 0),
          axis.title = element_blank(),
          axis.text = element_blank(),
          strip.text = element_text(hjust = 1),
          text = element_text(size = 10.5))
  
  pdf(paste0("mosaic_gap_simple_",species.i,"_new_dist.pdf"), width = 3.5, height = 1.5)
  print(p.dens)
  dev.off()
  
}

