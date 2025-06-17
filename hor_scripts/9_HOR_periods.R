setwd("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs/self_HOR_out_loose_settings") 

library(data.table)

ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}

hor_dirs <- list.dirs(path = ".", recursive = FALSE)

summary_df <- data.frame(assembly = hor_dirs, 
                         repeat_name = rep("", length(hor_dirs)),
                         repeat_count = rep("", length(hor_dirs)),
                         repeat_mean_width = rep("", length(hor_dirs)),
                         repeat_total_bp = rep("", length(hor_dirs)),
                         chr_no_analysed = rep("", length(hor_dirs)),
                         chr_names_analysed = rep("", length(hor_dirs)),
                         mean_repetitiveness = rep("", length(hor_dirs)))

for(i in seq_along(hor_dirs)) {
  print(i)
  summary_df$assembly[i] <- strsplit(summary_df$assembly[i], split = "/")[[1]][2]
  summary_df$assembly[i] <- strsplit(summary_df$assembly[i], split = ".fa")[[1]][1]
}

rep_dir_a <- "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs/"


max_hor_distance <- 60


taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 193
print(i)

# for every genome
# for(i in seq_along(hor_dirs)) 
{
  
  # if(grepl("GRCh38", summary_df$assembly[i])) {
  #   repeats <- read.csv(file = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs/GCA_000001405.29_GRCh38.p14_genomic.fa/GCA_000001405.29_GRCh38.p14_genomic.fa_repeats_with_seq.csv")
  #   
  # } else {
  #   repeats <- read.csv(file = paste0(rep_dir_a, summary_df$assembly[i], ".fa/", summary_df$assembly[i], "_repeats_filtered.csv"))
  #   
  # }
  
  cat(i, summary_df$assembly[i], "\n")
  
  hor_files <- list.files(path = hor_dirs[i])
  if(length(hor_files) == 0) {
    cat("No HOR files found\n")
    next
  }
  
  all_hors_files <- hor_files[grep("HORs_", hor_files)]
  all_hors_files <- all_hors_files[grep("csv", all_hors_files)]
  repeats_with_hors_files <- hor_files[grep("repeats_with_hors_", hor_files)]
  repeats_with_hors_files <- repeats_with_hors_files[!grepl("HOR_scored_repeats", repeats_with_hors_files)]
  if(length(repeats_with_hors_files) == 0) {
    cat("No repeats with HOR files found\n")
    next
  }
  if(length(all_hors_files) == 0) {
    cat("No HOR files found\n")
    next
  }
  
  for(j in seq_along(repeats_with_hors_files)) {
    repeat_class <- paste(strsplit(repeats_with_hors_files[j], split = "_")[[1]][4:5], collapse = "_")
    chromosome_name <- strsplit(repeats_with_hors_files[j], split = paste0(repeat_class, "_"))[[1]][2]
    chromosome_name <- strsplit(chromosome_name, split = ".csv")[[1]][1]
    cat(i, "/", length(hor_dirs), summary_df$assembly[i], j, "/", length(repeats_with_hors_files), 
        chromosome_name, repeat_class, "")
    
    
    # if(file.size(paste0(hor_dirs[i], "/HORs_", repeat_class, "_", chromosome_name, ".csv"))/1048576 > 5000) {
    #   cat(" file bigger than 5000 MB, will process later \n")
    #   next
    # }
    
    repeats_chr <- read.csv(file = paste0(hor_dirs[i], "/", repeats_with_hors_files[j]))
    
    if(file.exists(file.path("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs/HOR_periods",
                             paste0("new_HORs_lines_simple_", summary_df$assembly[i], "-", repeat_class, "-", chromosome_name, "-repeat_number_", nrow(repeats_chr), ".png", collapse = "")))) {
      next
    }
    
    if(!file.exists(paste0(hor_dirs[i], "/HORs_", repeat_class, "_", chromosome_name, ".csv"))) next
    hors <- read.csv(file = paste0(hor_dirs[i], "/HORs_", repeat_class, "_", chromosome_name, ".csv"))
    
    cat("plotting...")
    
    png(filename = file.path("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs/HOR_periods", 
                             paste0("new_HORs_lines_simple_", summary_df$assembly[i], "-", repeat_class, "-", chromosome_name, "-repeat_number_", nrow(repeats_chr), ".png", collapse = "")), 
        width = 6000, height = 10000, pointsize = 80)
    
    par(mar = c(4, 4, 4, 4), oma = c(1, 1, 1, 1))
    
    par(fig=c(0,1,0.4,1))
    ax.len <- nrow(repeats_chr)
    
    hors = hors[order(hors$SNV_per_kbp, decreasing = TRUE), ]
    ### FILTERING
    # hors <- hors[hors$SNV_per_kbp < 30,]
    
    
    unit.name <- "repeat ID"
    plot(x = NULL, y = NULL,
         xlab = paste0(chromosome_name, ", ", unit.name),
         ylab = paste0(chromosome_name, ", ", unit.name),
         xlim = c(0, ax.len),
         ylim = c(0, ax.len),
         pch = 19, cex = 0.1,
         main = paste0("HORs ", summary_df$assembly[i], " ", repeat_class, " ", chromosome_name, " ", nrow(repeats_chr), " repeats"))
    
    SNV_per_kbp_in_red = 20
    
    colours_SNV <- colorRampPalette(c("green", "yellow", "red"))(length(hors$SNV_per_kbp)) [findInterval(hors$SNV_per_kbp, seq(0, SNV_per_kbp_in_red, length.out = length(hors$SNV_per_kbp)))]
    while (TRUE) {
      if (ax.len >= 100000) {lwd_plot <- 2; break}
      if (ax.len >= 50000) {lwd_plot <- 3; break}
      if (ax.len >= 25000) {lwd_plot <- 4; break}
      if (ax.len >= 12500) {lwd_plot <- 5; break}
      lwd_plot <- 6
      break
    }
    for (k in seq_len(nrow(hors))) {
      lines(x = c(hors$start_A[k], hors$end_A[k]),
            y = c(hors$start_B[k], hors$end_B[k]),
            pch = 19, lwd = lwd_plot, col = colours_SNV[k])
    }
    
    
    ### check the most common HOR distance (between 3 and max_hor_distance)
    
    par(fig=c(0,1,0.2,0.4), new = TRUE)
    bins = 200
    
    hors_distance <- hors[(hors$start_B - hors$start_A) <= max_hor_distance,]
    
    hors_distances_counts <- rep(list(rep(0,bins)), max_hor_distance)
    hors_distances_bin_repeats <- ceiling((0:bins)*(nrow(repeats_chr)/bins))[-1] - ceiling((0:bins)*(nrow(repeats_chr)/bins))[-(bins+1)]
    
    hors_distance$dist <- hors_distance$start_B - hors_distance$start_A
    hors_distance <- hors_distance[order(hors_distance$start_A, decreasing = FALSE), ]
    for(k in seq_len(nrow(hors_distance))) {
      distance <- hors_distance$start_B[k] - hors_distance$start_A[k]
      for(l in 0 : (hors_distance$block.size.in.units[k]-1)) {
        which_bin <- ceiling((hors_distance$start_A[k] + l) / (nrow(repeats_chr) / bins))
        hors_distances_counts[[distance]][which_bin] = hors_distances_counts[[distance]][which_bin] + 1
        
        which_bin <- ceiling((hors_distance$start_B[k] + l) / (nrow(repeats_chr) / bins))
        hors_distances_counts[[distance]][which_bin] = hors_distances_counts[[distance]][which_bin] + 1
        
      }
      
    }
    plot(x = NULL, y = NULL, xlim = c(0,bins+1), ylim = c(0,max_hor_distance), 
         xlab = "repeat bin", ylab = "HOR distance frequency")
    for(k in 1 : max_hor_distance) {
      points(1:bins, rep(k,bins), cex = hors_distances_counts[[k]]/hors_distances_bin_repeats/3, pch = 15)
    }
    abline(h = 1:(max_hor_distance+1) - 0.5, v = 1:(bins+1)- 0.5, col = "gray")
    abline(h = seq(-0.5,max_hor_distance+0.5,by=5), col = "#0088ee", cex = 2)
    abline(h = seq(0.5,max_hor_distance+0.5,by=5), col = "#0088ee", cex = 2)
    abline(h = seq(-0.5,max_hor_distance+0.5,by=10), col = "#0000ee", cex = 4)
    abline(h = seq(0.5,max_hor_distance+0.5,by=10), col = "#0000ee", cex = 4)
    
    
    # Expand both block columns properly
    long_data1 <- data.table(
      element = unlist(lapply(1 : nrow(hors), function(X) rep(hors$start_A[X] : hors$end_A[X], hors$block.size.in.units[X])  )),
      partner = unlist(lapply(1 : nrow(hors), function(X) rep(hors$start_B[X] : hors$end_B[X], each = hors$block.size.in.units[X])   ))
    )
    # Expand both block columns properly
    long_data2 <- data.table(
      element = unlist(lapply(1 : nrow(hors), function(X) rep(hors$start_B[X] : hors$end_B[X], hors$block.size.in.units[X])  )),
      partner = unlist(lapply(1 : nrow(hors), function(X) rep(hors$start_A[X] : hors$end_A[X], each = hors$block.size.in.units[X])   ))
    )
    
    long_data1 <- unique(long_data1)
    long_data2 <- unique(long_data2)
    
    
    # Combine both data tables
    long_data <- rbind(long_data1, long_data2)
    long_data <- unique(long_data)
    
    remove(long_data1, long_data2)
    
    # Remove self-interactions
    long_data <- long_data[element != partner]
    # Compute unique interactors efficiently
    unique_interactions <- long_data[, .(num_interactors = uniqueN(partner)), by = element]
    
    # Create a result vector (0 for elements not in any block)
    result_vector <- integer(nrow(repeats_chr))
    # result_vector <- integer(nrow(repeats_chr))
    result_vector[unique_interactions$element] <- unique_interactions$num_interactors
    
    repeats_chr$HOR_score <- 100 * result_vector / nrow(repeats_chr)
    
    
    
    # add HOR score plot
    par(fig=c(0,1,0.0,0.2), new = TRUE)
    
    plot(y = repeats_chr$HOR_score, 
         x = seq_along(repeats_chr$HOR_score), type = "p", 
         pch = 16, col = "grey", ylim = c(0,100),
         xlab = "repeat ID", ylab = "HOR normalised count, %")
    points(y = ma(x = repeats_chr$HOR_score, n = round(nrow(repeats_chr) / bins)), 
           x = seq_along(repeats_chr$HOR_score), type = "l", ylim = c(0,100))
    
    dev.off()
    cat("\n")
    
    
    # calculate hors_formed_count
    repeats_chr$hors_formed_count = 0
    
    for(k in seq_len(nrow(hors))) {
      repeats_chr$hors_formed_count[hors$start_A[k] : hors$end_A[k]] = repeats_chr$hors_formed_count[hors$start_A[k] : hors$end_A[k]] + 1
      repeats_chr$hors_formed_count[hors$start_B[k] : hors$end_B[k]] = repeats_chr$hors_formed_count[hors$start_B[k] : hors$end_B[k]] + 1
    }
    repeats_chr$hors_formed_tot_rep_normalised = repeats_chr$hors_formed_count / nrow(repeats_chr)
    
    repeats_chr$hors_formed_tot_rep_normalised <- repeats_chr$hors_formed_tot_rep_normalised * 100
    
    write.csv(x = repeats_chr, file = paste0(hor_dirs[i], "/HOR_scored_", repeats_with_hors_files[j]), row.names = FALSE)
    
    
  }
  
  
  
  
}











