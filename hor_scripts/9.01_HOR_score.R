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
        chromosome_name, repeat_class, "\n")
    
    
    repeats_chr <- read.csv(file = paste0(hor_dirs[i], "/", repeats_with_hors_files[j]))
    
    if(file.exists(file.path("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs/HOR_periods",
                             paste0("new_HORs_lines_simple_", summary_df$assembly[i], "-", repeat_class, "-", chromosome_name, "-repeat_number_", nrow(repeats_chr), ".png", collapse = "")))) {
      next
    }
    
    if(!file.exists(paste0(hor_dirs[i], "/HORs_", repeat_class, "_", chromosome_name, ".csv"))) next
    hors <- read.csv(file = paste0(hor_dirs[i], "/HORs_", repeat_class, "_", chromosome_name, ".csv"))
    
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

