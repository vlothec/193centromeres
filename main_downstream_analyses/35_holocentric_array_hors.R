
library(seqinr)

chr_sizes <- read.csv(file = "/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full.Ian.csv")

species_to_analyse <- c("lpCarDepa1.1", "lpSchLacu1", "lpLuzSylv1.1", "iiLimLuna2.1", "iiLimMarm1.1", "iiLimRhom1.1")


repeats_to_analyse <- list("81_2", "183_4", c("124_1", "174_2"), "353_2", '166_3', "161_5")


min_distances <- c(70000, 60000, 50000, 111111, 70000, 70000, 25000)

library(foreach)
library(doParallel)
library(seqinr)
library(data.table)

# TODO: don't do scaffolds! For now, I removed them manually

#set up scripts

if(FALSE) {
  hor_scirpt = "/home/pwlodzimierz/ToL/git_ToL/35_islands_hors_submit_slurm.sh" 
  
  setwd("/home/pwlodzimierz/ToL/upload_files/35_holocentric_HOR_socres_per_array")
  if(file.exists("./hor_submission.sh")) {
    file.remove("./hor_submission.sh")
  }
  
  for(i in 1 : length(species_to_analyse)) {
    cat(i, "/", length(species_to_analyse), species_to_analyse[i],  "\n")
    
    
    for(irep in 1 : length(repeats_to_analyse[[i]])) {
      setwd(paste0("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs/", species_to_analyse[i], ".fa"))
      
      output_dir <- paste0(getwd(), "/array_HORs")
      
      
      repeats <- read.csv(list.files(path = ".", pattern = "repeats_filtered.csv", full.names = TRUE, recursive = FALSE))
      
      repeats <- repeats[repeats$new_class %in% repeats_to_analyse[[i]][irep], ]
      
      repeats <- repeats[order(repeats$start),]
      repeats <- repeats[order(repeats$seqID),]
      
      repeats <- repeats[repeats$seqID %in% chr_sizes$chromosome.name, ]
      
      dist_to_next <- repeats$start[2:nrow(repeats)] - repeats$end[1 : (nrow(repeats) - 1)]
      
      array_ends <- which(dist_to_next > min_distances[i])
      
      array_ends <- c(array_ends, which(dist_to_next < 0))
      
      array_ends <- c(array_ends, length(dist_to_next))
      
      repeats$end_array <- FALSE
      repeats$end_array[array_ends] <- TRUE
      repeats$end_array[nrow(repeats)] <- TRUE
      
      repeats$start_array <- FALSE
      repeats$start_array[array_ends + 1] <- TRUE
      repeats$start_array[1] <- TRUE
      
      repeats$width <- repeats$end - repeats$start + 1
      
      for(j in unique(repeats$seqID)) {
        repeats$start_array[repeats$seqID == j][1] = TRUE
        repeats$end_array[repeats$seqID == j][nrow(repeats[repeats$seqID == j,])] = TRUE
      }
      
      array_starts <- which(repeats$start_array)
      array_ends <- which(repeats$end_array)
      
      repeats$strand[repeats$strand == "+"] <- "1"
      repeats$strand[repeats$strand == "-"] <- "2"
      
      dir.create("./array_HORs")
      
      arrays <- NULL
      
      for(j in seq_along(array_starts)) {
        array_repeats <- repeats[array_starts[j] : array_ends[j], ]
          if(nrow(array_repeats) < 2) next
        
        arrays <- rbind(arrays, data.frame(chromosome = names(sort(table(array_repeats$seqID), decreasing = TRUE))[1], 
                                           start = min(array_repeats$start), 
                                           end = max(array_repeats$end), 
                                           rep_number = nrow(array_repeats), 
                                           length = max(array_repeats$end) - min(array_repeats$start), 
                                           repeat_total_length = sum(array_repeats$width), 
                                           plust_strand_perc = length(which(array_repeats$strand == "+")) / nrow(array_repeats)))
        
        array_repeats$seqID <- unlist(lapply(array_repeats$seqID, function(X) paste0(X, "_island", j)))

        repeats_file <- paste0(getwd(), "/array_HORs/", species_to_analyse[i], "_", repeats_to_analyse[[i]][irep], "_array_", j, ".csv")

        if(!file.exists(paste0(output_dir, "/repeats_with_hors_", repeats_to_analyse[[i]][irep], "_", array_repeats$seqID[1], ".csv"))) {
          write.csv(x = array_repeats, file = repeats_file, row.names = FALSE)
          
          cat(hor_scirpt, output_dir, repeats_file, array_repeats$seqID[1], repeats_to_analyse[[i]][irep], repeats_file, array_repeats$seqID[1], repeats_to_analyse[[i]][irep], paste0(species_to_analyse[i], ".fa"), paste0(species_to_analyse[i], ".fa"), "\n",
              file = "/home/pwlodzimierz/ToL/upload_files/35_holocentric_HOR_socres_per_array/hor_submission.sh", append = TRUE)
          cat("sleep 0.01", "\n",
              file = "/home/pwlodzimierz/ToL/upload_files/35_holocentric_HOR_socres_per_array/hor_submission.sh", append = TRUE)
          
        }

        
      }
      
      write.csv(file = "./holocentric_arrays_data.csv", x = arrays, row.names = FALSE)
      
      
    }
    
  }
}


if(TRUE) {
 
  for(i in 1 : length(species_to_analyse)) {
    cat(i, "/", length(species_to_analyse), species_to_analyse[i],  "\n")
    
    
    for(irep in 1 : length(repeats_to_analyse[[i]])) {
      setwd(paste0("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs/", species_to_analyse[i], ".fa/array_HORs/"))
      
      hor_files <- list.files(path = "./", pattern = paste0("HORs_", repeats_to_analyse[[i]][irep]))
      hor_files <- hor_files[!grepl("repeats", hor_files)]
      
      
      windows <- 30
      score_per_window <- rep(0, windows)
      repeats_per_window <- rep(0, windows)
      
      for(j in seq_along(hor_files)) {
        cat(i, j, "/", length(hor_files), "\n")
        hors <- read.csv(file = hor_files[j])
        
        arrayID <- strsplit(strsplit(hor_files[j], split = "island")[[1]][2], split = ".csv")[[1]][1]
        chrID <- strsplit(strsplit(hor_files[j], split = paste0(repeats_to_analyse[[i]][irep], "_"))[[1]][2], split = "_island")[[1]][1]
        
        repeats_file <- paste0(getwd(), "/repeats_with_hors_", repeats_to_analyse[[i]][irep], "_", chrID, "_island", arrayID, ".csv")
        
        if(file.exists(repeats_file)) {
          repeats <- read.csv(file = repeats_file)
          repeats$novel_score <- 0
          
          for(k in seq_len(nrow(hors))) {
            repeats$novel_score[hors$start_A[k] : hors$end_A[k]] = repeats$novel_score[hors$start_A[k] : hors$end_A[k]] + 1
            repeats$novel_score[hors$start_B[k] : hors$end_B[k]] = repeats$novel_score[hors$start_B[k] : hors$end_B[k]] + 1
          }
          repeats$novel_score = repeats$novel_score / nrow(repeats)
          
          repeats$novel_score <- repeats$novel_score * 100
          
          
          indices <- 1 : length(repeats$novel_score)
          bin_size <- floor(length(repeats$novel_score) / windows)
          remainder <- length(repeats$novel_score) %% windows
          bin_sizes <- rep(bin_size, windows)
          if (remainder > 0) {
            extra_bins <- sample(1 : windows, remainder, replace = FALSE)
            bin_sizes[extra_bins] <- bin_sizes[extra_bins] + 1
          }
          bin_assignments <- rep(1 : windows, times = bin_sizes)
          
          bin_assignments <- bin_assignments[1 : length(repeats$novel_score)]
          
          bins <- split(indices, bin_assignments)
          
          for(k in seq_along(bins)) {
            score_per_window[k] <- score_per_window[k] + sum(repeats$novel_score[bins[k][[1]]])
            repeats_per_window[k] <- repeats_per_window[k] + length(bins[k][[1]])
          }
          
        }
        
      }
      
      averaged_scores <- score_per_window / repeats_per_window
      
      pdf(file = paste0("/home/pwlodzimierz/ToL/upload_files/35_holocentric_HOR_socres_per_array/", species_to_analyse[i], "_", repeats_to_analyse[[i]][irep], "_HOR_island_averaged_values.pdf"), 
          width = 7, height = 3)
      plot(x = 1 : windows, y = averaged_scores, ylim = c(min(averaged_scores)*0.95, max(averaged_scores)*1.05),
           main = paste0(repeats_to_analyse[[i]][irep], " HOR scores per island"), 
           xlab = paste0(windows, " windows averaged scores of ", length(hor_files), " islands, totaling ", sum(repeats_per_window), " repeats"),
           ylab = "HOR score: % of other repeats a repeat forms a HOR with", 
           cex.lab = 0.7, cex.axis = 0.7, cex.main = 1, cex.sub = 0.5, 
           type = "b")
      dev.off()
      
    }
    
  }
}













  