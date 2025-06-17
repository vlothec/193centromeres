## 0. Prepare the data and choose which species to sue
## 1. Make histograms of distance between repeats sizes to find out the best array size definition
## 2. Calculate internal similarity measures within arrays, between arrays same chr, between arrays 
##    different chr: boxplot for each species
## 3. Scatter plots of repeat content vs chromosome size and array number vs chromosome size
# 4. HOR scores averaged within arrays
## 5. Strand switching statistics

array_pairs_to_sample_for_similarity <- 1500
repeats_to_sample_per_array <- 10

library(seqinr)

chr_sizes <- read.csv(file = "/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full.Ian.csv")

species_to_analyse <- c("lpCarDepa1.1", "lpSchLacu1", "lpLuzSylv1.1", "iiLimLuna2.1", "iiLimMarm1.1", "iiLimRhom1.1")


repeats_to_analyse <- list("81_2", "183_4", c("124_1", "174_2"), "353_2", '166_3', "161_5", "155_1")


min_distances <- c(70000, 60000, 50000, 111111, 70000, 70000, 25000)

for(i in seq_along(species_to_analyse)) {
  cat(i, "\n")
  setwd(paste0("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs/", species_to_analyse[i], ".fa"))
  
  repeats <- read.csv(list.files(path = ".", pattern = "repeats_filtered.csv", full.names = TRUE, recursive = FALSE))
  
  repeats <- repeats[repeats$new_class %in% repeats_to_analyse[[i]], ]
  
  repeats <- repeats[order(repeats$start),]
  repeats <- repeats[order(repeats$seqID),]
  
  dist_to_next <- repeats$start[2:nrow(repeats)] - repeats$end[1 : (nrow(repeats) - 1)]
  dist_to_next <- dist_to_next[dist_to_next > 10000]
  dist_to_next <- dist_to_next[dist_to_next < 3000000]
  
  pdf(paste0(species_to_analyse[i], "_distances_between_repeats_histogram.pdf"), 
      width = 12)
  hist(dist_to_next, breaks = seq(0,3000000, by = 10000), xlim = c(0,3000000))
  abline(v = min_distances[i], col = "red")
  dev.off()
  
  pdf(paste0("/home/pwlodzimierz/ToL/upload_files/33_histograms_holocentric_distances_between_repeats/", species_to_analyse[i], "_distances_between_repeats_histogram.pdf"), 
      width = 12)
  hist(dist_to_next, breaks = seq(0,3000000, by = 10000), xlim = c(0,3000000))
  abline(v = min_distances[i], col = "red")
  dev.off()
  
}



for(i in seq_along(species_to_analyse)) {
  cat(i, "\n")
  for(z in 1 : length(repeats_to_analyse[[i]])) {
    
    #i=3;z=1
    setwd(paste0("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs/", species_to_analyse[i], ".fa"))
    
    chr_sizes <- read.csv(file = "/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full.Ian.csv")
    chr_sizes <- chr_sizes[chr_sizes$assembly.name == paste0(species_to_analyse[i], ".fa"), ]
    chr_sizes <- chr_sizes[chr_sizes$is.chr == 1, ]
    
    repeats <- read.csv(list.files(path = ".", pattern = "repeats_filtered.csv", full.names = TRUE, recursive = FALSE))
    
    repeats <- repeats[repeats$new_class %in% repeats_to_analyse[[i]][z], ]
    
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
    
    repeats$array_ID <- NA
    
    arrays <- NULL
    
    for(j in seq_along(array_starts)) {
      array_repeats <- repeats[array_starts[j] : array_ends[j], ]
      repeats$array_ID[array_starts[j] : array_ends[j]] = j
      arrays <- rbind(arrays, data.frame(chromosome = names(sort(table(array_repeats$seqID), decreasing = TRUE))[1], 
                                         start = min(array_repeats$start), 
                                         end = max(array_repeats$end), 
                                         rep_number = nrow(array_repeats), 
                                         length = max(array_repeats$end) - min(array_repeats$start), 
                                         repeat_total_length = sum(array_repeats$width), 
                                         plust_strand_perc = length(which(array_repeats$strand == "+")) / nrow(array_repeats),
                                         chr_size = chr_sizes$size[chr_sizes$chromosome.name == names(sort(table(array_repeats$seqID), decreasing = TRUE))[1]]))
    }
    
    
    arrays <- arrays[arrays$rep_number > 2,]
    
    write.csv(arrays, file = paste0(species_to_analyse[i], "_", repeats_to_analyse[[i]][z], "_holocentric_arrays.csv"), row.names = FALSE)
    
    write.csv(arrays, file = paste0("/home/pwlodzimierz/ToL/upload_files/33_holocentric_arrays/", species_to_analyse[i], "_", repeats_to_analyse[[i]][z], "_holocentric_arrays.csv"), row.names = FALSE)
    
    ### similarity within and between arrays ###
    
    {
      sample_arrays <- sample(repeats$array_ID, size = array_pairs_to_sample_for_similarity, replace = TRUE)
      
      
      similarity_data_frame <- data.frame(method = rep("same_array", array_pairs_to_sample_for_similarity),
                                          array_1 = sample_arrays,
                                          array_2 = sample_arrays,
                                          similarity_score = rep(0, array_pairs_to_sample_for_similarity))
      
      for(j in 1 : array_pairs_to_sample_for_similarity) {
        cat("array", i, j, "/", array_pairs_to_sample_for_similarity, "\n")
        repeats_A <- repeats[repeats$array_ID == sample_arrays[j], ]
        repeats_A_sample <- repeats_A[sample(1 : nrow(repeats_A), replace = TRUE, size = repeats_to_sample_per_array),]
        dist_matrix <- adist(repeats_A_sample$sequence)
        similarity_data_frame$similarity_score[j] <- 100-100*mean(dist_matrix[upper.tri(dist_matrix)])/mean(repeats_A_sample$width)
      }
      
      
      # Get indices by chromosome
      chr_groups <- split(1:nrow(arrays), arrays$chromosome)
      
      # Initialize storage
      same_chr_pairs <- matrix(NA, nrow = array_pairs_to_sample_for_similarity, ncol = 2)
      pair_count <- 0
      
      # Keep sampling until we have array_pairs_to_sample_for_similarity valid pairs
      while (pair_count < array_pairs_to_sample_for_similarity) {
        # Randomly choose a chromosome group
        chr <- sample(names(chr_groups), 1)
        idxs <- chr_groups[[chr]]
        
        # Only sample if there are at least 2 elements
        if (length(idxs) >= 2) {
          pair <- sample(idxs, 2)
          pair_count <- pair_count + 1
          same_chr_pairs[pair_count, ] <- pair
        }
      }
      
      similarity_data_frame <- rbind(similarity_data_frame, data.frame(method = rep("different_array_same_chr", array_pairs_to_sample_for_similarity),
                                                                       array_1 = rep(0, array_pairs_to_sample_for_similarity),
                                                                       array_2 = rep(0, array_pairs_to_sample_for_similarity),
                                                                       similarity_score = rep(0, array_pairs_to_sample_for_similarity)))
      
      for(j in 1 : array_pairs_to_sample_for_similarity) {
        cat("chr", i, j, "/", array_pairs_to_sample_for_similarity, "\n")
        repeats_A <- repeats[repeats$array_ID == same_chr_pairs[j,1], ]
        repeats_A_sample <- repeats_A[sample(1 : nrow(repeats_A), replace = TRUE, size = repeats_to_sample_per_array),]
        repeats_B <- repeats[repeats$array_ID == same_chr_pairs[j,2], ]
        repeats_B_sample <- repeats_B[sample(1 : nrow(repeats_B), replace = TRUE, size = repeats_to_sample_per_array),]
        similarity_data_frame$similarity_score[array_pairs_to_sample_for_similarity + j] = 
          100-100*mean(adist(repeats_A_sample$sequence, repeats_B_sample$sequence))/mean(c(repeats_A_sample$width, repeats_B_sample$width))
        similarity_data_frame$array_1[array_pairs_to_sample_for_similarity + j] = same_chr_pairs[j,1]
        similarity_data_frame$array_2[array_pairs_to_sample_for_similarity + j] = same_chr_pairs[j,2]
      }
      
      
      
      
      diff_chr_pairs <- matrix(NA, nrow = array_pairs_to_sample_for_similarity, ncol = 2)
      pair_count <- 0
      n <- nrow(arrays)
      
      while (pair_count < array_pairs_to_sample_for_similarity) {
        pair <- sample(1:n, 2)
        chr1 <- arrays$chromosome[pair[1]]
        chr2 <- arrays$chromosome[pair[2]]
        
        if (chr1 != chr2) {
          pair_count <- pair_count + 1
          diff_chr_pairs[pair_count, ] <- pair
        }
      }
      
      similarity_data_frame <- rbind(similarity_data_frame, data.frame(method = rep("different_array_different_chr", array_pairs_to_sample_for_similarity),
                                                                       array_1 = rep(0, array_pairs_to_sample_for_similarity),
                                                                       array_2 = rep(0, array_pairs_to_sample_for_similarity),
                                                                       similarity_score = rep(0, array_pairs_to_sample_for_similarity)))
      
      for(j in 1 : array_pairs_to_sample_for_similarity) {
        cat("genome", i, j, "/", array_pairs_to_sample_for_similarity, "\n")
        repeats_A <- repeats[repeats$array_ID == diff_chr_pairs[j,1], ]
        repeats_A_sample <- repeats_A[sample(1 : nrow(repeats_A), replace = TRUE, size = repeats_to_sample_per_array),]
        repeats_B <- repeats[repeats$array_ID == diff_chr_pairs[j,2], ]
        repeats_B_sample <- repeats_B[sample(1 : nrow(repeats_B), replace = TRUE, size = repeats_to_sample_per_array),]
        similarity_data_frame$similarity_score[j + array_pairs_to_sample_for_similarity + array_pairs_to_sample_for_similarity] = 
          100-100*mean(adist(repeats_A_sample$sequence, repeats_B_sample$sequence))/mean(c(repeats_A_sample$width, repeats_B_sample$width))
        similarity_data_frame$array_1[j + array_pairs_to_sample_for_similarity + array_pairs_to_sample_for_similarity] = diff_chr_pairs[j,1]
        similarity_data_frame$array_2[j + array_pairs_to_sample_for_similarity + array_pairs_to_sample_for_similarity] = diff_chr_pairs[j,2]
      }
      print("A1")
      xmin = 0
      xmax = 100
      colors <- c("#0072B2", "#E69F00", "#009E73")
      for(ymax in c(100,150,200,250,300,350,400)) {
        
        similarity_data_frame$similarity_score[similarity_data_frame$similarity_score < 0] = 0
        similarity_data_frame$similarity_score[similarity_data_frame$similarity_score > 100] = 100
        
        pdf(file = paste0(species_to_analyse[i], "_", repeats_to_analyse[[i]][z], "_", ymax, "_histogrmas_similarity_holocentric_within_outside.pdf"), width = 12, height = 12)
        par(mfrow = c(2,1))
        hist(similarity_data_frame$similarity_score[similarity_data_frame$method == "same_array"], 
             breaks = seq(xmin,xmax, by = 1), xlim = c(xmin, xmax), ylim = c(0,ymax), border = "#0072B290", col = "#0072B290",
             xlab = "per chromosome mean pairwise distance of repeats", main = "Histogram of centromeric pairwise distances averaged per chromosome")
        hist(similarity_data_frame$similarity_score[similarity_data_frame$method == "different_array_same_chr"], 
             breaks = seq(xmin,xmax, by = 1), xlim = c(xmin, xmax), ylim = c(0,ymax), add = TRUE, border = "#E69F0090", col = "#E69F0090")
        hist(similarity_data_frame$similarity_score[similarity_data_frame$method == "different_array_different_chr"], 
             breaks = seq(xmin,xmax, by = 1), xlim = c(xmin, xmax), ylim = c(0,ymax), add = TRUE, border = "#009E7390", col = "#009E7390")
        legend(x = xmax * 0.9, y = ymax*0.95, legend = c("arrays", "chromosomes", "genomes"), fill = c("#0072B290", "#E69F0090", "#009E7390"))
        
        # boxplot(similarity_data_frame$similarity_score[similarity_data_frame$method == "same_array"], 
        #         similarity_data_frame$similarity_score[similarity_data_frame$method == "different_array_same_chr"], 
        #         similarity_data_frame$similarity_score[similarity_data_frame$method == "different_array_different_chr"], 
        #         names = c("arrays", "chromosomes", "genomes"),
        #         main = "per chromosome mean pairwise similarity of repeats within:",
        #         col = c("#0072B290", "#E69F0090","#009E7390", ylim = c(0,100)))
        
        arrays2      <- similarity_data_frame$similarity_score[similarity_data_frame$method == "same_array"]
        chromosomes <- similarity_data_frame$similarity_score[similarity_data_frame$method == "different_array_same_chr"]
        genomes     <- similarity_data_frame$similarity_score[similarity_data_frame$method == "different_array_different_chr"]
        
        # Plot the boxplot
        boxplot(arrays2, chromosomes, genomes,
                names = c("arrays", "chromosomes", "genomes"),
                main = "per chromosome mean pairwise similarity of repeats within:",
                col = c("#0072B290", "#E69F0090", "#009E7390"),
                ylim = c(0, 140))
        
        # Perform pairwise t-tests
        p1 <- t.test(arrays2, chromosomes)$p.value
        p2 <- t.test(arrays2, genomes)$p.value
        p3 <- t.test(chromosomes, genomes)$p.value
        
        # Add asterisks based on significance levels
        # Positioning settings
        y_max <- max(c(arrays2, chromosomes, genomes), na.rm = TRUE)
        offset <- 5
        
        # Function to determine asterisk level
        get_asterisks <- function(p) {
          if (p < 0.001) return("***")
          if (p < 0.01)  return("**")
          if (p < 0.05)  return("*")
          return("n.s.")  # not significant
        }
        
        # Draw lines and asterisks
        segments(1, y_max + offset, 2, y_max + offset)
        text(1.5, y_max + offset + 2, get_asterisks(p1))
        
        segments(1, y_max + 2*offset, 3, y_max + 2*offset)
        text(2, y_max + 2*offset + 2, get_asterisks(p2))
        
        segments(2, y_max + 3*offset, 3, y_max + 3*offset)
        text(2.5, y_max + 3*offset + 2, get_asterisks(p3))
        dev.off()
        
        print("A2")
        
        pdf(file = paste0("/home/pwlodzimierz/ToL/upload_files/33_holocentric_similarity/", species_to_analyse[i], "_", repeats_to_analyse[[i]][z], "_", ymax, "_histogrmas_similarity_holocentric_within_outside.pdf"), width = 12, height = 12)
        par(mfrow = c(2,1))
        hist(similarity_data_frame$similarity_score[similarity_data_frame$method == "same_array"], 
             breaks = seq(xmin,xmax, by = 1), xlim = c(xmin, xmax), ylim = c(0,ymax), border = "#0072B290", col = "#0072B290",
             xlab = "per chromosome mean pairwise distance of repeats", main = "Histogram of centromeric pairwise distances averaged per chromosome")
        hist(similarity_data_frame$similarity_score[similarity_data_frame$method == "different_array_same_chr"], 
             breaks = seq(xmin,xmax, by = 1), xlim = c(xmin, xmax), ylim = c(0,ymax), add = TRUE, border = "#E69F0090", col = "#E69F0090")
        hist(similarity_data_frame$similarity_score[similarity_data_frame$method == "different_array_different_chr"], 
             breaks = seq(xmin,xmax, by = 1), xlim = c(xmin, xmax), ylim = c(0,ymax), add = TRUE, border = "#009E7390", col = "#009E7390")
        legend(x = xmax * 0.9, y = ymax*0.95, legend = c("arrays", "chromosomes", "genomes"), fill = c("#0072B290", "#E69F0090", "#009E7390"))
        
        # boxplot(similarity_data_frame$similarity_score[similarity_data_frame$method == "same_array"], 
        #         similarity_data_frame$similarity_score[similarity_data_frame$method == "different_array_same_chr"], 
        #         similarity_data_frame$similarity_score[similarity_data_frame$method == "different_array_different_chr"], 
        #         names = c("arrays", "chromosomes", "genomes"),
        #         main = "per chromosome mean pairwise similarity of repeats within:",
        #         col = c("#0072B290", "#E69F0090","#009E7390", ylim = c(0,100)))
        
        arrays2      <- similarity_data_frame$similarity_score[similarity_data_frame$method == "same_array"]
        chromosomes <- similarity_data_frame$similarity_score[similarity_data_frame$method == "different_array_same_chr"]
        genomes     <- similarity_data_frame$similarity_score[similarity_data_frame$method == "different_array_different_chr"]
        
        # Plot the boxplot
        boxplot(arrays2, chromosomes, genomes,
                names = c("arrays", "chromosomes", "genomes"),
                main = "per chromosome mean pairwise similarity of repeats within:",
                col = c("#0072B290", "#E69F0090", "#009E7390"),
                ylim = c(0, 140))
        
        # Perform pairwise t-tests
        p1 <- t.test(arrays2, chromosomes)$p.value
        p2 <- t.test(arrays2, genomes)$p.value
        p3 <- t.test(chromosomes, genomes)$p.value
        
        # Add asterisks based on significance levels
        # Positioning settings
        y_max <- max(c(arrays2, chromosomes, genomes), na.rm = TRUE)
        offset <- 5
        
        # Function to determine asterisk level
        get_asterisks <- function(p) {
          if (p < 0.001) return("***")
          if (p < 0.01)  return("**")
          if (p < 0.05)  return("*")
          return("n.s.")  # not significant
        }
        
        # Draw lines and asterisks
        segments(1, y_max + offset, 2, y_max + offset)
        text(1.5, y_max + offset + 2, get_asterisks(p1))
        
        segments(1, y_max + 2*offset, 3, y_max + 2*offset)
        text(2, y_max + 2*offset + 2, get_asterisks(p2))
        
        segments(2, y_max + 3*offset, 3, y_max + 3*offset)
        text(2.5, y_max + 3*offset + 2, get_asterisks(p3))
        dev.off()
        
      }
      
    } ### similarity within and between arrays ###
    
    print("### repeat content vs chromosome size and array number vs chromosome size ###")
    
    ### repeat content vs chromosome size and array number vs chromosome size ###
    
    {
      chr_sizes$repeat_array_number <- NA
      chr_sizes$repeat_bp <- NA
      
      for(j in 1 : nrow(chr_sizes)) {
        chr_arrays <- arrays[arrays$chromosome == chr_sizes$chromosome.name[j], ]
        chr_sizes$repeat_array_number[j] = nrow(chr_arrays)
        chr_sizes$repeat_bp[j] = sum(chr_arrays$repeat_total_length)
      }
      
      pdf(paste0(species_to_analyse[i], "_", repeats_to_analyse[[i]][z], "_scatter_array_no_vs_chr_size.pdf"))
      
      plot(chr_sizes$repeat_array_number, chr_sizes$size,
           xlab = "Repeat array number", ylab = "Chromosome size",
           main = "Repeat array number vs Chromosome size")
      
      # Add trendline
      model1 <- lm(size ~ repeat_array_number, data = chr_sizes)
      abline(model1, col = "blue", lwd = 2)
      
      # Correlation test and annotate
      cor_result1 <- cor.test(chr_sizes$repeat_array_number, chr_sizes$size)
      r1 <- round(cor_result1$estimate, 2)
      p1 <- cor_result1$p.value
      sig1 <- if (p1 < 0.001) "***" else if (p1 < 0.01) "**" else if (p1 < 0.05) "*" else "n.s."
      
      legend("topleft", legend = paste0("r = ", r1, ", p = ", signif(p1, 2), " (", sig1, ")"),
             bty = "n", text.col = "blue")
      
      dev.off()
      
      
      pdf(paste0(species_to_analyse[i], "_", repeats_to_analyse[[i]][z], "_scatter_repeat_bp_vs_chr_size.pdf"))
      
      plot(chr_sizes$repeat_bp, chr_sizes$size,
           xlab = "Repeat base pairs", ylab = "Chromosome size",
           main = "Repeat base pairs vs Chromosome size")
      
      # Add trendline
      model2 <- lm(size ~ repeat_bp, data = chr_sizes)
      abline(model2, col = "blue", lwd = 2)
      
      # Correlation test and annotate
      cor_result2 <- cor.test(chr_sizes$repeat_bp, chr_sizes$size)
      r2 <- round(cor_result2$estimate, 2)
      p2 <- cor_result2$p.value
      sig2 <- if (p2 < 0.001) "***" else if (p2 < 0.01) "**" else if (p2 < 0.05) "*" else "n.s."
      
      legend("topleft", legend = paste0("r = ", r2, ", p = ", signif(p2, 2), " (", sig2, ")"),
             bty = "n", text.col = "blue")
      
      dev.off()
      
      
      
      
      
      pdf(paste0("/home/pwlodzimierz/ToL/upload_files/33_scatters_holocentrics_against_chr_size/", species_to_analyse[i], "_", repeats_to_analyse[[i]][z], "_scatter_array_no_vs_chr_size.pdf"))
      
      plot(chr_sizes$repeat_array_number, chr_sizes$size,
           xlab = "Repeat array number", ylab = "Chromosome size",
           main = "Repeat array number vs Chromosome size")
      
      # Add trendline
      model1 <- lm(size ~ repeat_array_number, data = chr_sizes)
      abline(model1, col = "blue", lwd = 2)
      
      # Correlation test and annotate
      cor_result1 <- cor.test(chr_sizes$repeat_array_number, chr_sizes$size)
      r1 <- round(cor_result1$estimate, 2)
      p1 <- cor_result1$p.value
      sig1 <- if (p1 < 0.001) "***" else if (p1 < 0.01) "**" else if (p1 < 0.05) "*" else "n.s."
      
      legend("topleft", legend = paste0("r = ", r1, ", p = ", signif(p1, 2), " (", sig1, ")"),
             bty = "n", text.col = "blue")
      
      dev.off()
      
      
      pdf(paste0("/home/pwlodzimierz/ToL/upload_files/33_scatters_holocentrics_against_chr_size/", species_to_analyse[i], "_", repeats_to_analyse[[i]][z], "_scatter_repeat_bp_vs_chr_size.pdf"))
      
      plot(chr_sizes$repeat_bp, chr_sizes$size,
           xlab = "Repeat base pairs", ylab = "Chromosome size",
           main = "Repeat base pairs vs Chromosome size")
      
      # Add trendline
      model2 <- lm(size ~ repeat_bp, data = chr_sizes)
      abline(model2, col = "blue", lwd = 2)
      
      # Correlation test and annotate
      cor_result2 <- cor.test(chr_sizes$repeat_bp, chr_sizes$size)
      r2 <- round(cor_result2$estimate, 2)
      p2 <- cor_result2$p.value
      sig2 <- if (p2 < 0.001) "***" else if (p2 < 0.01) "**" else if (p2 < 0.05) "*" else "n.s."
      
      legend("topleft", legend = paste0("r = ", r2, ", p = ", signif(p2, 2), " (", sig2, ")"),
             bty = "n", text.col = "blue")
      
      dev.off()
      
      
      
      
      
      
      
      
      
      
      
      
    } ### repeat content vs chromosome size and array number vs chromosome size
    
    
    print("### strand switching ###")
    
    ### strand switching ###
    
    ## arrays <- read.csv("C:\\Users\\Piotr Włodzimierz\\Desktop\\ToL\\temp_data/lpCarDepa1.1_holocentric_arrays.csv")
    
    {
      arrays$array_strands <- "mixed"
      
      arrays$array_strands[arrays$plust_strand_perc >= 0.6] <- "+"
      arrays$array_strands[arrays$plust_strand_perc <= 0.4] <- "-"
      
      strand_vec <- arrays$array_strands
      strand_clean <- strand_vec[strand_vec %in% c("+", "-")]
      
      transitions <- paste0(head(strand_clean, -1), "->", tail(strand_clean, -1))
      
      # Tabulate transitions
      transition_table <- table(transitions)
      
      # Reformat into matrix
      transition_matrix <- matrix(0, nrow = 2, ncol = 2,
                                  dimnames = list(from = c("+", "-"), to = c("+", "-")))
      for (tr in names(transition_table)) {
        parts <- strsplit(tr, "->")[[1]]
        transition_matrix[parts[1], parts[2]] <- transition_table[tr]
      }
      chisq.test(transition_matrix)
      
      # capture.output(chisq.test(transition_matrix), file = paste0(species_to_analyse[i], ".txt"))
      # 
      # capture.output(chisq.test(transition_matrix), file = paste0("/home/pwlodzimierz/ToL/upload_files/33_chisq_test_results/", species_to_analyse[i], ".txt"))
      
      arrays$array_strands <- "mixed"
      
      arrays$array_strands[arrays$plust_strand_perc >= 0.6] <- "+"
      arrays$array_strands[arrays$plust_strand_perc <= 0.4] <- "-"
      
      strand_vec <- arrays$array_strands
      strand_clean <- strand_vec[strand_vec %in% c("+", "-")]
      
      strand_binary <- ifelse(strand_clean == "+", 1, 0)
      
      # diff calculates difference between two elements of the vector, so it can be 0 if there's no
      # strand switch, -1 or 1 when there's 0-1 or 1-0 switches. Sum counts how many non-zero 
      # values are in the vector, so how many strand transistions happened
      count_switches <- function(vec) {
        sum(diff(vec) != 0)
      }
      
      # Count number of switches
      observed_switches <- count_switches(strand_binary)
      
      # permutate and get switches
      n_simulations <- 10000
      simulated_switches <- replicate(n_simulations, {
        random_vec <- sample(strand_binary)  # keep the 0/1 ratio
        count_switches(random_vec)
      })
      
      hist(simulated_switches, breaks = 50, col = "lightgray",
           main = "Monte Carlo Simulation of Strand Switching",
           xlab = "Number of Switches")
      abline(v = observed_switches, col = "red", lwd = 2)
      
      p_value <- mean(abs(simulated_switches - mean(simulated_switches)) >= abs(observed_switches - mean(simulated_switches)))
      
      cat("Observed switches:", observed_switches, "\n")
      cat("Mean simulated switches:", mean(simulated_switches), "\n")
      cat("P-value:", p_value, "\n")
      
      text(x = 390, y = 500, labels = paste0("P-value: ", p_value))
      
    } ### repeat content vs chromosome size and array number vs chromosome size
  }
  
  
}






for(i in seq_along(species_to_analyse)) {
  cat(i, "\n")
  for(z in 1 : length(repeats_to_analyse[[i]])) {
    setwd(paste0("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs/", species_to_analyse[i], ".fa"))
    
    
    arrays = read.csv(paste0("/home/pwlodzimierz/ToL/upload_files/33_holocentric_arrays/", species_to_analyse[i], "_", repeats_to_analyse[[i]][z], "_holocentric_arrays.csv"))
    
    pdf(paste0("/home/pwlodzimierz/ToL/upload_files/33_holocentric_arrays/", species_to_analyse[i], "_", repeats_to_analyse[[i]][z], "_array_size_histogram.pdf"))
    
    hist(arrays$length[(arrays$length/1000) < 500]/1000, breaks = (0:100)*5,
         xlab = "Array length, Kbp", main = paste0(species_to_analyse[i], "_", repeats_to_analyse[[i]][z], "_array_size_histogram"))
    
    dev.off()
    
    
  }
}





