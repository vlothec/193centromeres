.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))
suppressMessages(library(msa))
suppressMessages(library(seqinr))
suppressMessages(library(GenomicRanges))
suppressMessages(library(ggplot2))

setwd("/home/pwlodzimierz/ToL/git_ToL")
source("./aux_fun.R")
ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}
replace_existing_analysis = TRUE

data_directories <- list.dirs(path = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs", recursive = FALSE, full.names = TRUE)
data_directories <- data_directories[!grepl(pattern = "templated_", data_directories)]
data_directories <- data_directories[grepl(pattern = ".fa", data_directories)]
assembly_files <- list.files(path = "/home/pwlodzimierz/ToL/Assemblies/fastas_2021_Michael", recursive = FALSE, full.names = TRUE)
assembly_files <- assembly_files[!grepl(".fai", assembly_files)]

if(FALSE) { # dev settings
  setwd("C:\\Users\\Piotr Włodzimierz\\Documents\\GitHub\\ToL\\")
  repeats = read.csv(file = ".\\TRASH_data\\ddAraThal4.1.fa\\ddAraThal4.1_repeats_filtered.csv")
  arrays = read.csv(file = ".\\TRASH_data\\ddAraThal4.1.fa\\ddAraThal4.1_arrays_filtered.csv")
  classes = read.csv(file = ".\\TRASH_data\\ddAraThal4.1.fa\\ddAraThal4.1_classes_merged_filtered.csv")
  classes$num_ID <- 1 : nrow(classes)
  edta = read.csv(file = ".\\TRASH_data\\ddAraThal4.1.fa\\ddAraThal4.1.fa_edta_modified.csv")
  assembly_file <- "ddAraThal4.1.fa"
  genes <- read.table(file = ".\\TRASH_data\\ddAraThal4.1.fa\\ddAraThal4_helixer.gff", header = FALSE, sep = "\t", skip = 4)
  genes <- genes[genes$V3 == "gene", ]
}

taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 15
print(i)

print(paste0(i, " / ", length(data_directories)))
### Load data ================================================================
setwd(data_directories[i])
print(getwd())

assembly_name = strsplit(strsplit(data_directories[i], split = ".fa")[[1]][1], split = "v2_out_for_HORs/")[[1]][2]
assembly_file = grep(assembly_name, assembly_files)
print(assembly_file)
print(assembly_files[assembly_file])


if(!replace_existing_analysis) {
  if(file.exists(paste0("genome_classes_", assembly_name, ".csv"))) {print("Analysis finished"); quit(save = "no", status = 1)}
} else {
  if(file.exists(paste0("genome_classes_", assembly_name, ".csv"))) {
    file.remove(paste0("genome_classes_", assembly_name, ".csv"))
  }
}

helixer_file = list.files(pattern = "helixer_filtered.gff", full.names = TRUE)
if(length(helixer_file) != 1) {warning(paste0(i, "No genes!")); setwd(".."); quit(save = "no", status = 1)}

repeat_file = list.files(pattern = "_repeats_filtered.csv", full.names = TRUE)
if(length(repeat_file) != 1) {warning(paste0(i, "No repeats!")); setwd(".."); quit(save = "no", status = 1)}

array_file = list.files(pattern = "_arrays_filtered.csv", full.names = TRUE)
if(length(array_file) != 1) {warning(paste0(i, " no arrays!")); setwd(".."); quit(save = "no", status = 1)}

classes_file = list.files(pattern = "_classes_merged_filtered", full.names = TRUE)
if(length(classes_file) != 1) {warning(paste0(i, " no classes!")); setwd(".."); quit(save = "no", status = 1)}

edta_file = list.files(pattern = paste0(assembly_name, "_edta_filtered.csv"), full.names = TRUE)
edta_file = edta_file[!grepl("reassigned", edta_file)]
if(length(edta_file) != 1) {warning(paste0(i, " no edta!")); setwd(".."); quit(save = "no", status = 1)}



print("load annotations")
repeats = read.csv(file = repeat_file)
arrays = read.csv(file = array_file)
classes = read.csv(file = classes_file)
classes$num_ID <- 1 : nrow(classes)
edta = read.csv(file = edta_file)
genes <- read.table(file = helixer_file, header = FALSE, sep = "\t", skip = 4)
genes <- genes[genes$V3 == "CDS", ]
print("annotations loaded")

arrays$overlapping_bp = 0
edta$overlapping_bp = 0

edta_filtered_total = NULL

### filter data to include only chromosomes ====================================

print("filter chromosomes")

all_chromosomes <- read.csv("/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full.Ian.csv")
all_chromosomes <- all_chromosomes[all_chromosomes$assembly.name == strsplit(data_directories[i], split = "v2_out_for_HORs/")[[1]][2], ]
all_chromosomes <- all_chromosomes[all_chromosomes$is.chr == 1, ]
chromosomes <- all_chromosomes$chromosome.name
chromosomes_lengths <- all_chromosomes$size

# ==============================================================================

# Get top 5 classes of each chromosome =========================================
print("top 5 for chrs")
top_5_unique <- NULL
for (j in seq_along(chromosomes)) {  
  sequence_repeats = repeats[repeats$seqID == chromosomes[j], ]
  if(nrow(sequence_repeats) == 0) next
  chr_classes <- as.data.frame(table(sequence_repeats$new_class))
  names(chr_classes) <- c("class", "count")
  chr_classes$class <- as.character(chr_classes$class)
  chr_classes$mean_length <- 0
  chr_classes$total_bp <- 0
  for(k in seq_len(nrow(chr_classes))) {
    chr_classes$mean_length[k] <- mean(sequence_repeats$width[sequence_repeats$new_class == chr_classes$class[k]])
    chr_classes$total_bp[k] <- sum(sequence_repeats$width[sequence_repeats$new_class == chr_classes$class[k]])
  }
  chr_classes <- chr_classes[chr_classes$mean_length >= 8, ]
  if(nrow(chr_classes) == 0) next
  chr_classes <- chr_classes[order(chr_classes$total_bp, decreasing = TRUE), ]
  top_5_unique <- unique(c(top_5_unique, chr_classes$class[1 : ifelse(nrow(chr_classes) >= 5, yes = 5, no = nrow(chr_classes))]))
}
# ==============================================================================


genome_classes <- NULL
print("genome classes")
pdf(file = paste0(assembly_name, "_chrs_plots_gene_TE.pdf"), onefile = TRUE)
for (j in seq_along(chromosomes)) {  
  print(j)
  sequence_arrays = arrays[arrays$seqID == chromosomes[j], ]
  sequence_edta = edta[edta$V1 == chromosomes[j], ]
  sequence_repeats = repeats[repeats$seqID == chromosomes[j], ]
  sequence_genes <- genes[genes$V1 == chromosomes[j],]
  
  if((nrow(sequence_arrays) == 0) | (nrow(sequence_edta) == 0) | (nrow(sequence_genes) == 0)) next
  
  edta_filtered = sequence_edta
  
  gr2 <- with(sequence_edta, GRanges(chromosomes[j], IRanges(V4, V5)))
  
  ### Calculate TE and Gene landscape =====================================================
  print("calculate TE landscape")
  #   Calculation:  REMOVE ALL REPEAT COORDINATES FROM THE CHROMOSOME  -> 
  sequence_edta_no_rep <- edta_filtered
  sequence_edta_no_rep$legacy_V4 <- sequence_edta_no_rep$V4
  sequence_edta_no_rep$legacy_V5 <- sequence_edta_no_rep$V5
  sequence_edta_no_rep$adjustment <- 0
  
  sequence_genes_no_rep <- sequence_genes
  sequence_genes_no_rep$legacy_V4 <- sequence_genes_no_rep$V4
  sequence_genes_no_rep$legacy_V5 <- sequence_genes_no_rep$V5
  sequence_genes_no_rep$adjustment <- 0
  
  sequence_repeats_t <- sequence_repeats[order(sequence_repeats$start, decreasing = TRUE), ]
  sequence_repeats_t$adjustment <- 0
  
  for(repeat_id in seq_len(nrow(sequence_repeats_t))) {
    sequence_edta_no_rep$adjustment[sequence_edta_no_rep$V4 >= sequence_repeats_t$start[repeat_id]] <- sequence_edta_no_rep$adjustment[sequence_edta_no_rep$V4 >= sequence_repeats_t$start[repeat_id]] + sequence_repeats_t$width[repeat_id]
    sequence_genes_no_rep$adjustment[sequence_genes_no_rep$V4 >= sequence_repeats_t$start[repeat_id]] <- sequence_genes_no_rep$adjustment[sequence_genes_no_rep$V4 >= sequence_repeats_t$start[repeat_id]] + sequence_repeats_t$width[repeat_id]
    sequence_repeats_t$adjustment[sequence_repeats_t$start >= sequence_repeats_t$start[repeat_id]] <- sequence_repeats_t$adjustment[sequence_repeats_t$start >= sequence_repeats_t$start[repeat_id]] + sequence_repeats_t$width[repeat_id]
  }
  sequence_edta_no_rep$V4 <- sequence_edta_no_rep$V4 - sequence_edta_no_rep$adjustment
  sequence_edta_no_rep$V5 <- sequence_edta_no_rep$V5 - sequence_edta_no_rep$adjustment
  
  sequence_genes_no_rep$V4 <- sequence_genes_no_rep$V4 - sequence_genes_no_rep$adjustment
  sequence_genes_no_rep$V5 <- sequence_genes_no_rep$V5 - sequence_genes_no_rep$adjustment
  
  sequence_repeats_t$start_adj <- sequence_repeats_t$start - sequence_repeats_t$adjustment
  sequence_repeats_t$end_adj <- sequence_repeats_t$end - sequence_repeats_t$adjustment
  
  # Make a TE all bp position map
  TE_coordinates <- list()
  for(edta_id in seq_len(nrow(sequence_edta_no_rep))) {
    TE_coordinates <- append(TE_coordinates, list(sequence_edta_no_rep$V4[edta_id] : sequence_edta_no_rep$V5[edta_id]))
  }
  TE_coordinates <- unlist(TE_coordinates)
  length(TE_coordinates)
  length(unique(TE_coordinates))
  
  # find arithmetic mean of all scores as the “peak”
  hist_EDTA <- hist(TE_coordinates, breaks = seq(min(TE_coordinates), max(TE_coordinates), length.out = 25), plot = FALSE)
  counts <- c(hist_EDTA$counts[1], hist_EDTA$counts[1], hist_EDTA$counts, hist_EDTA$counts[length(hist_EDTA$counts)], hist_EDTA$counts[length(hist_EDTA$counts)])
  ma_values <- ma(counts)[3 : (length(counts) - 2)]
  edta_peak <- hist_EDTA$mids[which.max(ma_values)]
  plot(y = ma_values, x = hist_EDTA$mids, type = "h", lwd = 4, ylim = c(0, max(ma_values)),
       main = paste0(chromosomes[j], " EDTA (minus repeats) landscape"))
  abline(v = edta_peak, col = "red")
  
  # calculate TE bp positions density in 2% chromosome length bins
  # calculate window TE density vs distance to peak and test it’s correlation
  distance_to_mid <- TE_coordinates
  distance_to_mid[distance_to_mid <= edta_peak] = edta_peak - distance_to_mid[distance_to_mid <= edta_peak]
  distance_to_mid[distance_to_mid >= edta_peak] = distance_to_mid[distance_to_mid >= edta_peak] - edta_peak
  # double the values that are on the longer arm, to normalise for arm length discrepancy
  chromosome_no_rep_size <- chromosomes_lengths[j] - sum(sequence_repeats$width)
  if(edta_peak < (chromosome_no_rep_size / 2)) { # if right arm is longer
    distance_to_mid <- c(distance_to_mid, 
                         distance_to_mid[distance_to_mid > (edta_peak)])
  } else { # if left arm is longer
    distance_to_mid <- c(distance_to_mid, 
                         distance_to_mid[distance_to_mid > (chromosome_no_rep_size - edta_peak)])
  }
  
  dist_to_mid_hist <- hist(distance_to_mid, breaks = 25)
  dist_to_mid_hist <- data.frame(dist_to_mid = log10(100 * dist_to_mid_hist$mids / chromosome_no_rep_size),
                                 counts = (dist_to_mid_hist$counts) / (dist_to_mid_hist$breaks[2:length(dist_to_mid_hist$breaks)] - dist_to_mid_hist$breaks[1:(length(dist_to_mid_hist$breaks) - 1)]) / 2)
  dist_to_mid_hist <- dist_to_mid_hist[dist_to_mid_hist$counts != -Inf,]
  summary(lm(counts ~ dist_to_mid, dist_to_mid_hist))
  plot(dist_to_mid_hist, ylim = c(0, 5))
  abline(lm(counts ~ dist_to_mid, dist_to_mid_hist))
  lm_coef <- lm(counts ~ dist_to_mid, dist_to_mid_hist)
  lm_coef$coefficients[1] + lm_coef$model$dist_to_mid[1] * lm_coef$coefficients[2]
  
  # Now calculations for GENES
  
  gene_coordinates <- list()
  for(edta_id in seq_len(nrow(sequence_genes_no_rep))) {
    gene_coordinates <- append(gene_coordinates, list(sequence_genes_no_rep$V4[edta_id] : sequence_genes_no_rep$V5[edta_id]))
  }
  gene_coordinates <- unlist(gene_coordinates)
  length(gene_coordinates)
  length(unique(gene_coordinates))
  
  hist_gene <- hist(gene_coordinates, breaks = seq(min(gene_coordinates), max(gene_coordinates), length.out = 25), plot = FALSE)
  counts <- c(hist_gene$counts[1], hist_gene$counts[1], hist_gene$counts, hist_gene$counts[length(hist_gene$counts)], hist_gene$counts[length(hist_gene$counts)])
  ma_values <- ma(counts)[3 : (length(counts) - 2)]
  gene_valley <- hist_gene$mids[which.min(ma_values)]
  plot(y = ma_values, x = hist_gene$mids, type = "h", lwd = 4, ylim = c(0, max(ma_values)),
       main = paste0(chromosomes[j], " genes (minus repeats) landscape"))
  abline(v = gene_valley, col = "red")
  
  
  distance_to_mid <- gene_coordinates
  distance_to_mid[distance_to_mid <= gene_valley] = gene_valley - distance_to_mid[distance_to_mid <= gene_valley]
  distance_to_mid[distance_to_mid >= gene_valley] = distance_to_mid[distance_to_mid >= gene_valley] - gene_valley

  if(gene_valley < (chromosome_no_rep_size / 2)) { # if right arm is longer
    distance_to_mid <- c(distance_to_mid, 
                         distance_to_mid[distance_to_mid > (gene_valley)])
  } else { # if left arm is longer
    distance_to_mid <- c(distance_to_mid, 
                         distance_to_mid[distance_to_mid > (chromosome_no_rep_size - gene_valley)])
  }
  
  dist_to_mid_hist <- hist(distance_to_mid, breaks = 25)
  dist_to_mid_hist <- data.frame(dist_to_mid = log10(100 * dist_to_mid_hist$mids / chromosome_no_rep_size),
                                 counts = (dist_to_mid_hist$counts) / (dist_to_mid_hist$breaks[2:length(dist_to_mid_hist$breaks)] - dist_to_mid_hist$breaks[1:(length(dist_to_mid_hist$breaks) - 1)]) / 2)
  dist_to_mid_hist <- dist_to_mid_hist[dist_to_mid_hist$counts != -Inf,]
  summary(lm(counts ~ dist_to_mid, dist_to_mid_hist))
  plot(dist_to_mid_hist, ylim = c(0, 2))
  abline(lm(counts ~ dist_to_mid, dist_to_mid_hist))
  lm_coef_genes <- lm(counts ~ dist_to_mid, dist_to_mid_hist)
  
  
  
  # ==============================================================================
  
  
  ### Get top 10 repeats =======================================================
  print("get top 10 repeats")
  sequence_repeats = repeats[repeats$seqID == chromosomes[j], ]
  sequence_repeats <- sequence_repeats[sequence_repeats$width > 7, ] # this removes classes smaller than 8 bp
  chr_classes <- as.data.frame(table(sequence_repeats$new_class))
  names(chr_classes) <- c("class", "count")
  chr_classes$class <- as.character(chr_classes$class)
  chr_classes$mean_length <- 0
  chr_classes$total_bp <- 0
  chr_classes$new_class_num_ID <- 0
  for(k in seq_len(nrow(chr_classes))) {
    chr_classes$mean_length[k] <- mean(sequence_repeats$width[sequence_repeats$new_class == chr_classes$class[k]])
    chr_classes$total_bp[k] <- sum(sequence_repeats$width[sequence_repeats$new_class == chr_classes$class[k]])
    chr_classes$new_class_num_ID[k] <- classes$num_ID[classes$class == chr_classes$class[k]]
  }
  if(nrow(chr_classes) == 0) next
  chr_classes <- chr_classes[order(chr_classes$total_bp, decreasing = TRUE), ]
  top_10 <- chr_classes$class[1 : ifelse(nrow(chr_classes) >= 10, yes = 10, no = nrow(chr_classes))]
  
  ### Check top 10 contains top 5 of other chromosomes =========================
  top_classes_chr <- unique(c(top_10, top_5_unique))
  chr_classes <- chr_classes[chr_classes$class %in% top_classes_chr, ]
  
  ### Calculate predictor values ===============================================
  print("predictors")
  chr_classes$total_bp_norm_chr <- -1
  chr_classes$total_bp_norm_rep <- -1
  chr_classes$start_sd_norm_chr <- -1
  chr_classes$start_norm_chr_0_50 <- -1
  chr_classes$gaps_count <- -1
  chr_classes$gaps_with_TEs_fraction <- -1
  chr_classes$centre_array_edit <- -1
  chr_classes$centre_array_width_sd <- -1
  chr_classes$centre_chromosome_edit <- -1
  chr_classes$centre_chromosome_width_sd <- -1
  chr_classes$array_sizes_sd_norm_mean_arr_size <- -1
  chr_classes$array_count <- -1
  chr_classes$TE_prox_dist <- -1
  chr_classes$TE_prox_SD <- -1
  chr_classes$TE_lm_coef <- -1
  chr_classes$TE_prox_score <- -1
  chr_classes$gene_prox_dist <- -1
  chr_classes$gene_prox_SD <- -1
  chr_classes$gene_lm_coef <- -1
  chr_classes$gene_prox_score <- -1
  chr_classes$pred_centrophilic_TE <- "LTR/Gypsy"
  chr_classes$t_test_t_val <- -1
  chr_classes$t_test_p_val <- -1
  
  for(k in seq_len(nrow(chr_classes))) {
    cat(k, "/", nrow(chr_classes), "\n")
    chr_family_repeats <- sequence_repeats[sequence_repeats$new_class == chr_classes$class[k], ]
    chr_family_repeats_t <- sequence_repeats_t[sequence_repeats_t$new_class == chr_classes$class[k], ]
    chr_family_arrays <- arrays[arrays$new_class_num_ID == chr_classes$new_class_num_ID[k] & arrays$seqID == chromosomes[j], ]
    # What is the repeat size?
    #   Calculation: Mean repeat length
    chr_classes$mean_length[k] <- chr_classes$mean_length[k] ### PREDICTOR ###
    
    # How big is the family?
    #   Calculation 1: Total bp normalized by chromosome length
    #   Calculation 2: Total bp normalized by all chromosome repeats bp
    chr_classes$total_bp_norm_chr[k] <- chr_classes$total_bp[k] / chromosomes_lengths[j] ### PREDICTOR ###
    chr_classes$total_bp_norm_rep[k] <- chr_classes$total_bp[k] / sum(sequence_repeats$width) ### PREDICTOR ###
    
    # How concentrated is the family? Low score could be holocentric!
    #   Calculation: Start positions SD normalised by chromosome length
    chr_classes$start_sd_norm_chr[k] <- sd(chr_family_repeats$start) / chromosomes_lengths[j] ### PREDICTOR ###
    
    # Where are the repeats?
    #   Calculation: Mean start positions value normalized by chromosome length, presented as distance to the closest chromosome edge (values 0:50)
    chr_classes$start_norm_chr_0_50[k] <- mean(chr_family_repeats$start) / chromosomes_lengths[j] ### PREDICTOR ###
    if(chr_classes$start_norm_chr_0_50[k] > 0.5) chr_classes$start_norm_chr_0_50[k] <- 1 - chr_classes$start_norm_chr_0_50[k] ### PREDICTOR ###
    
    # What are the array interspersed elements (if any)?
    #   Calculation: Find array interspersed sequence (under 20 kbp gaps, more than 200 bp) and count what fraction of them are Tes or TE-derived sequences
    if(nrow(chr_family_repeats) > 1) {
      gaps_start <- chr_family_repeats$end[1 : (nrow(chr_family_repeats) - 1)] + 1
      gaps_size <- chr_family_repeats$start[2 : (nrow(chr_family_repeats) - 0)] - chr_family_repeats$end[1 : (nrow(chr_family_repeats) - 1)] - 1
      gaps_start <- gaps_start[gaps_size >= 200 & gaps_size <= 20000]
      gaps_size <- gaps_size[gaps_size >= 200 & gaps_size <= 20000]
      chr_classes$gaps_count[k] <- length(gaps_size) ### PREDICTOR ###
      if(length(gaps_start) != 0) {
        gaps_filled <- rep(FALSE, length(gaps_start))
        for(i_gap in seq_along(gaps_size)) {
          gr3 <- GRanges(chromosomes[j], IRanges(gaps_start[i_gap], (gaps_start[i_gap] + gaps_size[i_gap])))
          
          overlaps <- as.data.frame(findOverlaps(gr3, gr2)) # gr2 are EDTA annotations
          if(nrow(overlaps) != 0) {
            gap_starts <- sequence_edta$V4[overlaps$subjectHits]
            gap_ends <- sequence_edta$V5[overlaps$subjectHits]
            gap_overlaps <- unlist(lapply(seq_along(gap_starts), function(X) sum( (gap_starts[X] : gap_ends[X]) %in% (gaps_start[i_gap] : (gaps_start[i_gap] + gaps_size[i_gap])) )))
            gap_overlaps_fraction <- gap_overlaps / gaps_size[i_gap]
            if(sum(gap_overlaps_fraction) > 0.25) gaps_filled[i_gap] <- TRUE
          }
        
      }
      chr_classes$gaps_with_TEs_fraction[k] <- sum(gaps_filled) / chr_classes$gaps_count[k] ### PREDICTOR ###
      }
    }
    min_repeats_to_align <- 100
    max_repeats_to_align <- 10000
    desired_fraction_to_align <- 0.2
    
    number_of_repeats_scored <- 0
    cumulative_adist_score <- 0
    mean_centre_width_SD <- NULL
    centre_repeats_count <- NULL
    arrays_sizes <- NULL
    cat(nrow(chr_family_arrays), "arrays:")
    for(array_id in seq_len(nrow(chr_family_arrays))) {
      # What is the repeat divergence within the central parts of the array?
      #   Calculation: central 50% repeats edit distance to THEIR consensus normalized by repeat mean length
      array_repeats <- chr_family_repeats[chr_family_repeats$arrayID == chr_family_arrays$arrayID[array_id], ]
      if(nrow(array_repeats) == 0) {
        arrays_sizes <- c(arrays_sizes, 0)
      } else {
        arrays_sizes <- c(arrays_sizes, sum(array_repeats$width))
      }
      
      if(nrow(array_repeats) < 10) next
      central_repeats <- array_repeats[ceiling(nrow(array_repeats) / 4) : floor(3 * nrow(array_repeats) / 4), ]
      
      repeats_to_align_IDs <- sample(seq_len(nrow(central_repeats)), round(nrow(central_repeats) * desired_fraction_to_align))
      if(length(repeats_to_align_IDs) > max_repeats_to_align) {
        repeats_to_align_IDs <- sample(seq_len(nrow(central_repeats)), max_repeats_to_align)
      } else if(length(repeats_to_align_IDs) < min_repeats_to_align) {
        if(nrow(central_repeats) > min_repeats_to_align) {
          repeats_to_align_IDs <- sample(seq_len(nrow(central_repeats)), min_repeats_to_align)
        } else {
          repeats_to_align_IDs <- seq_len(nrow(central_repeats))
        }
      }
      
      sequences_to_align <- central_repeats$sequence[repeats_to_align_IDs]
      cat(array_id, " (", nrow(array_repeats), " reps), ", sep = "")
      if(length(sequences_to_align) < 2) {
        cat("\n\n\n\n")
        stop(paste0(assembly_file, " did not find repeats in one of the classes: ", classes$class[j], ", investigate"))
      }
      # a <- capture.output({alignment_matrix = tolower(as.matrix(msa(sequences_to_align, method = "ClustalOmega", type = "dna")))})
      a <- capture.output({alignment_matrix = msa(sequences_to_align, method = "ClustalOmega", type = "dna")})
      centre_consensus <- consensus_N(alignment_matrix, round(mean(nchar(sequences_to_align))))
      cumulative_adist_score = cumulative_adist_score + sum(adist(centre_consensus, sequences_to_align))
      number_of_repeats_scored <- number_of_repeats_scored + length(sequences_to_align)
      
      # What is the repeat length variation within the central parts of the array?
      #   Calculation: central 50% repeats length SD normalized by repeat mean length
      mean_centre_width_SD <- c(mean_centre_width_SD, sd(nchar(sequences_to_align)))
      centre_repeats_count <- c(centre_repeats_count, length(sequences_to_align))
    }
    cat("\n")
    if(number_of_repeats_scored != 0) {
      chr_classes$centre_array_edit[k] <- cumulative_adist_score / number_of_repeats_scored ### PREDICTOR ###
      chr_classes$centre_array_width_sd[k] <- sum(mean_centre_width_SD * centre_repeats_count / sum(centre_repeats_count)) ### PREDICTOR ###
    }
    
    if(nrow(chr_family_repeats) > 10) {
      # What is the repeat divergence within the central parts of the chromosome?
      #   Calculation: As above, to identify holocentrics
      central_repeats <- chr_family_repeats[ceiling(nrow(chr_family_repeats) / 4) : floor(3 * nrow(chr_family_repeats) / 4), ]
      repeats_to_align_IDs <- sample(seq_len(nrow(central_repeats)), round(nrow(central_repeats) * desired_fraction_to_align))
      if(length(repeats_to_align_IDs) > max_repeats_to_align) {
        repeats_to_align_IDs <- sample(seq_len(nrow(central_repeats)), max_repeats_to_align)
      } else if(length(repeats_to_align_IDs) < min_repeats_to_align) {
        if(nrow(central_repeats) > min_repeats_to_align) {
          repeats_to_align_IDs <- sample(seq_len(nrow(central_repeats)), min_repeats_to_align)
        } else {
          repeats_to_align_IDs <- seq_len(nrow(central_repeats))
        }
      }
      sequences_to_align <- central_repeats$sequence[repeats_to_align_IDs]
      cat("Chr ", chromosomes[j], ", class ", chr_classes$class[k],", central chromosome repeats (", nrow(chr_family_repeats), " reps)", sep = "")

      # a <- capture.output({alignment_matrix = tolower(as.matrix(msa(sequences_to_align, method = "ClustalOmega", type = "dna")))})
      a <- capture.output({alignment_matrix = msa(sequences_to_align, method = "ClustalOmega", type = "dna")})
      centre_consensus <- consensus_N(alignment_matrix, round(mean(nchar(sequences_to_align))))
      
      chr_classes$centre_chromosome_edit[k] <- sum(adist(centre_consensus, sequences_to_align)) / length(sequences_to_align) ### PREDICTOR ###
      
      # What is the repeat length variation within the central parts of the chromosome?
      #   Calculation: As above, to identify holocentrics
      chr_classes$centre_chromosome_width_sd[k] <- sd(nchar(sequences_to_align)) ### PREDICTOR ###
    }
    
    
    # Can we find the “main” array?
    #   Calculation: What is the SD of array sizes normalised by mean array size
    chr_classes$array_sizes_sd_norm_mean_arr_size[k] <- sd(arrays_sizes) / mean(arrays_sizes) ### PREDICTOR ###
    if(length(arrays_sizes) == 1) {
      chr_classes$array_sizes_sd_norm_mean_arr_size[k] = 0
    }
    chr_classes$array_count[k] <- length(arrays_sizes)
    
    #   What is the TE landscape in the proximity?
    chr_classes$TE_prox_dist[k] <- 100 * abs(mean(chr_family_repeats_t$start_adj) - edta_peak) / chromosome_no_rep_size # the lower the better
    
    chr_classes$TE_prox_SD[k] <- sd(100 * abs(chr_family_repeats_t$start_adj - edta_peak) / chromosome_no_rep_size) # the lower the better
    
    chr_classes$TE_lm_coef[k] <- lm_coef$coefficients[2] # the higher absolute the better
    
    # Score the repeat class based on it’s own mean start position distance to the TE “peak”, score accounting for the correlation calculated before
    
    chr_classes$TE_prox_score[k] <- abs(chr_classes$TE_lm_coef[k]) / (chr_classes$TE_prox_dist[k] + chr_classes$TE_prox_SD[k]) * 100
    
    # What is the gene landscape in the proximity?
    #   Calculation: Identical to the TE, but this time expect negative correlation
    
    chr_classes$gene_prox_dist[k] <- 100 * abs(mean(chr_family_repeats_t$start_adj) - gene_valley) / chromosome_no_rep_size # the lower the better
    
    chr_classes$gene_prox_SD[k] <- sd(100 * abs(chr_family_repeats_t$start_adj - gene_valley) / chromosome_no_rep_size) # the lower the better
    
    chr_classes$gene_lm_coef[k] <- lm_coef_genes$coefficients[2] # the higher the better
    
    chr_classes$gene_prox_score[k] <-  abs(chr_classes$gene_lm_coef[k]) / (chr_classes$gene_prox_dist[k] + chr_classes$gene_prox_SD[k]) * 100
    
    # What TE families can be identified in the proximity?
    #   Calculation: For each individual TE, calculate distance between it and the closest repeat unit of analysed family, then divide the mean distances of predicted centrophilic TEs by mean distances of other TEs, the lower the ratio, the higher the score
    
    sequence_edta$centrophilic_Classification = FALSE
    sequence_edta$centrophilic_Classification[grep(chr_classes$pred_centrophilic_TE[1], sequence_edta$Classification)] = TRUE
    
    sequence_edta$dist_to_closest_rep <- 0
    
    for(seq_edta_id in seq_len(nrow(sequence_edta))) {
      sequence_edta$dist_to_closest_rep[seq_edta_id] <- min(abs(sequence_edta$V4[seq_edta_id] - chr_family_repeats_t$start))
    }
    # plot(sequence_edta$V4, sequence_edta$dist_to_closest_rep)
    # abline(v = sequence_edta$V4[sequence_edta$centrophilic_Classification == TRUE])
    # mean(sequence_edta$dist_to_closest_rep[sequence_edta$centrophilic_Classification == TRUE])
    # mean(sequence_edta$dist_to_closest_rep[sequence_edta$centrophilic_Classification == FALSE])
    # 
    # ggplot(sequence_edta, aes(x=centrophilic_Classification, y=dist_to_closest_rep)) + 
    #   geom_violin() + geom_boxplot(width=0.1)
    # 
    chr_classes$t_test_p_val[k] <- -1
    chr_classes$t_test_t_val[k] <- -1
    if(length(sequence_edta$dist_to_closest_rep[sequence_edta$centrophilic_Classification]) < 6) next
    if(length(unique(sequence_edta$dist_to_closest_rep[sequence_edta$centrophilic_Classification])) == 1) next
    chr_classes$t_test_p_val[k] <-  t.test(sequence_edta$dist_to_closest_rep[sequence_edta$centrophilic_Classification], 
                                           sequence_edta$dist_to_closest_rep[!sequence_edta$centrophilic_Classification])$p.value
    
    chr_classes$t_test_t_val[k] <- t.test(sequence_edta$dist_to_closest_rep[sequence_edta$centrophilic_Classification], 
                                          sequence_edta$dist_to_closest_rep[!sequence_edta$centrophilic_Classification])$statistic
    
  }
  
  
  ### Analyse predictor values =================================================
  ### Centromeric repeat
  
  chr_classes$score_total <- 0
  # total_bp_norm_chr, 1 at 10%
  values <- chr_classes$total_bp_norm_chr
  values[values > 0.1] = 0.1
  chr_classes$score_total <- chr_classes$score_total + (values / 0.1)
  
  # total_bp_norm_rep, 1 at 50%
  values <- chr_classes$total_bp_norm_rep
  values[values > 0.5] = 0.5
  chr_classes$score_total <- chr_classes$score_total + (values / 0.5)
  
  # start_sd_norm_chr, 0 at 10%
  values <- chr_classes$start_sd_norm_chr
  values[values > 0.1] = 0.1
  chr_classes$score_total <- chr_classes$score_total + 1-(values / 0.1)
  
  # start_norm_chr_0_50, 1 at 0%, 0 at 25% and 1 at 50%
  values <- chr_classes$start_norm_chr_0_50
  values[values > 0.25] = abs(values[values > 0.25] - 0.5)
  chr_classes$score_total <- chr_classes$score_total + 1-(values / 0.25)
  
  # gaps_with_TEs_fraction, 1 at 75%
  values <- chr_classes$gaps_with_TEs_fraction
  values[values == -1] = 0
  values[values > 0.75] = 0.75
  chr_classes$score_total <- chr_classes$score_total + (values / 0.75)
  
  # centre_array_edit, 0 at 15%
  values <- chr_classes$centre_array_edit
  values[values == -1] = 15
  values[values > 15] = 15
  chr_classes$score_total <- chr_classes$score_total + 1 - (values / 15)
  
  # centre_array_width_sd, 0 at 15%
  values <- chr_classes$centre_array_width_sd
  values[values == -1] = 15
  values[values > 15] = 15
  chr_classes$score_total <- chr_classes$score_total + 1 - (values / 15)
  
  # centre_chromosome_edit, 0 at 15%
  values <- chr_classes$centre_chromosome_edit
  values[values == -1] = 15
  values[values > 15] = 15
  chr_classes$score_total <- chr_classes$score_total + 1 - (values / 15)
  
  # centre_chromosome_width_sd, 0 at 15%
  values <- chr_classes$centre_chromosome_width_sd
  values[values == -1] = 15
  values[values > 15] = 15
  chr_classes$score_total <- chr_classes$score_total + 1 - (values / 15)
  
  # TE_prox_score, 1 at 100
  values <- chr_classes$TE_prox_score
  values[values == -1] = 0
  values[values > 100] = 100
  chr_classes$score_total <- chr_classes$score_total + (values / 100)
  
  # gene_prox_score, 1 at 100
  values <- chr_classes$gene_prox_score
  values[values == -1] = 0
  values[values > 100] = 100
  chr_classes$score_total <- chr_classes$score_total + (values / 100)
  
  # t_test_p_val, 0 at 0.001; t_test_t_val must be negative
  values <- chr_classes$t_test_p_val
  values[values == -1] = 0.001
  values[values > 0.001] = 0.001
  values[chr_classes$t_test_t_val > 0] = 0.001
  chr_classes$score_total <- chr_classes$score_total + 1 - (values / 0.001)
  
  
  
  ### Analyse predictor values =================================================
  ### Holo/monocentricity
  chr_classes$holo_signatures <- 0
  
  # start_sd_norm_chr, 3 at 30%
  values <- chr_classes$start_sd_norm_chr
  values[values > 0.3] = 0.3
  chr_classes$holo_signatures <- chr_classes$holo_signatures + (values / 0.3)
  
  # start_norm_chr_0_50, 1 at 50%
  values <- chr_classes$start_norm_chr_0_50
  values[values > 0.5] = 0.5
  chr_classes$holo_signatures <- chr_classes$holo_signatures + (values / 0.5)
  
  # (centre_chromosome_edit / centre_array_edit), 2 at 2, 0 at 1
  values <- chr_classes$centre_chromosome_edit / chr_classes$centre_array_edit
  values[chr_classes$centre_chromosome_edit == -1] = 0
  values = values - 1
  values[values < 0] = 0
  values[values > 1] = 1
  chr_classes$holo_signatures <- chr_classes$holo_signatures + 2 * (values / 1)
  
  # array_count, 10 at 100
  values <- chr_classes$array_count
  values[values > 100] = 100
  chr_classes$holo_signatures <- chr_classes$holo_signatures + (values / 10)
  
  # array_sizes_sd_norm_mean_arr_size, 1 at 0, 0 at 1.5, multiplied by array_count - 5 (0 if negative)
  values <- chr_classes$array_sizes_sd_norm_mean_arr_size
  values[values > 1.5] = 1.5
  values <- (1.5 - values) * (chr_classes$array_count - 5)
  values[values < 0] = 0
  chr_classes$holo_signatures <- chr_classes$holo_signatures + (values / 0.5)
  
  # TE_lm_coef, 0 at -1, 3 at 0, use power to make logarithmic scoring
  values <- chr_classes$TE_lm_coef
  values[values < -1] <- -1
  chr_classes$holo_signatures <- chr_classes$holo_signatures + 3 * values^5
  
  # gene_lm_coef, 0 at 1, 3 at 0, use power to make logarithmic scoring
  values <- chr_classes$gene_lm_coef
  values[values > 1] <- 1
  chr_classes$holo_signatures <- chr_classes$holo_signatures + 3 * (1-values)^5
  
  
  ### Merge into big file
  chr_classes$chromosome <- chromosomes[j]
  genome_classes <- rbind(genome_classes, chr_classes)
}
dev.off()

write.csv(x = genome_classes, file = paste0("/home/pwlodzimierz/ToL/centromeric_class_identification/genome_classes_", assembly_name, ".csv"))
write.csv(x = genome_classes, file = paste0("genome_classes_", assembly_name, ".csv"))

print("Finished successfully")




