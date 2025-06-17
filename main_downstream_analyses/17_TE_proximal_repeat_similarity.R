#!/usr/bin/env Rscript
taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 193
print(i)

# g. make a matrix of gap similarity scores, make a tree with gap lengths and TE annotation ploted alongside
# --> 17_TE_...R matrix csv 1
# h. add to the values calculated in g: merge up and downstream repeats and make matrix of those repeat regions similarity (direct alignment). 
#    Plot a scatter of gap (g) vs repeat similarity scores. 
# --> 17_TE_...R matrix csv 2
# i. measure similarity between 1 000 randomly picked repeat region pairs using the method as in g. plot the histogram of those similarities, marking 
#    where the values from (g) fall along the distribution
# --> 17_TE_...R two-column df csv

.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))


suppressMessages(library(seqinr))
suppressMessages(library(Biostrings))
suppressMessages(library(text2vec))
suppressMessages(library(dplyr))
suppressMessages(library(Matrix))
suppressMessages(library(pheatmap))


colors_vector <- rep(c("#E69F00", "#0072B2", "#A6D854", "#E7298A", "#F0E442", "#7570B3", "#66C2A5" ,
                       "#D55E00", "#56B4E9", "#CC79A7", "#FFD92F" ,"#8DA0CB", "#009E73", "#FC8D62", 
                       "#B3B3CC", "#66A61E", "#E5C494"), 100)

surrounding_sequence_length <- 50
k <- 6 
min_gap_size <- 750
min_gap_size_2 <- 250

get_kmer_freq <- function(sequence, k) {
  seq <- DNAString(as.character(sequence))
  kmer_counts <- oligonucleotideFrequency(seq, width = k, as.prob = FALSE)
  kmer_freq <- kmer_counts / sum(kmer_counts)
  return(kmer_freq)
}

compute_similarity_matrix <- function(sequences, k, metric = "cosine") {
  # Generate k-mer frequency matrix for all sequences
  kmer_freq_matrix <- t(sapply(sequences, get_kmer_freq, k = k))
  
  # Compute pairwise similarities
  if (metric == "cosine") {
    sim_matrix <- as.matrix(kmer_freq_matrix %*% t(kmer_freq_matrix) /
                              (sqrt(rowSums(kmer_freq_matrix^2)) %*%
                                 t(sqrt(rowSums(kmer_freq_matrix^2)))))
  } else if (metric == "euclidean") {
    dist_matrix <- as.matrix(dist(kmer_freq_matrix, method = "euclidean"))
    sim_matrix <- 1 / (1 + dist_matrix)
  } else {
    stop("Unsupported metric. Use 'cosine' or 'euclidean'.")
  }
  
  rownames(sim_matrix) <- names(sequences)
  colnames(sim_matrix) <- names(sequences)
  return(sim_matrix)
}


### load data
data_directories <- list.dirs(path = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs", recursive = FALSE, full.names = TRUE)
# 14, 16, 19
assembly_name_long <- strsplit(data_directories[i], split = "/")[[1]][7]
assembly_name_short <- strsplit(assembly_name_long, split = ".fa")[[1]][1]

cat(data_directories[i], "\n")

assembly_files <- list.files(path = "/home/pwlodzimierz/ToL/Assemblies/fastas_2021_Michael", recursive = FALSE, full.names = TRUE)
assembly_files <- assembly_files[!grepl(".fai", assembly_files)]

assembly_file <- assembly_files[grep(assembly_name_short, assembly_files)]

setwd(data_directories[i])
repeats_file <- list.files(path = ".", pattern = "repeats_filtered.csv", full.names = TRUE)
if(length(repeats_file) != 1) {cat("issue with repeats files: \n", assembly_name_short, repeats_file); stop()}

edta_file <- list.files(path = ".", pattern = "_edta_filtered.csv.reassigned", full.names = TRUE)
if(length(edta_file) != 1) {cat("issue with edta files: \n", assembly_name_short, edta_file); stop()}

arrays_file <- list.files(path = ".", pattern = "_centromeric_arrays.csv", full.names = TRUE)
if(length(arrays_file) != 1) {cat("issue with arrays files: \n", assembly_name_short, arrays_file); stop()}

gaps_file <- list.files(path = ".", pattern = "_centromeric_gaps.csv", full.names = TRUE)
if(length(gaps_file) != 1) {cat("issue with gaps files: \n", assembly_name_short, gaps_file); stop()}

repeats <- read.csv(repeats_file); cat("repeats in\n")
edta <- read.csv(edta_file); cat("edta in\n")
arrays <- read.csv(arrays_file); cat("arrays in\n")
gaps <- read.csv(gaps_file); cat("gaps in\n")
fasta <- read.fasta(assembly_file); cat("fasta in\n")

satellites_metadata <- read.csv(file = "/home/pwlodzimierz/ToL/curated_satellites_repDec24_jan2025.csv")
chr_sizes_metadata <- read.csv(file = "/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full.Ian.csv")


if(assembly_name_short == "daMisOron1.1") {
  gaps_of_interest <- read.table("/home/pwlodzimierz/ToL/temp_data/daMisOron1_mosaic_gap_grouped.txt", sep = "\t", header = TRUE)
  gaps <- gaps[which(gaps$start %in% gaps_of_interest$start), ]
} else if(assembly_name_short == "daLinVulg1.1") {
  gaps_of_interest <- read.table("/home/pwlodzimierz/ToL/temp_data/daLinVulg1_mosaic_gap_grouped.txt", sep = "\t", header = TRUE)
  gaps <- gaps[which(gaps$start %in% gaps_of_interest$start), ]
} else if(assembly_name_short == "daSheArve1.1") {
  gaps_of_interest <- read.table("/home/pwlodzimierz/ToL/temp_data/daSheArve1_mosaic_gap_grouped.txt", sep = "\t", header = TRUE)
  gaps <- gaps[which(gaps$start %in% gaps_of_interest$start), ]
} else stop()

### filter the data

satellites_metadata <- satellites_metadata[satellites_metadata$Genome == assembly_name_long, ]

satellite_TRASH_names <- strsplit(satellites_metadata$TRASH_name_dec2024runs[1], split = ";")[[1]]

chr_sizes_metadata <- chr_sizes_metadata[chr_sizes_metadata$assembly.name == assembly_name_long, ]
chr_sizes_metadata <- chr_sizes_metadata[chr_sizes_metadata$is.chr == 1, ]

repeats <- repeats[repeats$new_class %in% satellite_TRASH_names, ]

# samples <- sample(1 : nrow(repeats), 100)
# write.fasta(sequences = as.list(repeats$sequence[samples]), 
#             names = paste0(rep(satellite_TRASH_names,100), rep("_",100), repeats$seqID[samples], rep("_",100), repeats$start[samples], rep(":",100), repeats$end[samples]), 
#             file.out = paste0("/home/pwlodzimierz/ToL/upload_files/sample_repeats_for_Alex/", assembly_name_short, "_sample_repeats.fasta"))

###
gaps_seq <- NULL
gap_surrounding_seq <- NULL
chr_name <- NULL
annotation_colors <- NULL

gaps <- gaps[order(gaps$chromosome, decreasing = FALSE), ]
gaps$strand <- rep("+", nrow(gaps))


for(j in seq_len(nrow(gaps))) {
  cat(j, "/", nrow(gaps), "\n")
  
  if((gaps$end[j] - gaps$start[j] + 1) < min_gap_size_2) {
    gaps$strand[j] <- ""
    next
  }
  
  repeats_around_gaps <- NULL
  repeats_around_gaps <- c(repeats$strand[repeats$seqID == gaps$chromosome[j] & repeats$start >= (gaps$start[j] - 10000) & repeats$start <= gaps$start[j]], 
                           repeats$strand[repeats$seqID == gaps$chromosome[j] & repeats$start >= gaps$end[j] & repeats$start <= gaps$end[j] + 10000])
  if(length(repeats_around_gaps) != 0) if(sum(repeats_around_gaps == "+") / length(repeats_around_gaps) < 0.5) {
    gaps$strand[j] = "-"
  }
  
  if((gaps$end[j] - gaps$start[j] + 1) < min_gap_size) next
  
  chr_len <- length(fasta[[which(names(fasta) == gaps$chromosome[j])]])

  seq_gap <- fasta[[which(names(fasta) == gaps$chromosome[j])]][gaps$start[j] : gaps$end[j]]
  
  ext_start <- gaps$start[j] - surrounding_sequence_length
  if(ext_start <= 0) ext_start <- 1
  ext_end <- gaps$start[j] - 1
  if(ext_end <= 0) ext_end <- 1
  seq_up <- fasta[[which(names(fasta) == gaps$chromosome[j])]][ext_start : ext_end]
  
  ext_start <- gaps$end[j] + 1
  if(ext_start >= chr_len) ext_start <- chr_len
  ext_end <- gaps$end[j] + (surrounding_sequence_length)
  if(ext_end >= chr_len) ext_end <- chr_len
  seq_down <- fasta[[which(names(fasta) == gaps$chromosome[j])]][ext_start : ext_end]

  seq_round <- c(seq_up, seq_down)

  seq_round = seq_round[seq_round %in% c("a","c","t","g")]
  seq_gap = seq_gap[seq_gap %in% c("a","c","t","g")]
  
  if(gaps$strand[j] == "-") {
    seq_gap <- rev(comp(seq_gap)) # if most of the repeats around are on - strand, reverse the sequence
    seq_round <- rev(comp(seq_round)) # if most of the repeats around are on - strand, reverse the gap
  }
  
  
  gaps_seq[length(gaps_seq) + 1] <- list(seq_gap)
  gap_surrounding_seq[length(gap_surrounding_seq) + 1] <- list(seq_round)
  chr_name <- c(chr_name, gaps$chromosome[j])
  annotation_colors <- c(annotation_colors, colors_vector[which(chr_sizes_metadata$chromosome.name == gaps$chromosome[j])])
  
}

gaps$width <- gaps$end - gaps$start + 1

write.csv(x = gaps, file = paste0(assembly_name_short, "_centromeric_gaps_with_strand.csv"), row.names = FALSE)
write.csv(x = gaps, file = paste0("/home/pwlodzimierz/ToL/upload_files/17_gap_files_with_strand/", assembly_name_short, "_centromeric_gaps_with_strand.csv"), row.names = FALSE)

### similarity of surrounding sequence

sequences <- DNAStringSet(sapply(gap_surrounding_seq, paste, collapse = ""))
names(sequences) <- paste0(chr_name, "_rep_", 1:length(gap_surrounding_seq))
num_sequences <- length(sequences)
annotation_df <- data.frame(
  row.names = names(sequences),
  Group = chr_name
)
annotation_colors <- colors_vector[1 : length(unique(chr_name))]
names(annotation_colors) <- unique(chr_name)
annotation_colors <- list(Group = annotation_colors)


sim_matrix_surrounding_seq <- compute_similarity_matrix(sequences, k = k, metric = "cosine")
sim_matrix_surrounding_seq[is.na(sim_matrix_surrounding_seq)] = 0
rownames(annotation_df) <- rownames(sim_matrix_surrounding_seq)

write.csv(x = sim_matrix_surrounding_seq, file = paste0(assembly_name_short, "_gaps_surrounding_seq_", surrounding_sequence_length, "_bp_", k, "_kmer.csv"))

pdf(file = paste0(assembly_name_short, "_gaps_surrounding_seq_", surrounding_sequence_length, "_bp_", k, "_kmer.pdf"), width = 12, height = 12, onefile = T)
pheatmap(sim_matrix_surrounding_seq, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", 
         annotation_colors = annotation_colors,
         annotation_row = annotation_df,
         annotation_col = annotation_df,  
         fontsize_row = 6,  
         fontsize_col = 6,
         main = paste0(assembly_name_short, " K-mer (k=", k, ") Similarity Matrix repeats around gaps"))
pheatmap(sim_matrix_surrounding_seq, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         annotation_colors = annotation_colors,
         annotation_row = annotation_df,
         annotation_col = annotation_df,  
         fontsize_row = 6,  
         fontsize_col = 6,
         main = paste0(assembly_name_short, " K-mer (k=", k, ") Similarity Matrix repeats around gaps"))
dev.off()


file.copy(from = paste0(assembly_name_short, "_gaps_surrounding_seq_", surrounding_sequence_length, "_bp_", k, "_kmer.pdf"), 
          to = paste0("/home/pwlodzimierz/ToL/upload_files/17_matrices_simialrity_gaps/", assembly_name_short, "_gaps_surrounding_seq_", surrounding_sequence_length, "_bp_", k, "_kmer.pdf"), overwrite = T)



### similarity of gaps

sequences <- DNAStringSet(sapply(gaps_seq, paste, collapse = ""))
names(sequences) <- paste0(chr_name, "_gap_", 1:length(gap_surrounding_seq))
num_sequences <- length(sequences)
annotation_df <- data.frame(
  row.names = names(sequences),
  Group = chr_name
)
annotation_colors <- colors_vector[1 : length(unique(chr_name))]
names(annotation_colors) <- unique(chr_name)
annotation_colors <- list(Group = annotation_colors)



sim_matrix_gaps <- compute_similarity_matrix(sequences, k = k, metric = "cosine")
sim_matrix_gaps[is.na(sim_matrix_gaps)] = 0

write.csv(x = sim_matrix_gaps, file = paste0(assembly_name_short, "_gaps_seq_", k, "_kmer.csv"))

pdf(file = paste0(assembly_name_short, "_gaps_seq_", k, "_kmer.pdf"), width = 12, height = 12, onefile = T)
pheatmap(sim_matrix_gaps, clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", 
         annotation_colors = annotation_colors,
         annotation_row = annotation_df,
         annotation_col = annotation_df,  
         fontsize_row = 6,  
         fontsize_col = 6,
         main = paste0(assembly_name_short, " K-mer (k=", k, ") Similarity Matrix gaps"))
pheatmap(sim_matrix_gaps, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         annotation_colors = annotation_colors,
         annotation_row = annotation_df,
         annotation_col = annotation_df,  
         fontsize_row = 6,  
         fontsize_col = 6,
         main = paste0(assembly_name_short, " K-mer (k=", k, ") Similarity Matrix gaps"))
dev.off()

file.copy(from = paste0(assembly_name_short, "_gaps_seq_", k, "_kmer.pdf"), to = paste0("/home/pwlodzimierz/ToL/upload_files/17_matrices_simialrity_gaps/", assembly_name_short, "_gaps_seq_", k, "_kmer.pdf"), overwrite = T)


# 
# pdf(file = paste0(assembly_name_short, "_scatter_gaps_vs_repeats_around_", k, "_kmer.pdf"), width = 12, height = 12, onefile = T)
# plot(y = sim_matrix_surrounding_seq, x = sim_matrix_gaps, pch = 16, col = "#00000060", 
#      main = paste0(assembly_name_short, " K-mer (k=", k, ") scatter gaps vs repeats around"))
# dev.off()
# pdf(file = paste0("/home/pwlodzimierz/ToL/upload_files/17_matrices_simialrity_gaps/", assembly_name_short, "_scatter_gaps_vs_repeats_around_", k, "_kmer.pdf"), width = 12, height = 12)
# plot(y = sim_matrix_surrounding_seq, x = sim_matrix_gaps, pch = 16, col = "#00000060", 
#      main = paste0(assembly_name_short, " K-mer (k=", k, ") scatter gaps vs repeats around"))
# dev.off()
# 
# pdf(file = paste0(assembly_name_short, "_scatter_gaps_vs_repeats_around_over_80perc_", k, "_kmer.pdf"), width = 12, height = 12, onefile = T)
# plot(y = sim_matrix_surrounding_seq[sim_matrix_gaps >= 0.8], x = sim_matrix_gaps[sim_matrix_gaps >= 0.8], pch = 16, col = "#00000060", 
#      main = paste0(assembly_name_short, " K-mer (k=", k, ") scatter gaps vs repeats around"))
# dev.off()
# pdf(file = paste0("/home/pwlodzimierz/ToL/upload_files/17_matrices_simialrity_gaps/", assembly_name_short, "_scatter_gaps_vs_repeats_around_over_80perc_", k, "_kmer.pdf"), width = 12, height = 12)
# plot(y = sim_matrix_surrounding_seq[sim_matrix_gaps >= 0.8], x = sim_matrix_gaps[sim_matrix_gaps >= 0.8], pch = 16, col = "#00000060", 
#      main = paste0(assembly_name_short, " K-mer (k=", k, ") scatter gaps vs repeats around"))
# dev.off()

library(ggplot2)

for(j in nrow(sim_matrix_gaps)) {
  sim_matrix_gaps[j,j] = NA
}
for(j in nrow(sim_matrix_surrounding_seq)) {
  sim_matrix_surrounding_seq[j,j] = NA
}

# Flatten matrices into vectors
gaps_vec <- as.vector(sim_matrix_gaps)
surrounding_vec <- as.vector(sim_matrix_surrounding_seq)

gaps_vec[gaps_vec > 1] = 1
gaps_vec[gaps_vec < 0] = 0
surrounding_vec[surrounding_vec > 1] = 1
surrounding_vec[surrounding_vec < 0] = 0
 
# Check for NA values and remove them
valid_idx <- !is.na(gaps_vec) & !is.na(surrounding_vec)
gaps_vec <- gaps_vec[valid_idx]
surrounding_vec <- surrounding_vec[valid_idx]

# Define bin edges and labels
bins_lower <- seq(0, 0.9, by = 0.1)  # 0, 0.1, ..., 0.9
bins_upper <- seq(0.9, 1.01, by = 0.01)  # 0.9, 0.91, ..., 1, 1.01
bins <- c(bins_lower[-length(bins_lower)], bins_upper)  # Combine, avoiding duplicate 0.9

# Define bin labels
bin_labels <- paste0("<", bins[-length(bins)], ",", bins[-1], ")")
bin_labels[length(bin_labels)] <- "<1,1.01]"  # Last bin includes up to 1.01

# Bin the sim_matrix_gaps values, ensuring 1 is included
gap_bins <- cut(gaps_vec, breaks = bins, include.lowest = TRUE, right = TRUE, labels = bin_labels)

# Create a data frame for ggplot2
plot_data <- data.frame(
  Gaps_Bin = gap_bins,
  Surrounding_Seq = surrounding_vec
)

plot_data <- plot_data[!is.na(plot_data$Gaps_Bin),]

# Create the violin plot
p <- ggplot(plot_data, aes(x = Gaps_Bin, y = Surrounding_Seq)) +
  geom_violin(fill = "skyblue", color = "black") +
  geom_boxplot(width = 0.2, fill = "white") +
  labs(
    x = "sim_matrix_gaps Bins",
    y = "sim_matrix_surrounding_seq",
    title = paste0(assembly_name_short, " Violin Plot of Surrounding Sequence Similarity by Gap Similarity Bins")
  ) +
  theme_minimal() +
  ylim(c(0,1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability

ggsave(filename = paste0(assembly_name_short, "volin_box_gaps_against_surrounding_seq.pdf"), dpi = 320)

ggsave(filename = paste0("/home/pwlodzimierz/ToL/upload_files/17_matrices_simialrity_gaps/", assembly_name_short, "volin_box_gaps_against_surrounding_seq.pdf"), dpi = 320)






























