taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 193
print(i)

suppressMessages(library(seqinr))
library(ggplot2)
library(pheatmap)


samples = 10000 # PCA


### load data
data_directories <- list.dirs(path = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs", recursive = FALSE, full.names = TRUE)

setwd(data_directories[i])
repeats_file <- list.files(path = ".", pattern = "repeats_filtered.csv", full.names = TRUE)
if(length(repeats_file) != 1) {cat("issue with repeats files: \n", repeats_file); stop()}

edta_file <- list.files(path = ".", pattern = "_edta_filtered.csv.reassigned", full.names = TRUE)
if(length(edta_file) != 1) {cat("issue with edta files: \n", edta_file); stop()}

arrays_file <- list.files(path = ".", pattern = "_centromeric_arrays.csv", full.names = TRUE)
if(length(arrays_file) != 1) {cat("issue with gaps files: \n", arrays_file); stop()}

gaps_file <- list.files(path = ".", pattern = "_centromeric_gaps.csv", full.names = TRUE)
if(length(gaps_file) != 1) {cat("issue with gaps files: \n", gaps_file); stop()}

repeats <- read.csv(repeats_file)
edta <- read.csv(edta_file)
arrays <- read.csv(arrays_file)
gaps <- read.csv(gaps_file)

satellites_metadata <- read.csv(file = "/home/pwlodzimierz/ToL/curated_satellites_repDec24_jan2025.csv")
chr_sizes_metadata <- read.csv(file = "/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full.Ian.csv")

### filter the data

assembly_name_long <- strsplit(data_directories[i], split = "/")[[1]][7]
assembly_name_short <- strsplit(assembly_name_long, split = ".fa")[[1]][1]

satellites_metadata <- satellites_metadata[satellites_metadata$Genome == assembly_name_long, ]


chr_sizes_metadata <- chr_sizes_metadata[chr_sizes_metadata$assembly.name == assembly_name_long, ]
chr_sizes_metadata <- chr_sizes_metadata[chr_sizes_metadata$is.chr == 1, ]


repeats <- repeats[repeats$seqID %in% chr_sizes_metadata$chromosome.name, ]

###

for(h in seq_len(nrow(satellites_metadata))) {
  satellite_TRASH_names <- strsplit(satellites_metadata$TRASH_name_dec2024runs[h], split = ";")[[1]]
  repeat_trash <- satellites_metadata$Satellite_name[h]
  
  sample_repeats <- repeats[repeats$new_class %in% satellite_TRASH_names, ]
  repeat_sample <- sort(sample(1 : nrow(sample_repeats), samples, replace = T))
  sample_repeats <- sample_repeats[repeat_sample, ]
  
  # Step 1: Align repeats
  write.fasta(sequences = as.list(sample_repeats$sequence), names = c(paste(seq_len(nrow(sample_repeats)))), file.out = paste0(assembly_name_short, "_", repeat_trash, "_repeats.fasta"))
  
  system(paste0("mafft --retree 2 --thread -1 --inputorder ", assembly_name_short, "_", repeat_trash,  "_repeats.fasta > ", assembly_name_short, "_", repeat_trash, "_aligned.fasta"))
  
  fasta_file <- paste0(assembly_name_short, "_", repeat_trash, "_aligned.fasta")
  sequences <- read.fasta(fasta_file, seqtype = "DNA", as.string = TRUE)
  
  # Convert to character vector
  seq_vector <- sapply(sequences, function(x) toupper(as.character(x)))
  
  # Check sequence length (should all be equal after MAFFT alignment)
  seq_lengths <- nchar(seq_vector)
  summary(seq_lengths)
  if (length(unique(seq_lengths)) > 1) {
    stop("Sequences have varying lengths after alignment. Ensure MAFFT output is properly aligned.")
  }
  seq_length <- unique(seq_lengths)
  cat("Aligned sequence length:", seq_length, "bp\n")
  
  # Step 2: Feature Extraction
  # Toggle between k-mer and positional features
  feature_type <- "positional"  # Change to "kmer" to use k-mer frequencies
  
  if (feature_type == "kmer") {
    # K-mer frequencies (4-mers)
    k <- 4
    kmers <- oligonucleotideFrequency(DNAStringSet(seq_vector), width = k, as.prob = FALSE)
    feature_matrix <- as.matrix(kmers)  # 1M rows x 256 columns (4^4)
    cat("Using 4-mer frequencies with", ncol(feature_matrix), "features\n")
  } else if (feature_type == "positional") {
    # Positional nucleotide frequencies
    seq_matrix <- do.call(rbind, strsplit(seq_vector, ""))
    feature_matrix <- matrix(0, nrow = length(seq_vector), ncol = seq_length * 4)
    for (i in 1:seq_length) {
      # if(i%/%100 == 0) cat(i,"\n")
      idx <- (i - 1) * 4 + 1:4
      feature_matrix[, idx[1]] <- seq_matrix[, i] == "A"
      feature_matrix[, idx[2]] <- seq_matrix[, i] == "C"
      feature_matrix[, idx[3]] <- seq_matrix[, i] == "G"
      feature_matrix[, idx[4]] <- seq_matrix[, i] == "T"
      # Gaps (-) are left as 0,0,0,0
    }
    cat("Using positional frequencies with", ncol(feature_matrix), "features\n")
  } else {
    stop("Invalid feature_type. Use 'kmer' or 'positional'.")
  }
  
  # Remove constant/zero-variance columns
  variance <- apply(feature_matrix, 2, var)
  constant_cols <- which(variance == 0)
  if (length(constant_cols) > 0) {
    cat("Removing", length(constant_cols), "constant/zero-variance columns from feature_matrix\n")
    feature_matrix <- feature_matrix[, variance > 0]
  }
  
  # Check if any columns remain
  if (ncol(feature_matrix) == 0) {
    stop("All columns have zero variance. Try a different feature extraction method.")
  }
  
  
  
  sample_idx <- sort(sample(nrow(feature_matrix), 500))
  downsampled_matrix <- feature_matrix[sample_idx, ]
  
  group_info <- sample_repeats$seqID[sample_idx]
  
  # Create similarity matrix with explicit row names
  sim_matrix <- cor(t(downsampled_matrix))
  rownames(sim_matrix) <- paste0("Seq_", 1:500)  # Simple unique names
  colnames(sim_matrix) <- rownames(sim_matrix)
  
  # Create annotation with matching row names
  annotation_row <- data.frame(Group = group_info, row.names = rownames(sim_matrix))
  
  # Define colors for unique groups
  unique_groups <- unique(group_info)
  group_colors <- rainbow(length(unique_groups))
  names(group_colors) <- unique_groups
  annotation_colors <- list(Group = group_colors)
  
  
  pdf(file = paste0( assembly_name_short, "_pheatmap_", repeat_trash, "_repeats.pdf"), width = 12, height = 12, onefile = T)
  # Create heatmap with group-based coloring
  pheatmap(sim_matrix,
           clustering_method = "complete",
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           annotation_row = annotation_row,
           annotation_colors = annotation_colors,
           main = "Repeat Sequence Similarities (500 sequences)")
  pheatmap(sim_matrix,
           cluster_rows = FALSE, 
           cluster_cols = FALSE,
           annotation_row = annotation_row,
           annotation_colors = annotation_colors,
           main = "Repeat Sequence Similarities (500 sequences)")
  
  dev.off()
  
  
  pdf(file = paste0("/home/pwlodzimierz/ToL/upload_files/20_PCA/", assembly_name_short, "_pheatmap_", repeat_trash, "_repeats.pdf"), width = 12, height = 12, onefile = T)
  # Create heatmap with group-based coloring
  pheatmap(sim_matrix,
           clustering_method = "complete",
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           annotation_row = annotation_row,
           annotation_colors = annotation_colors,
           main = "Repeat Sequence Similarities (500 sequences)")
  pheatmap(sim_matrix,
           cluster_rows = FALSE, 
           cluster_cols = FALSE,
           annotation_row = annotation_row,
           annotation_colors = annotation_colors,
           main = "Repeat Sequence Similarities (500 sequences)")
  
  dev.off()
}
