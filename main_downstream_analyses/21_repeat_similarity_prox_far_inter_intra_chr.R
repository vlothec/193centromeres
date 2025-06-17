#!/usr/bin/env Rscript
.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))

# c. sample 10 000 pairs of repeats, calculate their similarity, plot simialrity vs physical distance (within chromosomes)
# --> 21_repeat_similarity_prox_far_inter_intra_chr.R two-column df csv 1 
# d. sample 10 000 repeats, calculate the similarity with:
#     10 000 repeats in a 100-repeat distance 
#     10 000 repeats in a > 1000-repeat distance 
#     10 000 repeats on different chromosomes
#    and make a 3-bar bar-plot for those (try to keep the y axis on the same range as the c plot for easier comparison)
# --> 21_repeat_similarity_prox_far_inter_intra_chr.R two-column df csv 2 3 and 4
taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 193
print(i)

suppressMessages(library(seqinr))
library(ggplot2)
library(dplyr)
library(cowplot)

### settings

sample_within_chromosomes <- 10000
sample_between_chromosomes <- 10000



### load data
data_directories <- list.dirs(path = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs", recursive = FALSE, full.names = TRUE)

setwd(data_directories[i])
repeats_file <- list.files(path = ".", pattern = "repeats_filtered.csv", full.names = TRUE)
if(length(repeats_file) != 1) {cat("issue with repeats files: \n", repeats_file); stop()}

repeats <- read.csv(repeats_file)

satellites_metadata <- read.csv(file = "/home/pwlodzimierz/ToL/curated_satellites_repDec24_jan2025.csv")
chr_sizes_metadata <- read.csv(file = "/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full.Ian.csv")

### filter the data

assembly_name_long <- strsplit(data_directories[i], split = "/")[[1]][7]
assembly_name_short <- strsplit(assembly_name_long, split = ".fa")[[1]][1]

satellites_metadata <- satellites_metadata[satellites_metadata$Genome == assembly_name_long, ]

satellite_TRASH_names <- NULL
for(j in seq_len(nrow(satellites_metadata))) {
  satellite_TRASH_names <- c(satellite_TRASH_names, strsplit(satellites_metadata$TRASH_name_dec2024runs[j], split = ";")[[1]])
}

if(length(satellite_TRASH_names) == 0) stop()

chr_sizes_metadata <- chr_sizes_metadata[chr_sizes_metadata$assembly.name == assembly_name_long, ]
chr_sizes_metadata <- chr_sizes_metadata[chr_sizes_metadata$is.chr == 1,]

repeats <- repeats[repeats$new_class %in% satellite_TRASH_names, ]


###

# TODO: account for multiple classes per species of satellite repeats

### get similarities within and between  chromosomes

# Sample within chromosomes
sample_within <- function(repeats, n) {
  # Initialize vectors for pairs
  idx1 <- integer(n)
  idx2 <- integer(n)
  
  # Get chromosome frequencies
  chroms <- unique(repeats$seqID)
  
  # Sample n pairs
  for(i in 1:n) {
    # Randomly select a chromosome
    chr <- sample(chroms, 1)
    
    # Get indices for this chromosome
    chr_idx <- which(repeats$seqID == chr)
    
    # Sample two different indices from same chromosome
    pair <- sample(chr_idx, 2, replace = FALSE)
    
    idx1[i] <- pair[1]
    idx2[i] <- pair[2]
  }
  
  return(list(idx1 = idx1, idx2 = idx2))
}

# Sample between chromosomes
sample_between <- function(repeats, n) {
  # Initialize vectors for pairs
  idx1 <- integer(n)
  idx2 <- integer(n)
  
  for(i in 1:n) {
    # Sample first index
    idx1[i] <- sample(1:nrow(repeats), 1)
    
    # Get chromosome of first index
    chr1 <- repeats$seqID[idx1[i]]
    
    # Get indices from different chromosomes
    diff_chr_idx <- which(repeats$seqID != chr1)
    
    # Sample second index from different chromosome
    idx2[i] <- sample(diff_chr_idx, 1)
  }
  
  return(list(idx1 = idx1, idx2 = idx2))
}

# Generate samples
within_pairs <- sample_within(repeats, sample_within_chromosomes)
between_pairs <- sample_between(repeats, sample_between_chromosomes)

within_pairs_df <- data.frame(
  idx1 = within_pairs$idx1,
  idx2 = within_pairs$idx2,
  chr1 = repeats$seqID[within_pairs$idx1],
  chr2 = repeats$seqID[within_pairs$idx2]
)

between_pairs_df <- data.frame(
  idx1 = between_pairs$idx1,
  idx2 = between_pairs$idx2,
  chr1 = repeats$seqID[between_pairs$idx1],
  chr2 = repeats$seqID[between_pairs$idx2]
)

within_pairs_df$within_similarities <- 0
within_pairs_df$within_distances <- 0

for(j in 1 : sample_within_chromosomes) {
  cat(j, "")
  within_pairs_df$within_similarities[j] <- adist(repeats$sequence[within_pairs_df$idx1[j]], 
                                                      repeats$sequence[within_pairs_df$idx2[j]])
  within_pairs_df$within_distances[j] <- abs(repeats$start[within_pairs_df$idx1[j]] - 
                                              repeats$start[within_pairs_df$idx2[j]])
}
cat("\n")

between_pairs_df$between_similarities <- 0

for(j in 1 : sample_between_chromosomes) {
  cat(j, "")
  between_pairs_df$between_similarities[j] <- adist(repeats$sequence[between_pairs_df$idx1[j]], 
                                                       repeats$sequence[between_pairs_df$idx2[j]])
}
cat("\n")
within_pairs_df$within_similarities <- within_pairs_df$within_similarities / mean(c(repeats$width[within_pairs_df$idx1], 
                                                                                    repeats$width[within_pairs_df$idx2]))
between_pairs_df$between_similarities <- between_pairs_df$between_similarities / mean(c(repeats$width[between_pairs_df$idx1], 
                                                                                    repeats$width[between_pairs_df$idx2]))

within_pairs_df$within_similarities <- 100 * (1 - within_pairs_df$within_similarities)
between_pairs_df$between_similarities <- 100 * (1 - between_pairs_df$between_similarities)

within_pairs_df$within_similarities[within_pairs_df$within_similarities < 0] = 0
between_pairs_df$between_similarities[between_pairs_df$between_similarities < 0] = 0



write.csv(within_pairs_df, file = paste0(assembly_name_short, "_within_pairs_df.csv"), row.names = FALSE)
write.csv(between_pairs_df, file = paste0(assembly_name_short, "_between_pairs_df.csv"), row.names = FALSE)


write.csv(within_pairs_df, file = paste0("/home/pwlodzimierz/ToL/upload_files/21_histograms_similarity_within_between_csvs/", assembly_name_short, "_within_pairs_df.csv"), row.names = FALSE)
write.csv(between_pairs_df, file = paste0("/home/pwlodzimierz/ToL/upload_files/21_histograms_similarity_within_between_csvs/", assembly_name_short, "_between_pairs_df.csv"), row.names = FALSE)








### PLOTTING

# Calculate distance percentiles for color scaling
dist_percentiles <- quantile(within_pairs_df$within_distances, probs = c(0.10, 0.90), na.rm = TRUE)
green_threshold <- dist_percentiles[1]  # 10th percentile
red_threshold <- dist_percentiles[2]    # 90th percentile


# Function to map distance to color (green to red)
distance_to_color <- function(distance) {
  norm_dist <- pmin(pmax((distance - green_threshold) / (red_threshold - green_threshold), 0), 1)
  red <- round(255 * norm_dist)
  green <- round(255 * (1 - norm_dist))
  return(rgb(red, green, 0, maxColorValue = 255))
}

# Calculate x-axis limits (same for both histograms)
x_min <- min(c(within_pairs_df$within_similarities, between_pairs_df$between_similarities), na.rm = TRUE)
x_max <- max(c(within_pairs_df$within_similarities, between_pairs_df$between_similarities), na.rm = TRUE)

# Calculate mean similarities
mean_within <- mean(within_pairs_df$within_similarities, na.rm = TRUE)
mean_between <- mean(between_pairs_df$between_similarities, na.rm = TRUE)

# Prepare data for within histogram
within_bins <- within_pairs_df %>%
  mutate(bin = cut(within_similarities, breaks = seq(x_min, x_max, length.out = 20))) %>%
  group_by(bin) %>%
  summarise(mean_distance = mean(within_distances, na.rm = TRUE),
            count = n()) %>%
  mutate(color = sapply(mean_distance, distance_to_color)) %>%
  mutate(bin_mid = as.numeric(gsub("\\((.*),.*\\]", "\\1", bin)) + 
           (as.numeric(gsub(".*,(.+)\\]", "\\1", bin)) - 
              as.numeric(gsub("\\((.*),.*\\]", "\\1", bin))) / 2)

# Prepare data for between histogram
between_bins <- between_pairs_df %>%
  mutate(bin = cut(between_similarities, breaks = seq(x_min, x_max, length.out = 20))) %>%
  group_by(bin) %>%
  summarise(count = n()) %>%
  mutate(color = "grey50") %>%
  mutate(bin_mid = as.numeric(gsub("\\((.*),.*\\]", "\\1", bin)) + 
           (as.numeric(gsub(".*,(.+)\\]", "\\1", bin)) - 
              as.numeric(gsub("\\((.*),.*\\]", "\\1", bin))) / 2)

# Create within histogram
p_within <- ggplot(within_bins, aes(x = bin_mid, y = count, fill = color)) +
  geom_bar(stat = "identity", width = (x_max - x_min) / 20 * 0.45) +
  geom_vline(xintercept = mean_within, color = "blue", linetype = "dashed", size = 1) +
  scale_fill_identity() +
  labs(title = "Within vs Between Chromosome Similarities",
       x = "Similarity (Edit Distance)",
       y = "Count (Within)") +
  theme_minimal() +
  theme(plot.margin = margin(b = 0)) +
  coord_cartesian(xlim = c(x_min, x_max))

# Create between histogram
p_between <- ggplot(between_bins, aes(x = bin_mid, y = count, fill = color)) +
  geom_bar(stat = "identity", width = (x_max - x_min) / 20 * 0.45) +
  geom_vline(xintercept = mean_between, color = "blue", linetype = "dashed", size = 1) +
  scale_fill_identity() +
  labs(x = "Similarity (Edit Distance)",
       y = "Count (Between)") +
  theme_minimal() +
  theme(plot.title = element_blank(),
        plot.margin = margin(t = 0)) +
  coord_cartesian(xlim = c(x_min, x_max))

# Combine plots vertically
pdf(paste0(assembly_name_short, "_hist_similarity_within_between.pdf"), width = 8, height = 6)
plot_grid(p_within, p_between, ncol = 1, align = "v", rel_heights = c(1, 1))
dev.off()

pdf(paste0("/home/pwlodzimierz/ToL/upload_files/21_histograms_similarity_within_between/", assembly_name_short, "_hist_similarity_within_between.pdf"), width = 8, height = 6)
plot_grid(p_within, p_between, ncol = 1, align = "v", rel_heights = c(1, 1))
dev.off()








