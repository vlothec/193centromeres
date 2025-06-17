#!/usr/bin/env Rscript
.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))

# f. HOR score in a 100 kbp distance from gaps. For each repeat, check which TEs are in 100 kbp distance, for each case save into a df the two
#    parameters: HOR score and physical distance
# --> 18_hor_score...R two-column df csv 
taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 193
print(i)

suppressMessages(library(seqinr))
library(ggplot2)
library(dplyr)

max_dist <- 100000



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

repeats <- read.csv(repeats_file); cat("repeats in\n")
edta <- read.csv(edta_file); cat("edta in\n")
arrays <- read.csv(arrays_file); cat("arrays in\n")
gaps <- read.csv(gaps_file); cat("gaps in\n")

satellites_metadata <- read.csv(file = "/home/pwlodzimierz/ToL/curated_satellites_repDec24_jan2025.csv")
chr_sizes_metadata <- read.csv(file = "/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full.Ian.csv")

assembly_name_long <- strsplit(data_directories[i], split = "/")[[1]][7]
assembly_name_short <- strsplit(assembly_name_long, split = ".fa")[[1]][1]

hor_repeats_file <- list.files(path = paste0("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs/self_HOR_out_loose_settings/", assembly_name_long), 
                               pattern = "HOR_scored_repeats_with_hors", full.names = TRUE)
if(length(repeats_file) == 0) {cat("issue with HOR repeats files: \n", repeats_file); stop()}

hor_repeats <- NULL
for(j in seq_along(hor_repeats_file)) {
  hor_repeats <- rbind(hor_repeats, read.csv(file = hor_repeats_file[j]))
}

### filter the data


satellites_metadata <- satellites_metadata[satellites_metadata$Genome == assembly_name_long, ]

satellite_TRASH_names <- NULL
for(j in seq_len(nrow(satellites_metadata))) {
  satellite_TRASH_names <- c(satellite_TRASH_names, strsplit(satellites_metadata$TRASH_name_dec2024runs[j], split = ";")[[1]])
}
chr_sizes_metadata <- chr_sizes_metadata[chr_sizes_metadata$assembly.name == assembly_name_long, ]

repeats <- repeats[repeats$new_class %in% satellite_TRASH_names, ]

###
gaps_distances <- NULL
HOR_scores <- NULL

for(j in seq_len(nrow(chr_sizes_metadata))) {
  cat(j, "/", nrow(chr_sizes_metadata), "\n")
  
  chr_repeats <- hor_repeats[hor_repeats$seqID == chr_sizes_metadata$chromosome.name[j], ]
  repeat_mids <- round(chr_repeats$start + (chr_repeats$end - chr_repeats$start) / 2)
  
  
  chr_gaps <- gaps[gaps$chromosome == chr_sizes_metadata$chromosome.name[j], ]
  chr_gaps$mid <- round(chr_gaps$start + (chr_gaps$end - chr_gaps$start) / 2)
  
  for(k in seq_along(chr_gaps)) {
    local_dists <- abs(chr_gaps$start[k] - repeat_mids)
    local_dists2 <- abs(chr_gaps$end[k] - repeat_mids)
    local_dists <- unlist(lapply(X = seq_along(local_dists), function(X) min(local_dists[X], local_dists2[X])))
    gaps_distances <- c(gaps_distances, local_dists[which(local_dists <= max_dist)])
    HOR_scores <- c(HOR_scores, chr_repeats$HOR_score[which(local_dists <= max_dist)])
  }
  
}

write.csv(x = data.frame(gaps_distances = gaps_distances, HOR_scores = HOR_scores), file = paste0(assembly_name_short, "_hor_score_vs_gaps_dist.csv"), row.names = FALSE)



data <- data.frame(gaps_distances = gaps_distances, HOR_scores = HOR_scores)

# Define the bin breaks (0 to 100,000 in steps of 10,000)
breaks <- c(0,100,500,1000,2500,5000,seq(10000, 100000, by = 10000))

all_bins <- paste0(breaks[-length(breaks)] + 1, "-", breaks[-1])

data <- data %>%
  mutate(distance_bin = cut(gaps_distances, 
                            breaks = breaks, 
                            labels = paste0(breaks[-length(breaks)] + 1, "-", breaks[-1]),
                            include.lowest = TRUE))

# Convert distance_bin to factor with all levels
data$distance_bin <- factor(data$distance_bin, levels = all_bins)

# Recreate the box plot
pdf(paste0(assembly_name_short, "_boxplot_by_distance_bins.pdf"))
ggplot(data, aes(x = distance_bin, y = HOR_scores)) +
  geom_boxplot(fill = "lightblue", color = "black", outlier.color = "red") +
  labs(title = "Box Plot of HOR Scores by Gaps Distance Bins",
       x = "Gaps Distance Bins", y = "HOR Scores") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


pdf(paste0("/home/pwlodzimierz/ToL/upload_files/18_hor_vs_gap_dist/", assembly_name_short, "_boxplot_by_distance_bins.pdf"))
ggplot(data, aes(x = distance_bin, y = HOR_scores)) +
  geom_boxplot(fill = "lightblue", color = "black", outlier.color = "red") +
  labs(title = "Box Plot of HOR Scores by Gaps Distance Bins",
       x = "Gaps Distance Bins", y = "HOR Scores") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


















