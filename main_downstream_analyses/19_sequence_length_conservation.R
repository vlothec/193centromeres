#!/usr/bin/env Rscript
.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))

# e. sequence length conservation: sample 10 000 repeat pairs, align them, check how many insertions, deletions and substitutions there are from
#    the perspective of one of the repeats. Make bar plots for each insertion:deletion combination with y axis being substitutions number: 0-0, 0-1, 
#    1-0, 1-1, ect. The idea is to show that it's more common to have the 1-1, 2-2, 3-3 etc. arrangement, that preserve the repeat length, let's call
#    it "compensating mutations"
# --> 19_sequence_...R occurances of each mutation csv 
taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 193
print(i)

suppressMessages(library(seqinr))
library(Biostrings)
library(stringdist)
library(parallel)


samples = 1000


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


satellite_TRASH_names <- NULL
for(j in seq_len(nrow(satellites_metadata))) {
  satellite_TRASH_names <- c(satellite_TRASH_names, strsplit(satellites_metadata$TRASH_name_dec2024runs[j], split = ";")[[1]])
}

chr_sizes_metadata <- chr_sizes_metadata[chr_sizes_metadata$assembly.name == assembly_name_long, ]

repeats <- repeats[repeats$new_class %in% satellite_TRASH_names, ]

###

# Sample random pairs
pairs <- data.frame(seq1 = sample(repeats$sequence, samples, replace = TRUE),
                    seq2 = sample(repeats$sequence, samples, replace = TRUE))


align_and_count <- function(seq1, seq2) {
  s1 <- as.character(seq1)
  s2 <- as.character(seq2)
  alignment <- pairwiseAlignment(s1, s2, type = "global")
  aln1 <- as.character(alignedPattern(alignment))
  aln2 <- as.character(alignedSubject(alignment))
  ins <- sum(charToRaw(aln1) == charToRaw("-"))
  del <- sum(charToRaw(aln2) == charToRaw("-"))
  sub <- sum(charToRaw(aln1) != charToRaw(aln2)) - ins - del
  c(insertions = ins, deletions = del, substitutions = sub)
}

# Apply to all pairs
results <- t(apply(pairs, 1, function(x) align_and_count(x[1], x[2])))

# Combine results
final_results <- cbind(pairs, results)
final_results$indels <- final_results$insertions + final_results$deletions
final_results$ins_minus_dels <- final_results$insertions - final_results$deletions
final_results$size_dif <- abs(nchar(final_results$seq1) - nchar(final_results$seq2))

write.csv(x = final_results[,c(3,4,5,6,7)], file = paste0(assembly_name_short, "_sequence_length_data.csv"), row.names = F)

pdf(paste0(assembly_name_short, "_sequence_indels_substitutions_histogram.pdf"), onefile = T)
hist(final_results$indels[final_results$indels < 25], breaks = 0:25, 
     main = paste0(assembly_name_short, " indels histogram"), right = F)
hist(final_results$indels[final_results$indels < 25] - final_results$size_dif[final_results$indels < 25], breaks = -5:25, 
     main = paste0(assembly_name_short, " indels minus size dif histogram"), right = F)
hist(final_results$substitutions[final_results$substitutions < 50], breaks = 0:50, 
     main = paste0(assembly_name_short, " substitutions histogram"), right = F)
dev.off()


pdf(paste0("/home/pwlodzimierz/ToL/upload_files/19_seq_length_conservation/", assembly_name_short, "_sequence_indels_substitutions_histogram.pdf"), onefile = T)
hist(final_results$indels[final_results$indels < 25], breaks = 0:25, 
     main = paste0(assembly_name_short, " indels histogram"), right = F)
hist(final_results$ins_minus_dels[final_results$ins_minus_dels < 15 & final_results$ins_minus_dels > -15], breaks = -15:15, 
     main = paste0(assembly_name_short, " ins minus dels histogram"), right = F)
hist(final_results$indels[final_results$indels < 25] - final_results$size_dif[final_results$indels < 25], breaks = -5:25, 
     main = paste0(assembly_name_short, " indels minus size dif histogram"), right = F)
hist(final_results$substitutions[final_results$substitutions < 50], breaks = 0:50, 
     main = paste0(assembly_name_short, " substitutions histogram"), right = F)
dev.off()






