min_gap_distance <- 1


###########


.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))
suppressMessages(library(msa))
suppressMessages(library(seqinr))
suppressMessages(library(GenomicRanges))
suppressMessages(library(ggplot2))

setwd("/home/pwlodzimierz/ToL/git_ToL")
source("./aux_fun.R")
ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}

data_directories <- list.dirs(path = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs", recursive = FALSE, full.names = TRUE)
data_directories <- data_directories[!grepl(pattern = "templated_", data_directories)]
data_directories <- data_directories[grepl(pattern = ".fa", data_directories)]
assembly_files <- list.files(path = "/home/pwlodzimierz/ToL/Assemblies/fastas_2021_Michael", recursive = FALSE, full.names = TRUE)
assembly_files <- assembly_files[!grepl(".fai", assembly_files)]

satellite_metadata <- read.csv("/home/pwlodzimierz/ToL/curated_satellites_metadata_on_chromosomes_only_may.csv")

taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 192
print(i)

print(paste0(i, " / ", length(data_directories)))
### Load data ================================================================
setwd(data_directories[i])
print(getwd())

assembly_name = strsplit(strsplit(data_directories[i], split = ".fa")[[1]][1], split = "v2_out_for_HORs/")[[1]][2]
fasta_name = strsplit(data_directories[i], split = "v2_out_for_HORs/")[[1]][2]
assembly_file = grep(assembly_name, assembly_files)
print(assembly_file)
print(assembly_files[assembly_file])

if(fasta_name %in% satellite_metadata$Genome) {
  cat("Analysing the ", fasta_name, "genome\n")
} else {
  cat("The genome is not satellite-based, exiting\n")
  quit(save = "no", status = 1)
}



repeat_file = list.files(pattern = "_repeats_filtered_array_assigned.csv", full.names = TRUE)
if(length(repeat_file) != 1) {warning(paste0(i, "No repeats!")); setwd(".."); quit(save = "no", status = 1)}


print("load annotations")
repeats = read.csv(file = repeat_file)
print(nrow(repeats))
print(names(repeats))
print("annotations loaded")

unique(repeats$is_centromeric)

repeats <- repeats[repeats$is_centromeric == "TRUE",]
repeats <- repeats[order(repeats$start, decreasing = FALSE), ]
repeats <- repeats[order(repeats$seqID, decreasing = FALSE), ]
print(nrow(repeats))

if(nrow(repeats) == 0) {
  stop("no centromeric repeats")
}


array_ids <- unique(repeats$arrayID)
array_ids <- array_ids[array_ids != ""]

gaps_total <- NULL

for(j in seq_along(array_ids)) {
  cat(array_ids[j], "\n")
  array_repeats <- repeats[repeats$arrayID == array_ids[j], ]
  
  if(nrow(array_repeats) < 3) next
  
  array_repeats$dist_to_next = 0
  array_repeats$dist_to_next[1:(nrow(array_repeats) - 1)] = 
    array_repeats$start[2:nrow(array_repeats)] - 
    array_repeats$end[1:(nrow(array_repeats) - 1)] - 1
  array_repeats$dist_to_next[nrow(array_repeats)] = 0
  
  gaps_starts <- array_repeats$end[array_repeats$dist_to_next >= min_gap_distance] + 1
  gaps_ends <- array_repeats$end[array_repeats$dist_to_next >= min_gap_distance] + array_repeats$dist_to_next[array_repeats$dist_to_next >= min_gap_distance] 
  gaps_repeat_array_id <- rep(array_ids[j], length(gaps_starts))
  gaps_species <- rep(fasta_name, length(gaps_starts))
  gaps_chromosomes <- rep(array_repeats$seqID[1], length(gaps_starts))
  
  if(length(gaps_starts) > 0) {
    
    gaps <- data.frame(start = gaps_starts, end = gaps_ends, repeat_array_id = gaps_repeat_array_id, 
                       species = gaps_species, chromosome = gaps_chromosomes)
    
    gaps_total <- rbind(gaps_total, gaps)
  }
  
}



write.csv(x = gaps_total, file = paste0("./", fasta_name, "_centromeric_gaps.csv"), row.names = FALSE )



write.csv(x = gaps_total, file = paste0("/home/pwlodzimierz/ToL/upload_files/centromeric_gaps/", fasta_name, "_centromeric_gaps.csv"), row.names = FALSE )












