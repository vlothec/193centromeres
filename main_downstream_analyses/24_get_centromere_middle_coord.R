#!/usr/bin/env Rscript
.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))

taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 193
print(i)


suppressMessages(library(Biostrings))
suppressMessages(library(GenomicRanges))
suppressMessages(library(Rsamtools))
suppressMessages(library(ggplot2))


### load data
data_directories <- list.dirs(path = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs", recursive = FALSE, full.names = TRUE)

assembly_name_long <- strsplit(data_directories[i], split = "/")[[1]][7]
assembly_name_short <- strsplit(assembly_name_long, split = ".fa")[[1]][1]

assembly_files <- list.files(path = "/home/pwlodzimierz/ToL/Assemblies/fastas_2021_Michael", recursive = FALSE, full.names = TRUE)
assembly_files <- assembly_files[-grep("fa.fai", assembly_files)]

assembly_file <- assembly_files[grep(assembly_name_short, assembly_files)]

setwd(data_directories[i])

repeats_file <- list.files(path = ".", pattern = "repeats_filtered.csv", full.names = TRUE)
if(length(repeats_file) != 1) {cat("issue with repeats files: \n", repeats_file); stop()}

repeats <- read.csv(repeats_file); cat("repeats in\n")

satellites_metadata <- read.csv(file = "/home/pwlodzimierz/ToL/curated_satellites_repDec24_jan2025.csv")
chr_sizes_metadata <- read.csv(file = "/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full.Ian.csv")

satellites_metadata <- satellites_metadata[satellites_metadata$Genome == assembly_name_long, ]


chr_sizes_metadata <- chr_sizes_metadata[chr_sizes_metadata$assembly.name == assembly_name_long, ]
chr_sizes_metadata <- chr_sizes_metadata[chr_sizes_metadata$is.chr == 1, ]

TRASH_names <- NULL
if(nrow(satellites_metadata) != 0) {
  for(j in 1 : nrow(satellites_metadata)) {
    TRASH_names <- c(TRASH_names, strsplit(satellites_metadata$TRASH_name_dec2024runs[j], split = ";")[[1]])
  }
}
if(length(TRASH_names) == 0) stop(paste0("there are no centromeric repeats on ", assembly_name_short))

repeats_full <- repeats 

table_classes <- sort(table(repeats$new_class), decreasing = T)
table_classes <- table_classes[table_classes >= 50]

repeats_new <- NULL
for(j in seq_along(table_classes)) {
  repeats_new <- rbind(repeats_new, repeats[repeats$new_class == names(table_classes)[j],][sample(1 : nrow(repeats[repeats$new_class == names(table_classes)[j], ]), 50), ])
}

repeats <- repeats_new

if(nrow(repeats) == 0) stop("WTF")

chr_sizes_metadata$repeats_mean <- NA

for(j in 1 : length(chr_sizes_metadata$chromosome.name)) {
  chr_repeats <- repeats[repeats$seqID == chr_sizes_metadata$chromosome.name[j],]
  
  chr_sizes_metadata$repeats_mean[j] = mean(chr_repeats$start)
  
}

write.csv(x = chr_sizes_metadata, file = paste0(assembly_name_short, "_chr_sizes_metadata_with_cen_middle.csv"), row.names = FALSE)
write.csv(x = chr_sizes_metadata, file = paste0("/home/pwlodzimierz/ToL/upload_files/24_array_middle_coordinates/", assembly_name_short, "_chr_sizes_metadata_with_cen_middle.csv"), row.names = FALSE)
































