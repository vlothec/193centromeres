
taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 193
print(i)

redo = F

.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))
suppressMessages(library(msa))
suppressMessages(library(seqinr)) 
suppressMessages(library(GenomicRanges))
suppressMessages(library(ggplot2))

ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}


setwd("/home/pwlodzimierz/ToL/git_ToL")
source("./aux_fun.R")

data_directories <- list.dirs(path = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs", recursive = FALSE, full.names = TRUE)
data_directories <- data_directories[!grepl(pattern = "templated_", data_directories)]
data_directories <- data_directories[grepl(pattern = ".fa", data_directories)]

assembly_files <- list.files(path = "/home/pwlodzimierz/ToL/Assemblies/fastas_2021_Michael", recursive = FALSE, full.names = TRUE)
assembly_files <- assembly_files[!grepl(".fai", assembly_files)]

if(i > length(data_directories)) stop()
if(grepl("Rosa_agrestis_DTOL_chrs_5n.fasta", data_directories[i])) stop()
if(grepl("Rosa_canina_DTOL_chrs_5n.fasta", data_directories[i])) stop()
if(grepl("drGeuRiva1.hap1.1.fa", data_directories[i])) stop()
if(grepl("drGeuRiva1.hap2.1.fa", data_directories[i])) stop()
if(grepl("drRosRugo1.1.fa", data_directories[i])) stop()
if(grepl("drRosSpin1.hap1.1.fa", data_directories[i])) stop()
if(grepl("GCA_000001405.29_GRCh38.p14_genomic.fa", data_directories[i])) stop()

print(paste0(i, " / ", length(data_directories)))
### Load data ================================================================
assembly_name = strsplit(strsplit(data_directories[i], split = ".fa")[[1]][1], split = "v2_out_for_HORs/")[[1]][2]
fasta_name = strsplit(data_directories[i], split = "v2_out_for_HORs/")[[1]][2]
assembly_file = grep(assembly_name, assembly_files)
print(assembly_file)
print(assembly_files[assembly_file])

genomes_organisation_data <- read.csv("/home/pwlodzimierz/ToL/genomes_organisation_type.csv")
genomes_organisation_data <- genomes_organisation_data[genomes_organisation_data$fasta == fasta_name,]

satellite_metadata <- read.csv("/home/pwlodzimierz/ToL/curated_satellites_metadata_on_chromosomes_only_may.csv")
centromeric_stallite_for_species <- satellite_metadata[satellite_metadata$Genome == fasta_name,]

setwd(data_directories[i])
print(getwd())


if(!(genomes_organisation_data$Centromere.architecture %in% c("Satellite", "Holocentric satellite"))) {
  print("Not a satellite species")
  quit(save = "no")
}


for(satID in 1 : nrow(centromeric_stallite_for_species)) {
  
  satellite_name <- centromeric_stallite_for_species$Satellite_name_current[satID]
  trash_satellite_names <- strsplit(centromeric_stallite_for_species$TRASH_name_dec2024runs[satID], split = ";")[[1]]
  
  if(!redo) {
    if(file.exists(paste0("/home/pwlodzimierz/ToL/upload_files/13_hor_grand_tables/", fasta_name, "_", genomes_organisation_data$Centromere.architecture[1], "_", centromeric_stallite_for_species$Satellite_name_current[satID], "_HORs_grand_table.csv"))) {
      print("Already done, not redoing")
      next
    }
  }
  
  setwd("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs/self_HOR_out_loose_settings") 
  
  hor_dirs <- list.dirs(path = ".", recursive = FALSE)
  hor_dirs_full <- list.dirs(path = ".", recursive = FALSE, full.names = T)
  
  which.dir <- grep(fasta_name, hor_dirs)
  setwd(hor_dirs_full[which.dir])
  
  
  
  hor_files <- list.files(path = ".")
  
  all_hors_files <- hor_files[grep("HORs_", hor_files)]
  all_hors_files <- all_hors_files[grep("csv", all_hors_files)]
  
  repeats_with_hors_files <- hor_files[grep("HOR_scored_", hor_files)]
  
  all_hors_files_single_class <- NULL 
  
  for(k in seq_along(trash_satellite_names)) {
    all_hors_files_single_class <- c(all_hors_files_single_class, all_hors_files[grepl(trash_satellite_names[k], all_hors_files)])
  }
  
  
  
  repeats_with_hors_files_single_class <- NULL 
  
  for(k in seq_along(trash_satellite_names)) {
    repeats_with_hors_files_single_class <- c(repeats_with_hors_files_single_class, repeats_with_hors_files[grepl(trash_satellite_names[k], repeats_with_hors_files)])
  }
  
  
  
  if(length(repeats_with_hors_files_single_class) == 0) {
    cat("No repeats with HOR files found\n")
    next
  }
  if(length(all_hors_files) == 0) {
    cat("No HOR files found\n")
    next
  }
  
  
  repeats_all <- NULL
  hors_no = 0
  hors_mean_block_size = 0
  hors_mean_block_distance = 0
  
  for(j in seq_along(repeats_with_hors_files_single_class)) {
    print(paste0(j, "/", length(repeats_with_hors_files_single_class)))
    repeat_class <- paste(strsplit(repeats_with_hors_files_single_class[j], split = "_")[[1]][6:7], collapse = "_")
    
    chromosome_name <- strsplit(repeats_with_hors_files_single_class[j], split = paste0(repeat_class, "_"))[[1]][2]
    chromosome_name <- strsplit(chromosome_name, split = ".csv")[[1]][1]
    
    
    repeats_chr <- read.csv(file = paste0(repeats_with_hors_files_single_class[j]))
    if(!("HOR_score" %in% names(repeats_chr))) repeats_chr$HOR_score = 0
    repeats_all <- rbind(repeats_all, repeats_chr)
    
    
    if(file.exists(paste0("HORs_", repeat_class, "_", chromosome_name, ".csv"))) {
      hors <- read.csv(file = paste0("HORs_", repeat_class, "_", chromosome_name, ".csv"))
      
      hors_no <- hors_no + nrow(hors)
      hors_mean_block_size <- hors_mean_block_size + sum(hors$block.size.in.units)
      hors_mean_block_distance <- hors_mean_block_distance + sum(abs(hors$start.A.bp - hors$start.B.bp))
      
    }
    
  }
  
  hors_mean_block_size <- hors_mean_block_size / hors_no
  hors_mean_block_distance <- hors_mean_block_distance / hors_no
  #
  repeats_no <- nrow(repeats_all)
  repeats_total_bp <- sum(repeats_all$width)
  repeats_mean_HOR_score <- mean(repeats_all$HOR_score)
  
  
  write.csv(x = data.frame(fasta_name = fasta_name,
                           group = genomes_organisation_data$Group,
                           repeats_no = repeats_no,
                           repeats_total_bp = repeats_total_bp,
                           repeats_mean_HOR_score = repeats_mean_HOR_score,
                           hors_no = hors_no,
                           hors_mean_block_size = hors_mean_block_size,
                           hors_mean_block_distance = hors_mean_block_distance), 
            file = paste0("/home/pwlodzimierz/ToL/upload_files/13_hor_grand_tables/", fasta_name, "_", genomes_organisation_data$Centromere.architecture[1], "_", centromeric_stallite_for_species$Satellite_name_current[satID], "_HORs_grand_table.csv"), 
            row.names = F)
  
  setwd(data_directories[i])
  
  write.csv(x = data.frame(fasta_name = fasta_name,
                           group = genomes_organisation_data$Group,
                           repeats_no = repeats_no,
                           repeats_total_bp = repeats_total_bp,
                           repeats_mean_HOR_score = repeats_mean_HOR_score,
                           hors_no = hors_no,
                           hors_mean_block_size = hors_mean_block_size,
                           hors_mean_block_distance = hors_mean_block_distance), 
            file = paste0(fasta_name, "_", genomes_organisation_data$Centromere.architecture[1], "_", centromeric_stallite_for_species$Satellite_name_current[satID], "_HORs_grand_table.csv"), 
            row.names = F)
  
  
  
}









