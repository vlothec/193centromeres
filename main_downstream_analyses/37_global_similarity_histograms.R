
taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 193
print(i)


.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))

library(plot.matrix)

repeat_pairs_to_sample <- 1000
repeat_pairs_to_sample_long_repeats <- 100

setwd("/home/pwlodzimierz/ToL/git_ToL")
source("./aux_fun.R")

data_directories <- list.dirs(path = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs", recursive = FALSE, full.names = TRUE)
data_directories <- data_directories[!grepl(pattern = "templated_", data_directories)]
data_directories <- data_directories[grepl(pattern = ".fa", data_directories)]

assembly_files <- list.files(path = "/home/pwlodzimierz/ToL/Assemblies/fastas_2021_Michael", recursive = FALSE, full.names = TRUE)
assembly_files <- assembly_files[!grepl(".fai", assembly_files)]


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

genome_metadata <- read.csv("/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full.Ian.csv")
genome_metadata <- genome_metadata[genome_metadata$assembly.name == fasta_name, ]
chromosomes <- genome_metadata$chromosome.name[genome_metadata$is.chr == 1]


if(nrow(centromeric_stallite_for_species) == 0) stop("No centromeric satellites")

setwd(data_directories[i])
print(getwd())


repeat_file = list.files(pattern = "_repeats_filtered.csv", full.names = TRUE)
if(length(repeat_file) != 1) {warning(paste0(i, "No repeats!")); setwd(".."); quit(save = "no", status = 1)}


print("load annotations")
repeats = read.csv(file = repeat_file)

if(fasta_name %in% c("mRhiSin1.1.fa", "iyVesVulg1.1.fa")) { 
  repeats$new_class <- repeats$class
  repeat_pairs_to_sample <- repeat_pairs_to_sample_long_repeats
}


for(j in 1 : nrow(centromeric_stallite_for_species)) {

  cen_repeats <- repeats[repeats$new_class == strsplit(centromeric_stallite_for_species$TRASH_name_dec2024runs[j], split = ";")[[1]][1],]
  
  # set up pairs to compare within and between chromosomes
  
  same_seqid_pairs <- data.frame(index1 = integer(), index2 = integer())
  diff_seqid_pairs <- data.frame(index1 = integer(), index2 = integer())
  
  seq_ids <- unique(cen_repeats$seqID)
  
  chromosomes <- genome_metadata$chromosome.name[genome_metadata$is.chr == 1]
  
  seq_ids <- seq_ids[seq_ids %in% chromosomes]
  
  while (nrow(same_seqid_pairs) < repeat_pairs_to_sample) {
    selected_seqid <- sample(seq_ids, 1)
    seqid_indices <- which(cen_repeats$seqID == selected_seqid)
    
    if (length(seqid_indices) >= 2) {
      sampled_indices <- sample(seqid_indices, 2, replace = FALSE)
      same_seqid_pairs <- rbind(same_seqid_pairs, 
                                data.frame(index1 = sampled_indices[1], 
                                           index2 = sampled_indices[2]))
    }
  }
  
  while (nrow(diff_seqid_pairs) < repeat_pairs_to_sample) {
    print(nrow(diff_seqid_pairs))
    selected_seqids <- sample(seq_ids, 2, replace = FALSE)
    
    index1 <- sample(which(cen_repeats$seqID == selected_seqids[1]), 1)
    index2 <- sample(which(cen_repeats$seqID == selected_seqids[2]), 1)
    
    diff_seqid_pairs <- rbind(diff_seqid_pairs,
                              data.frame(index1 = index1,
                                         index2 = index2))
  }
  
  # calculate similarity
  same_seqid_pairs$width1 <- cen_repeats$width[same_seqid_pairs$index1]
  same_seqid_pairs$width2 <- cen_repeats$width[same_seqid_pairs$index2]
  
  same_seqid_pairs$similarity <- unlist(lapply(1 : nrow(same_seqid_pairs), function(X) {
    100-100*adist(cen_repeats$sequence[same_seqid_pairs$index1[X]], cen_repeats$sequence[same_seqid_pairs$index2[X]])/mean(c(same_seqid_pairs$width1[X], same_seqid_pairs$width2[X]))
  }))
  
  diff_seqid_pairs$width1 <- cen_repeats$width[diff_seqid_pairs$index1]
  diff_seqid_pairs$width2 <- cen_repeats$width[diff_seqid_pairs$index2]
  
  diff_seqid_pairs$similarity <- unlist(lapply(1 : nrow(same_seqid_pairs), function(X) {
    100-100*adist(cen_repeats$sequence[diff_seqid_pairs$index1[X]], cen_repeats$sequence[diff_seqid_pairs$index2[X]])/mean(c(diff_seqid_pairs$width1[X], diff_seqid_pairs$width2[X]))
  }))
  
  same_seqid_pairs$method <- "same_chr"
  diff_seqid_pairs$method <- "diff_chr"
  
  merged_df <- rbind(same_seqid_pairs, diff_seqid_pairs)
  
  merged_df$similarity[merged_df$similarity < 0] = 0
  merged_df$similarity[merged_df$similarity > 100] = 100
  
  # get similarity between specific chromosomes
  
  merged_df$chr1 <- cen_repeats$seqID[merged_df$index1]
  merged_df$chr2 <- cen_repeats$seqID[merged_df$index2]
  
  write.csv(x = merged_df,
            file = paste0("/home/pwlodzimierz/ToL/upload_files/37_similarity_values_within_between_chr/similarities_within_between_chrs_", fasta_name, "_", centromeric_stallite_for_species$Satellite_name[j], ".csv"), row.names = F)
  
  
  
  chrs <- sort(unique(c(merged_df$chr1, merged_df$chr2)))
  
  between_chr_matrix <- matrix(nrow = length(chrs), ncol = length(chrs))
  
  for(idY in 1 : nrow(between_chr_matrix)) {
    for(idX in 1 : ncol(between_chr_matrix)) {
      between_chr_matrix[idX,idY] = mean(merged_df$similarity[(merged_df$chr1 == chrs[idY] & merged_df$chr2 == chrs[idX]) | (merged_df$chr1 == chrs[idX] & merged_df$chr2 == chrs[idY])])
    }
  }
  
  rownames(between_chr_matrix) <- chrs
  colnames(between_chr_matrix) <- chrs
  
  
  xmin = 0
  xmax = 100
  
  pdf(file = paste0("/home/pwlodzimierz/ToL/upload_files/37_similarity_histograms_within_between_chr/", fasta_name, "_", centromeric_stallite_for_species$Satellite_name[j], "_histogrmas_similarity_within_between_chrs.pdf"), width = 12, height = 18)
  par(mfrow = c(3,1))
  hist(merged_df$similarity[merged_df$method == "same_chr"], 
       breaks = seq(xmin,xmax, by = 1), xlim = c(xmin, xmax), border = "#eeee0080", col = "#eeee0080",
       xlab = "per chromosome mean pairwise similarity of repeats", main = "Histogram of centromeric pairwise similarities")
  hist(merged_df$similarity[merged_df$method == "diff_chr"], 
       breaks = seq(xmin,xmax, by = 1), xlim = c(xmin, xmax), add = TRUE, border = "#00ee0080", col = "#00ee0080")
  # legend(x = 1, y = 50, legend = c("chromosomes", "genomes"), fill = c("#eeee0080", "#00ee0080"))
  
  chromosomes <- merged_df$similarity[merged_df$method == "same_chr"]
  genomes     <- merged_df$similarity[merged_df$method == "diff_chr"]
  
  # Plot the boxplot
  boxplot(chromosomes, genomes,
          names = c("chromosomes", "genomes"),
          main = "pairwise similarity of repeats within:",
          col = c("#eeee0080", "#00ee0080"),
          ylim = c(0, 140))
  
  # Perform pairwise t-tests
  p1 <- t.test(chromosomes, genomes)$p.value
  
  # Add asterisks based on significance levels
  # Positioning settings
  y_max <- max(c(chromosomes, genomes), na.rm = TRUE)
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
  
  
  color_func <- colorRampPalette(c("red", "yellow", "green"))
  
  plot(between_chr_matrix, 
       main = "Chromosome Similarity Matrix",
       col = color_func(18),
       # breaks = c(0,seq(80, 100, length.out = 101)),
       breaks = c(0,10,20,30,40,50,60,70,80,82,84,86,88,90,92,94,96,98,100),
       fmt.cell = "%.1f",
       text.cell = list(cex = 1),
       axis.col = list(cex.axis = 0.8),
       axis.row = list(cex.axis = 0.8),
       xlab = "Chromosome",
       ylab = "Chromosome")
  
  
  dev.off()
  
}
  
  

















