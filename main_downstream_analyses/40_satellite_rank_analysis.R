

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


satellite_metadata <- read.csv("/home/pwlodzimierz/ToL/curated_satellites_metadata_on_chromosomes_only_may.csv")

all.proper.sats <- NULL

for(i in 1 : nrow(satellite_metadata)) { # 1 to 91

  centromeric_satellite_for_species <- satellite_metadata[i,]
  
  fasta_name <- centromeric_satellite_for_species$Genome
  
  setwd(data_directories[grep(fasta_name, data_directories)])
  
  print(paste0(i, " / ", nrow(satellite_metadata)))
  
  ### Load data ================================================================
  
  repeat_file = list.files(pattern = "_repeats_filtered.csv", full.names = TRUE)
  if(length(repeat_file) != 1) {warning(paste0(i, "No repeats!")); setwd(".."); quit(save = "no", status = 1)}
  
  repeats = read.csv(file = repeat_file)
  
  repeats$family_name <- satellite_metadata$Satellite_name_current[i]
  
  repeats <- repeats[repeats$new_class %in% strsplit(centromeric_satellite_for_species$TRASH_name_dec2024runs, split = ";")[[1]], ]
  
  repeats <- repeats[,c("sequence", "family_name")]
  
  all.proper.sats <- rbind(all.proper.sats, repeats)
  
}

write.csv(x = all.proper.sats, file = "/home/pwlodzimierz/ToL/upload_files/40_all_sat_seqs/all.sat.seqs.csv", row.names = FALSE)
# all.proper.sats = read.csv( "/home/pwlodzimierz/ToL/upload_files/40_all_sat_seqs/all.sat.seqs.csv")


all.proper.sats$fam <- ""
for(i in 1 : nrow(satellite_metadata)) { # 1 to 91
  print(i)
  all.proper.sats$fam[all.proper.sats$family_name == satellite_metadata$Satellite_name_current[i]] <- satellite_metadata$group2[i]
}

table_invertebrate <- table(all.proper.sats$sequence[all.proper.sats$fam == "invertebrate"])
table_vert <- table(all.proper.sats$sequence[all.proper.sats$fam == "vertebrate"])
table_plant <- table(all.proper.sats$sequence[all.proper.sats$fam == "plant"])


to_rank_freq <- function(freq_table) {
  freq <- sort(as.numeric(freq_table), decreasing = TRUE)
  if (length(freq) == 0) {
    warning("Empty frequency table")
    return(data.frame(rank = numeric(0), freq = numeric(0), log.rank = numeric(0), log.freq = numeric(0)))
  }
  valid <- freq > 0
  ranks <- seq_along(freq)
  log_ranks <- log(ranks[valid])
  log_freq <- log(freq[valid])
  data.frame(
    rank = ranks,
    freq = freq,
    log.rank = c(log_ranks, rep(NA, sum(!valid))),
    log.freq = c(log_freq, rep(NA, sum(!valid)))
  )
}

# Prepare data for each group
invert_data <- to_rank_freq(table_invertebrate)
vert_data <- to_rank_freq(table_vert)
plant_data <- to_rank_freq(table_plant)

invert_fit <- lm(log.freq ~ log.rank, data = invert_data[!is.na(invert_data$log.freq), ])
vert_fit <- lm(log.freq ~ log.rank, data = vert_data[!is.na(vert_data$log.freq), ])
plant_fit <- lm(log.freq ~ log.rank, data = plant_data[!is.na(plant_data$log.freq), ])


# pdf(file = "/home/pwlodzimierz/ToL/upload_files/40_all_sat_seqs/satellite_sequence_rank_analysis.pdf", width = 12, height = 12)
png("/home/pwlodzimierz/ToL/upload_files/40_all_sat_seqs/satellite_sequence_rank_analysis_ln.png", width = 900, height = 900)


plot(invert_data$log.rank, invert_data$log.freq, type = "o", 
     xlab = "log(Rank)", ylab = "log(Frequency)", main = "Invertebrates", col = "#3f37c9", pch = 1)
abline(invert_fit, col = "#3f37c9", lty = 2)

points(vert_data$log.rank, vert_data$log.freq, type = "o", 
     col = "#f72585", pch = 2)
abline(vert_fit, col = "#f72585", lty = 2)

points(plant_data$log.rank, plant_data$log.freq, type = "o", 
     col = "#8ac926", pch = 3)
abline(plant_fit, col = "#8ac926", lty = 2)

# legend("topright", 
#        legend = c(sprintf("Invertebrates (R²=%s)", format(summary(invert_fit)$r.squared, digits = 4) %||% "NA"),
#                   sprintf("Vertebrates (R²=%s)", format(summary(vert_fit)$r.squared, digits = 4) %||% "NA"),
#                   sprintf("Plants (R²=%s)", format(summary(plant_fit)$r.squared, digits = 4) %||% "NA")),
#        col = c("#3f37c9", "#f72585", "#8ac926"), 
#        pch = c(1, 2, 3), lty = 2)

dev.off()
















repeat_sequence_table_vert <- NULL  # RANK, divide by no of vert
repeat_sequence_perc_table_vert <- NULL # % of satellite sequence at rank, divide by no of vert

repeat_sequence_table_invert <- NULL
repeat_sequence_perc_table_invert <- NULL

repeat_sequence_table_plant <- NULL
repeat_sequence_perc_table_plant <- NULL

all_chromosome_sizes <- read.csv(file = "/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full.Ian.csv")

for(i in 1 : nrow(satellite_metadata)) { # 1 to 91
  setwd(paste0("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs/self_HOR_out_loose_settings/", satellite_metadata$Genome[i]))
  
  
  centromeric_satellite_for_species <- satellite_metadata[i,]
  
  fasta_name <- centromeric_satellite_for_species$Genome
  
  chromosome_sizes <- all_chromosome_sizes[grepl(fasta_name, all_chromosome_sizes$assembly.name), ]
  
  setwd(data_directories[grep(fasta_name, data_directories)])
  
  print(paste0(i, " / ", nrow(satellite_metadata)))
  
  ### Load data ================================================================
  
  repeat_file = list.files(pattern = "_repeats_filtered.csv", full.names = TRUE)
  if(length(repeat_file) != 1) {warning(paste0(i, "No repeats!")); setwd(".."); quit(save = "no", status = 1)}
  
  repeats = read.csv(file = repeat_file)
  
  repeats$family_name <- satellite_metadata$Satellite_name_current[i]
  
  repeats <- repeats[repeats$new_class %in% strsplit(centromeric_satellite_for_species$TRASH_name_dec2024runs, split = ";")[[1]], ]
  
  repeats <- repeats[repeats$seqID %in% chromosome_sizes$chromosome.name[chromosome_sizes$is.chr == 1], ]
  
  
  
  mono_temp_table <- sort(table(repeats$sequence), decreasing = TRUE)
  names(mono_temp_table) <- 1 : length(mono_temp_table)
  
  mono_temp_table_perc <- 100 * mono_temp_table / sum(mono_temp_table)
  
  
  if(satellite_metadata$group2[i] == "vertebrate") {
    repeat_sequence_table_vert <- c(repeat_sequence_table_vert, mono_temp_table[setdiff(names(mono_temp_table), names(repeat_sequence_table_vert))])
    repeat_sequence_table_vert[names(mono_temp_table)] <- repeat_sequence_table_vert[names(mono_temp_table)] + mono_temp_table
    
    
    repeat_sequence_perc_table_vert <- c(repeat_sequence_perc_table_vert, mono_temp_table_perc[setdiff(names(mono_temp_table_perc), names(repeat_sequence_perc_table_vert))])
    repeat_sequence_perc_table_vert[names(mono_temp_table)] <- repeat_sequence_perc_table_vert[names(mono_temp_table_perc)] + mono_temp_table_perc
    
  }
  if(satellite_metadata$group2[i] == "invertebrate") {
    repeat_sequence_table_invert <- c(repeat_sequence_table_invert, mono_temp_table[setdiff(names(mono_temp_table), names(repeat_sequence_table_invert))])
    repeat_sequence_table_invert[names(mono_temp_table)] <- repeat_sequence_table_invert[names(mono_temp_table)] + mono_temp_table
    
    
    repeat_sequence_perc_table_invert <- c(repeat_sequence_perc_table_invert, mono_temp_table_perc[setdiff(names(mono_temp_table_perc), names(repeat_sequence_perc_table_invert))])
    repeat_sequence_perc_table_invert[names(mono_temp_table)] <- repeat_sequence_perc_table_invert[names(mono_temp_table_perc)] + mono_temp_table_perc
    
  }
  if(satellite_metadata$group2[i] == "plant") {
    repeat_sequence_table_plant <- c(repeat_sequence_table_plant, mono_temp_table[setdiff(names(mono_temp_table), names(repeat_sequence_table_plant))])
    repeat_sequence_table_plant[names(mono_temp_table)] <- repeat_sequence_table_plant[names(mono_temp_table)] + mono_temp_table
    
    
    repeat_sequence_perc_table_plant <- c(repeat_sequence_perc_table_plant, mono_temp_table_perc[setdiff(names(mono_temp_table_perc), names(repeat_sequence_perc_table_plant))])
    repeat_sequence_perc_table_plant[names(mono_temp_table)] <- repeat_sequence_perc_table_plant[names(mono_temp_table_perc)] + mono_temp_table_perc
    
  }
  
  
  
  
}

write.csv(x = repeat_sequence_table_vert, file = "repeat_sequence_table_vert.csv")
write.csv(x = repeat_sequence_perc_table_vert, file = "repeat_sequence_perc_table_vert.csv")
write.csv(x = repeat_sequence_table_invert, file = "repeat_sequence_table_invert.csv")
write.csv(x = repeat_sequence_perc_table_invert, file = "repeat_sequence_perc_table_invert.csv")
write.csv(x = repeat_sequence_table_plant, file = "repeat_sequence_table_plant.csv")
write.csv(x = repeat_sequence_perc_table_plant, file = "repeat_sequence_perc_table_plant.csv")

repeat_sequence_table_vert <- repeat_sequence_table_vert / nrow(satellite_metadata[satellite_metadata$group2 == "vertebrate",])
repeat_sequence_perc_table_vert <- repeat_sequence_perc_table_vert / nrow(satellite_metadata[satellite_metadata$group2 == "vertebrate",])

repeat_sequence_table_invert <- repeat_sequence_table_invert / nrow(satellite_metadata[satellite_metadata$group2 == "invertebrate",])
repeat_sequence_perc_table_invert <- repeat_sequence_perc_table_invert / nrow(satellite_metadata[satellite_metadata$group2 == "invertebrate",])

repeat_sequence_table_plant <- repeat_sequence_table_plant / nrow(satellite_metadata[satellite_metadata$group2 == "plant",])
repeat_sequence_perc_table_plant  <- repeat_sequence_perc_table_plant / nrow(satellite_metadata[satellite_metadata$group2 == "plant",])

repeat_sequence_table_vert <- repeat_sequence_table_vert[!is.na(repeat_sequence_table_vert)]
repeat_sequence_perc_table_vert <- repeat_sequence_perc_table_vert[!is.na(repeat_sequence_perc_table_vert)]

repeat_sequence_table_invert <- repeat_sequence_table_invert[!is.na(repeat_sequence_table_invert)]
repeat_sequence_perc_table_invert <- repeat_sequence_perc_table_invert[!is.na(repeat_sequence_perc_table_invert)]

repeat_sequence_table_plant <- repeat_sequence_table_plant[!is.na(repeat_sequence_table_plant)]
repeat_sequence_perc_table_vert <- repeat_sequence_perc_table_vert[!is.na(repeat_sequence_perc_table_vert)]



setwd("/home/pwlodzimierz/ToL/upload_files/40_all_sat_seqs/")

# - HOR monomer length rank analysis
freq <- sort(as.numeric(repeat_sequence_table_vert), decreasing = TRUE)
valid <- freq > 0
ranks <- seq_along(freq)
log_ranks <- log(ranks[valid])
log_freq <- log(freq[valid])
vert_data <- data.frame(
  rank = ranks,
  freq = freq,
  log.rank = c(log_ranks, rep(NA, sum(!valid))),
  log.freq = c(log_freq, rep(NA, sum(!valid)))
)
vert_fit <- lm(log.freq ~ log.rank, data = vert_data[!is.na(vert_data$log.freq), ])

freq <- sort(as.numeric(repeat_sequence_table_invert), decreasing = TRUE)
valid <- freq > 0
ranks <- seq_along(freq)
log_ranks <- log(ranks[valid])
log_freq <- log(freq[valid])
invert_data <- data.frame(
  rank = ranks,
  freq = freq,
  log.rank = c(log_ranks, rep(NA, sum(!valid))),
  log.freq = c(log_freq, rep(NA, sum(!valid)))
)
invert_fit <- lm(log.freq ~ log.rank, data = invert_data[!is.na(invert_data$log.freq), ])

freq <- sort(as.numeric(repeat_sequence_table_plant), decreasing = TRUE)
valid <- freq > 0
ranks <- seq_along(freq)
log_ranks <- log(ranks[valid])
log_freq <- log(freq[valid])
plant_data <- data.frame(
  rank = ranks,
  freq = freq,
  log.rank = c(log_ranks, rep(NA, sum(!valid))),
  log.freq = c(log_freq, rep(NA, sum(!valid)))
)
plant_fit <- lm(log.freq ~ log.rank, data = plant_data[!is.na(plant_data$log.freq), ])



invert_data <- invert_data[-sample(which(invert_data$log.freq == min(invert_data$log.freq)), round(length(which(invert_data$log.freq == min(invert_data$log.freq)))/20)), ]
vert_data <- vert_data[-sample(which(vert_data$log.freq == min(vert_data$log.freq)), round(length(which(vert_data$log.freq == min(vert_data$log.freq)))/20)), ]
plant_data <- plant_data[-sample(which(plant_data$log.freq == min(plant_data$log.freq)), round(length(which(plant_data$log.freq == min(plant_data$log.freq)))/20)), ]



# png("./satellite_sequence_rank_analysis_ln.png", width = 900, height = 900)
pdf("./satellite_sequence_rank_analysis_ln.pdf", width = 5, height = 5)

plot(invert_data$log.rank, invert_data$log.freq, type = "o", 
     xlab = "log(Rank)", ylab = "log(Frequency)", main = "satellite_sequence_rank_analysis_ln", col = "#3f37c9", pch = 1, cex = 1.25,
     xlim = c(0,14), ylim = c(-3,9))
abline(invert_fit, col = "#3f37c9", lty = 2)

points(vert_data$log.rank, vert_data$log.freq, type = "o", cex = 1.25,
       col = "#f72585", pch = 1)
abline(vert_fit, col = "#f72585", lty = 2)

points(plant_data$log.rank, plant_data$log.freq, type = "o", cex = 1.25,
       col = "#8ac926", pch = 1)
abline(plant_fit, col = "#8ac926", lty = 2)

legend("topright",
       legend = c("Plant", "Invertebrate", "Chordate"),
       text.col = c("#8ac926", "#3f37c9", "#f72585"), bty = "n")

dev.off()



# - HOR monomer length histogram (use % instead of frequency) 

pdf("./satellite_sequence_histogram.pdf", width = 5, height = 5)

plot(x = as.numeric(names(repeat_sequence_perc_table_vert)), y = repeat_sequence_perc_table_vert, type = "o", 
     xlab = "HOR monomer length", ylab = "% of HORs", main = "satellite_sequence_histogram", col = "#3f37c9", pch = 1, cex = 1.25,
     xlim = c(0,100), ylim = c(0,2.4)
)
points(x = as.numeric(names(repeat_sequence_perc_table_invert)), y = repeat_sequence_perc_table_invert, type = "o", cex = 1.25,
       col = "#f72585", pch = 1)
points(x = as.numeric(names(repeat_sequence_perc_table_plant)), y = repeat_sequence_perc_table_plant, type = "o", cex = 1.25,
       col = "#8ac926", pch = 1)

legend("topright",
       legend = c("Plant", "Invertebrate", "Chordate"),
       text.col = c("#8ac926", "#3f37c9", "#f72585"), bty = "n")

dev.off()













