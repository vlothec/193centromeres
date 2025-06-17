

taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 193
print(i)



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

print(paste0(i, " / ", nrow(satellite_metadata)))

hors_monomer_length_table <- NULL
hors_distances_table_kbp <- NULL

# for(i in 1 : nrow(satellite_metadata)) { # 1 to 91
setwd(paste0("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs/self_HOR_out_loose_settings/", satellite_metadata$Genome[i]))
  
  if(file.exists(paste0("hors_monomer_length_table_", satellite_metadata$Genome[i], ".csv"))) {
    file.remove(paste0("hors_monomer_length_table_", satellite_metadata$Genome[i], ".csv"))
    file.remove(paste0("hors_distances_table_kbp_", satellite_metadata$Genome[i], ".csv"))
  }

hor_files <- list.files(path = ".", pattern = "HORs_")
hor_files <- hor_files[grepl(".csv", hor_files)]
hor_files <- hor_files[!grepl("grand", hor_files)]

for(j in seq_along(hor_files)) {
  cat(i, "/", nrow(satellite_metadata), "", j, "/", length(hor_files), "\n")
  hors <- read.csv(file = hor_files[j])
  mono_temp_table <- sort(table(hors$block.size.in.units), decreasing = T)
  
  hors_monomer_length_table <- c(hors_monomer_length_table, mono_temp_table[setdiff(names(mono_temp_table), names(hors_monomer_length_table))])
  hors_monomer_length_table[names(mono_temp_table)] <- hors_monomer_length_table[names(mono_temp_table)] + mono_temp_table
  
  dist_temp_table <- sort(table(round(abs(hors$start.B.bp - hors$start.A.bp)/1000)), decreasing = T)
  hors_distances_table_kbp <- c(hors_distances_table_kbp, dist_temp_table[setdiff(names(dist_temp_table), names(hors_distances_table_kbp))])
  hors_distances_table_kbp[names(dist_temp_table)] <- hors_distances_table_kbp[names(dist_temp_table)] + dist_temp_table
}

write.csv(hors_monomer_length_table, file = paste0("hors_monomer_length_table_", satellite_metadata$Genome[i], ".csv"), row.names = T)

write.csv(hors_distances_table_kbp, file = paste0("hors_distances_table_kbp_", satellite_metadata$Genome[i], ".csv"), row.names = T)


# }






if(FALSE) ### Do this manually after all the scripts above finished
{
  satellite_metadata <- read.csv("/home/pwlodzimierz/ToL/curated_satellites_metadata_on_chromosomes_only_may.csv")
  
  hors_monomer_length_table_vert <- NULL # 5h. RANK, row count, divide by number of vertebrate species to get mean value per species
  hors_monomer_length_perc_table_vert <- NULL # 5g. Block lengths in monomer, divide by number of vertebrate species to make sure the vector adds up to 1
  hors_distances_table_perc_kbp_vert <- NULL # 5i. Distances, percentages, divide by number of vertebrate species to make sure the vector adds up to 1
  
  hors_monomer_length_table_invert <- NULL
  hors_monomer_length_perc_table_invert <- NULL
  hors_distances_table_perc_kbp_invert <- NULL
  
  hors_monomer_length_table_plant <- NULL
  hors_monomer_length_perc_table_plant <- NULL
  hors_distances_table_perc_kbp_plant <- NULL
  
  for(i in 1 : nrow(satellite_metadata)) { # 1 to 91
    setwd(paste0("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs/self_HOR_out_loose_settings/", satellite_metadata$Genome[i]))
    
    
    mono_temp_table_file <- list.files(path = ".", pattern = "hors_monomer_length_table_")
    if(length(mono_temp_table_file) == 0) {
      cat(i, "No file from", satellite_metadata$Genome[i], "\n")
      next
    }
    mono_temp_table <- read.csv(file = mono_temp_table_file, row.names = 1)
    spare_names_mono <- rownames(mono_temp_table)
    row.names_mono <- 1 : nrow(mono_temp_table)
    mono_temp_table <- as.vector(t(mono_temp_table))
    names(mono_temp_table) <- row.names_mono
    mono_temp_table_perc <- 100 * mono_temp_table / sum(mono_temp_table)
    # format
    dist_temp_table_file <- list.files(path = ".", pattern = "hors_distances_table_kbp_")
    dist_temp_table <- read.csv(file = dist_temp_table_file, row.names = 1)
    row.names_dist <- rownames(dist_temp_table)
    dist_temp_table <- as.vector(t(dist_temp_table))
    names(dist_temp_table) <- row.names_dist
    dist_temp_table_perc <- 100 * dist_temp_table / sum(dist_temp_table)
    # format
    
    
    if(satellite_metadata$group2[i] == "vertebrate") {
      hors_monomer_length_table_vert <- c(hors_monomer_length_table_vert, mono_temp_table[setdiff(names(mono_temp_table), names(hors_monomer_length_table_vert))])
      hors_monomer_length_table_vert[names(mono_temp_table)] <- hors_monomer_length_table_vert[names(mono_temp_table)] + mono_temp_table
      
      hors_monomer_length_perc_table_vert <- c(hors_monomer_length_perc_table_vert, mono_temp_table_perc[setdiff(spare_names_mono, names(hors_monomer_length_perc_table_vert))])
      hors_monomer_length_perc_table_vert[spare_names_mono] <- hors_monomer_length_perc_table_vert[spare_names_mono] + mono_temp_table_perc
      
      hors_distances_table_perc_kbp_vert <- c(hors_distances_table_perc_kbp_vert, dist_temp_table_perc[setdiff(names(dist_temp_table_perc), names(hors_distances_table_perc_kbp_vert))])
      hors_distances_table_perc_kbp_vert[names(dist_temp_table_perc)] <- hors_distances_table_perc_kbp_vert[names(dist_temp_table_perc)] + dist_temp_table_perc
    }
    if(satellite_metadata$group2[i] == "invertebrate") {
      hors_monomer_length_table_invert <- c(hors_monomer_length_table_invert, mono_temp_table[setdiff(names(mono_temp_table), names(hors_monomer_length_table_invert))])
      hors_monomer_length_table_invert[names(mono_temp_table)] <- hors_monomer_length_table_invert[names(mono_temp_table)] + mono_temp_table
      
      hors_monomer_length_perc_table_invert <- c(hors_monomer_length_perc_table_invert, mono_temp_table_perc[setdiff(spare_names_mono, names(hors_monomer_length_perc_table_invert))])
      hors_monomer_length_perc_table_invert[spare_names_mono] <- hors_monomer_length_perc_table_invert[spare_names_mono] + mono_temp_table_perc
      
      hors_distances_table_perc_kbp_invert <- c(hors_distances_table_perc_kbp_invert, dist_temp_table_perc[setdiff(names(dist_temp_table_perc), names(hors_distances_table_perc_kbp_invert))])
      hors_distances_table_perc_kbp_invert[names(dist_temp_table_perc)] <- hors_distances_table_perc_kbp_invert[names(dist_temp_table_perc)] + dist_temp_table_perc
    }
    if(satellite_metadata$group2[i] == "plant") {
      hors_monomer_length_table_plant <- c(hors_monomer_length_table_plant, mono_temp_table[setdiff(names(mono_temp_table), names(hors_monomer_length_table_plant))])
      hors_monomer_length_table_plant[names(mono_temp_table)] <- hors_monomer_length_table_plant[names(mono_temp_table)] + mono_temp_table
      
      hors_monomer_length_perc_table_plant <- c(hors_monomer_length_perc_table_plant, mono_temp_table_perc[setdiff(spare_names_mono, names(hors_monomer_length_perc_table_plant))])
      hors_monomer_length_perc_table_plant[spare_names_mono] <- hors_monomer_length_perc_table_plant[spare_names_mono] + mono_temp_table_perc
      
      hors_distances_table_perc_kbp_plant <- c(hors_distances_table_perc_kbp_plant, dist_temp_table_perc[setdiff(names(dist_temp_table_perc), names(hors_distances_table_perc_kbp_plant))])
      hors_distances_table_perc_kbp_plant[names(dist_temp_table_perc)] <- hors_distances_table_perc_kbp_plant[names(dist_temp_table_perc)] + dist_temp_table_perc
    }
    
    
    
    
  }
  
  # normalise by the count
  hors_monomer_length_table_vert <- hors_monomer_length_table_vert / nrow(satellite_metadata[satellite_metadata$group2 == "vertebrate",])
  hors_monomer_length_perc_table_vert <- hors_monomer_length_perc_table_vert / nrow(satellite_metadata[satellite_metadata$group2 == "vertebrate",])
  hors_distances_table_perc_kbp_vert  <- hors_distances_table_perc_kbp_vert / nrow(satellite_metadata[satellite_metadata$group2 == "vertebrate",])
  
  hors_monomer_length_table_invert <- hors_monomer_length_table_invert / nrow(satellite_metadata[satellite_metadata$group2 == "invertebrate",])
  hors_monomer_length_perc_table_invert <- hors_monomer_length_perc_table_invert / nrow(satellite_metadata[satellite_metadata$group2 == "invertebrate",])
  hors_distances_table_perc_kbp_invert <- hors_distances_table_perc_kbp_invert / nrow(satellite_metadata[satellite_metadata$group2 == "invertebrate",])
  
  hors_monomer_length_table_plant <- hors_monomer_length_table_plant / nrow(satellite_metadata[satellite_metadata$group2 == "plant",])
  hors_monomer_length_perc_table_plant  <- hors_monomer_length_perc_table_plant / nrow(satellite_metadata[satellite_metadata$group2 == "plant",])
  hors_distances_table_perc_kbp_plant <- hors_distances_table_perc_kbp_plant / nrow(satellite_metadata[satellite_metadata$group2 == "plant",])
  
  
  # clean up NA values
  hors_monomer_length_table_vert <- hors_monomer_length_table_vert[!is.na(hors_monomer_length_table_vert)]
  hors_monomer_length_perc_table_vert <- hors_monomer_length_perc_table_vert[!is.na(hors_monomer_length_perc_table_vert)]
  hors_distances_table_perc_kbp_vert  <- hors_distances_table_perc_kbp_vert[!is.na(hors_distances_table_perc_kbp_vert)]
  
  hors_monomer_length_table_invert <- hors_monomer_length_table_invert[!is.na(hors_monomer_length_table_invert)]
  hors_monomer_length_perc_table_invert <- hors_monomer_length_perc_table_invert[!is.na(hors_monomer_length_perc_table_invert)]
  hors_distances_table_perc_kbp_invert <- hors_distances_table_perc_kbp_invert[!is.na(hors_distances_table_perc_kbp_invert)]
  
  hors_monomer_length_table_plant <- hors_monomer_length_table_plant[!is.na(hors_monomer_length_table_plant)]
  hors_monomer_length_perc_table_plant  <- hors_monomer_length_perc_table_plant[!is.na(hors_monomer_length_perc_table_plant)]
  hors_distances_table_perc_kbp_plant <- hors_distances_table_perc_kbp_plant[!is.na(hors_distances_table_perc_kbp_plant)]
  
  
  hors_distances_table_perc_kbp_vert <- hors_distances_table_perc_kbp_vert[order(as.numeric(names(hors_distances_table_perc_kbp_vert)))]
  hors_distances_table_perc_kbp_invert <- hors_distances_table_perc_kbp_invert[order(as.numeric(names(hors_distances_table_perc_kbp_invert)))]
  hors_distances_table_perc_kbp_plant <- hors_distances_table_perc_kbp_plant[order(as.numeric(names(hors_distances_table_perc_kbp_plant)))]
  
  hors_distances_table_perc_kbp_vert <- hors_distances_table_perc_kbp_vert[-1]
  hors_distances_table_perc_kbp_invert <- hors_distances_table_perc_kbp_invert[-1]
  hors_distances_table_perc_kbp_plant <- hors_distances_table_perc_kbp_plant[-1]
  
  setwd("/home/pwlodzimierz/ToL/upload_files/41_hor_data/")
  
  # - HOR monomer length rank analysis
  freq <- sort(as.numeric(hors_monomer_length_table_vert), decreasing = TRUE)
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
  
  freq <- sort(as.numeric(hors_monomer_length_table_invert), decreasing = TRUE)
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
  
  freq <- sort(as.numeric(hors_monomer_length_table_plant), decreasing = TRUE)
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
  
  summary(vert_fit)$r.squared
  summary(invert_fit)$r.squared
  summary(plant_fit)$r.squared
  
  # png("./HOR_monomer_length_rank_analysis_ln.png", width = 900, height = 900)
  pdf("./HOR_monomer_length_rank_analysis_ln.pdf", width = 5, height = 5)
  
  plot(invert_data$log.rank, invert_data$log.freq, type = "o", 
       xlab = "log(Rank)", ylab = "log(Frequency)", main = "HOR_monomer_length_rank_analysis_ln", col = "#3f37c9", pch = 1, cex = 1.25,
       xlim = c(0,8), ylim = c(-2,15))
  abline(invert_fit, col = "#3f37c9", lty = 2)
  
  points(vert_data$log.rank, vert_data$log.freq, type = "o", 
         col = "#f72585", cex = 1.25, pch = 1)
  abline(vert_fit, col = "#f72585", lty = 2)
  
  points(plant_data$log.rank, plant_data$log.freq, type = "o", 
         col = "#8ac926", cex = 1.25, pch = 1)
  abline(plant_fit, col = "#8ac926", lty = 2)
  
  legend("topright",
         legend = c("Plant", "Invertebrate", "Chordate"),
         text.col = c("#8ac926", "#3f37c9", "#f72585"), bty = "n")
  
  dev.off()
  
  
  
  # - HOR monomer length histogram (use % instead of frequency) 
  
  pdf("./HOR_monomer_length_histogram.pdf", width = 5, height = 5)
  
  plot(x = as.numeric(names(hors_monomer_length_perc_table_vert)), y = hors_monomer_length_perc_table_vert, type = "o", 
       xlab = "HOR monomer length", ylab = "% of HORs", main = "HOR_monomer_length_histogram", col = "#3f37c9", pch = 1, cex = 1.25,
       xlim = c(4,34), ylim = c(0,45)
  )
  points(x = as.numeric(names(hors_monomer_length_perc_table_invert)), y = hors_monomer_length_perc_table_invert, type = "o", cex = 1.25,
         col = "#f72585", pch = 1)
  points(x = as.numeric(names(hors_monomer_length_perc_table_plant)), y = hors_monomer_length_perc_table_plant, type = "o", cex = 1.25,
         col = "#8ac926", pch = 1)
  legend("topright",
         legend = c("Plant", "Invertebrate", "Chordate"),
         text.col = c("#8ac926", "#3f37c9", "#f72585"), bty = "n")
  
  dev.off()
  
  
  
  # - inter-HOR distances histogram (use % instead of frequency)
  
  
  pdf("./HOR_block_distances_histogram.pdf", width = 5, height = 5)
  
  plot(x = as.numeric(names(hors_distances_table_perc_kbp_vert))/1000, y = hors_distances_table_perc_kbp_vert, type = "o", 
       xlab = "Inter-HOR distances (Mb)", ylab = "% of HORs", main = "HOR_block_distances_histogram", col = "#3f37c9", pch = 1, cex = 1.25,
       xlim = c(0,1.5), ylim = c(0,2.3)
       )
  points(x = as.numeric(names(hors_distances_table_perc_kbp_invert))/1000, y = hors_distances_table_perc_kbp_invert, type = "o", cex = 1.25,
         col = "#f72585", pch = 1)
  points(x = as.numeric(names(hors_distances_table_perc_kbp_plant))/1000, y = hors_distances_table_perc_kbp_plant, type = "o", cex = 1.25,
         col = "#8ac926", pch = 1)
  legend("topright",
         legend = c("Plant", "Invertebrate", "Chordate"),
         text.col = c("#8ac926", "#3f37c9", "#f72585"), bty = "n")
 
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}











