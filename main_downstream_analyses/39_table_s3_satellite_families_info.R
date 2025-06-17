


satellite_metadata <- read.csv("/home/pwlodzimierz/ToL/curated_satellites_repDec24_jan2025_may25.csv")

satellite_metadata$mean_repeat_size <- 0
satellite_metadata$repeat_no <- 0
satellite_metadata$repeat_bp <- 0
satellite_metadata$repeat_size_SD_in_perc <- 0
satellite_metadata$consensus <- NA
satellite_metadata$mean_similairty_within_chromosomes <- NA
satellite_metadata$mean_similairty_between_chromosomes <- NA
satellite_metadata$holocentrics.mean_similairty_between_arrays <- NA
satellite_metadata$holocentrics.mean_array_size <- NA

for(i in 1 : nrow(satellite_metadata)) {
  cat(i, satellite_metadata$Species[i], "\n")
  
  setwd(paste0("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs/", satellite_metadata$Genome[i]))
  
  repeats <- read.csv(file = paste0("./", satellite_metadata$Genome[i], "_repeats_filtered_array_assigned.csv"))
  
  
  trash_names <- strsplit(satellite_metadata$TRASH_name_dec2024runs[i], split = ";")[[1]]
  
  repeats <- repeats[repeats$new_class %in% trash_names, ]
  
  # this only includes chromosome-only repeats!!!
  
  satellite_metadata$mean_repeat_size[i] <- mean(repeats$width)
  satellite_metadata$repeat_no[i] <- length(repeats$width)
  satellite_metadata$repeat_bp[i] <- sum(repeats$width)
  satellite_metadata$repeat_size_SD_in_perc[i] <- 100 * sd(repeats$width) / mean(repeats$width)
  
  
  
}

satellite_metadata$Satellite_name = ""
for(i in 1 : nrow(satellite_metadata)) {
  cat(i, satellite_metadata$Species[i], "\n")
  
  shorter_name <- strsplit(satellite_metadata$Genome[i], split = "[.]")[[1]][1]
  remove_last_char <- paste0(strsplit(shorter_name, split = "")[[1]][1 : (nchar(shorter_name) - 1)], collapse = "")
  
  satellite_metadata$Satellite_name[i] <- paste0(remove_last_char, ".", round(satellite_metadata$mean_repeat_size[i]))
  
  
  
}

satellite_metadata$repeat_no_incl_non_chr_seq <- 0

for(i in 1 : nrow(satellite_metadata)) {
  cat(i, satellite_metadata$Species[i], "\n")
  
  setwd(paste0("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs/", satellite_metadata$Genome[i]))
  
  short_name <- strsplit(satellite_metadata$Genome[i], split = ".fa")[[1]][1]
  
  repeats <- read.csv(file = paste0("./", short_name, "_repeats_filtered.csv"))
  
  
  trash_names <- strsplit(satellite_metadata$TRASH_name_dec2024runs[i], split = ";")[[1]]
  
  repeats <- repeats[repeats$new_class %in% trash_names, ]
  
  satellite_metadata$repeat_no_incl_non_chr_seq[i] <- length(repeats$width)
  
}

names(satellite_metadata)[names(satellite_metadata) == "repeat_no_incl_non_chr_seq"] <- "repeat_no.including_non_chr_seq"


satellite_metadata$Satellite_name[satellite_metadata$Satellite_name == "rosCan_S27_v.141"] = "drRosCani.141"

write.csv(x = satellite_metadata, file = "/home/pwlodzimierz/ToL/curated_satellites_metadata_on_chromosomes_only_may.csv")





