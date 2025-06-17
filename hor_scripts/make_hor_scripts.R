setwd("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs")

satellites_to_run <- read.csv(file = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/curated_satellites_repDec24_jan2025.csv")
chromosomes_table <- read.csv(file = "/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full.Ian.csv")

# out_dir="$1"
# basename="repeats_with_hors_"
# class="$3"
# chr="$2"
# end=".csv"

commands <- NULL
commands_first2 <- NULL
for(i in seq_len(nrow(satellites_to_run))) {
  genome_name <- satellites_to_run$Genome[i]
  
  cat(genome_name, "\n")
  
  repeat_names <- strsplit(satellites_to_run$TRASH_name_dec2024runs[i], split = ";")[[1]]
  
  
  repeats_file <- list.files(path = paste0("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs/", genome_name), pattern = "_repeats_filtered.csv", full.names = TRUE)
  
  if(length(repeats_file) != 1) stop(paste0("No repeats could be found for genome_name"))
  
  chromosomes <- chromosomes_table$chromosome.name[which(chromosomes_table$assembly.name == genome_name & chromosomes_table$is.chr == TRUE)]
  
  for(j in seq_along(repeat_names)) {
    for(k in seq_along(chromosomes)) {
      
      out_dir <- paste0("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs/self_HOR_out/", genome_name, "/")
      
      command <- paste("./hors_slurm.sh", out_dir, chromosomes[k], repeat_names[j], repeats_file, paste = " ")
      commands <- c(commands, command)
      
      if(k > 2 & k < 7) {
        commands_first2 <- c(commands_first2, command)
      }
      
    }
    
    
  }
  
}

write.table(x = commands, file = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs/commands_same_chr.txt",
            quote = F, sep = "\n", append = FALSE, col.names = FALSE, row.names = FALSE)

write.table(x = commands_first2, file = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs/commands_same_chr_first2chr.txt",
            quote = F, sep = "\n", append = FALSE, col.names = FALSE, row.names = FALSE)



#####
#
# Make new commands only for chromosomes that didn't finish previously
#
#####



setwd("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs")

satellites_to_run <- read.csv(file = "/home/pwlodzimierz/ToL/curated_satellites_repDec24_jan2025_may25.csv")
chromosomes_table <- read.csv(file = "/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full.Ian.csv")

# out_dir="$1"
# basename="repeats_with_hors_"
# class="$3"
# chr="$2"
# end=".csv"

commands <- NULL
commands_first2 <- NULL
for(i in seq_len(nrow(satellites_to_run))) {
  genome_name <- satellites_to_run$Genome[i]
  
  cat(genome_name, "\n")
  
  repeat_names <- strsplit(satellites_to_run$TRASH_name_dec2024runs[i], split = ";")[[1]]
  
  
  repeats_file <- list.files(path = paste0("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs/", genome_name), pattern = "_repeats_filtered.csv", full.names = TRUE)
  
  if(length(repeats_file) != 1) stop(paste0("No repeats could be found for genome_name"))
  
  chromosomes <- chromosomes_table$chromosome.name[which(chromosomes_table$assembly.name == genome_name & chromosomes_table$is.chr == TRUE)]
  
  for(j in seq_along(repeat_names)) {
    for(k in seq_along(chromosomes)) {
      
      
      out_dir <- paste0("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs/self_HOR_out_loose_settings/", genome_name, "/")
      
      if(file.exists(paste0(out_dir, "HORs_lines_", repeat_names[j], "_", chromosomes[k], ".png"))) next
      
      command <- paste("./hors_slurm_loose_settings.sh", out_dir, chromosomes[k], repeat_names[j], repeats_file, paste = " ")
      commands <- c(commands, command)
      
      if(k > 0 & k < 3) {
        commands_first2 <- c(commands_first2, command)
      }
      
    }
    
    
  }
  
}

write.table(x = commands, file = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs/loose_commands_same_chr_fill_undone.txt",
            quote = F, sep = "\n", append = FALSE, col.names = FALSE, row.names = FALSE)
# 
# write.table(x = commands_first2, file = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs/loose_commands_same_chr_first2chr_fill_undone.txt",
#             quote = F, sep = "\n", append = FALSE, col.names = FALSE, row.names = FALSE)


















