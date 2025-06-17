#!/usr/bin/env Rscript
.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))
suppressMessages(library(seqinr))

source("./aux_fun.R")
replace_existing_analysis = FALSE
no_edta = FALSE # if TRUE, this assumes scripts 2 and 3 were skipped, no edta included, no cen templates recalculated



data_directories <- list.dirs(path = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs", recursive = FALSE, full.names = TRUE)
data_directories <- data_directories[!grepl(pattern = "templated_", data_directories)]
data_directories <- data_directories[grepl(pattern = ".fa", data_directories)]
assembly_files <- list.files(path = "/home/pwlodzimierz/ToL/Assemblies/fastas_2021_Michael", recursive = FALSE, full.names = TRUE)
assembly_files <- assembly_files[!grepl(".fai", assembly_files)]

taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 15
print(i)

# for(i in seq_along(data_directories)) 
{
  print(paste0(i, " / ", length(data_directories)))
  ### Load data ================================================================
  setwd(data_directories[i])
  print(data_directories[i])
  # TODO: restore the repeat and array load to templated directory
  # repeat_file = list.files(pattern = "_repeats_with_seq.csv", full.names = TRUE, recursive = TRUE)
  repeat_file = list.files(pattern = "_repeats_filtered.csv", full.names = TRUE, recursive = TRUE)
  if(length(repeat_file) != 1) {print(paste0(i, "No repeats!")); setwd(".."); next}
  #repeat_file = repeat_file[grep("template_run", repeat_file)]
  
  array_file = list.files(pattern = "_arrays_filtered.csv", full.names = TRUE, recursive = TRUE)
  if(length(array_file) != 1) {print(paste0(i, " no arrays!")); setwd(".."); next}
  #array_file = array_file[grep("template_run", array_file)]
  
  if(!no_edta) {
    templates_file = list.files(pattern = "_classes_merged_filtered.csv", full.names = TRUE, recursive = TRUE)
    if(length(templates_file) != 1) {print(paste0(i, " no templates!")); setwd(".."); next}
      
    edta_file = list.files(pattern = "_edta_filtered.csv", full.names = TRUE, recursive = TRUE)
    if(length(edta_file) != 1) {print(paste0(i, " no edta!")); setwd(".."); next}
    
    templates = read.csv(file = templates_file)
    names(templates)[1] = "class.Var1"
    names(templates)[3] = "template"
    templates$length = unlist(lapply(seq_along(templates$template), function(X) nchar(templates$template[X])))
    edta = read.csv(file = edta_file)
  } else {
    templates_file = list.files(pattern = "_classes_merged_filtered.csv", full.names = TRUE, recursive = TRUE)
    if(length(templates_file) != 1) {print(paste0(i, " no templates!")); setwd(".."); next}
    templates = read.csv(file = templates_file)
    names(templates)[1] = "class.Var1"
    names(templates)[3] = "template"
    templates$length = unlist(lapply(seq_along(templates$template), function(X) nchar(templates$template[X])))
  }
  
  
  assembly_name = strsplit(strsplit(data_directories[i], split = ".fa")[[1]][1], split = "v2_out_for_HORs/")[[1]][2]
  assembly_file = grep(strsplit(strsplit(repeat_file, split = "_repeats")[[1]][1], split = "/")[[1]][2], assembly_files)
  print(assembly_file)
  print(assembly_files[assembly_file])
  
  if(!replace_existing_analysis) {
    if(file.exists(paste0("genome_summary_", assembly_name, ".csv"))) quit(save = "no", status = 1)
  }
  
  if(length(assembly_file) == 1) {
    if(i != 178) {
      assembly = read.fasta(assembly_files[assembly_file], seqtype = "DNA", forceDNAtolower = TRUE)
    } 
    
  } else {
    print(" No assembly found, cannot proceed")
    setwd(def_wd)
    next
  }
  
  repeats = read.csv(file = repeat_file)
  arrays = read.csv(file = array_file)
  
  ### filter data to include only chromosomes
  
  print("filter chromosomes")
  
  all_chromosomes <- read.csv("/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full.Ian.csv")
  all_chromosomes <- all_chromosomes[all_chromosomes$assembly.name == strsplit(data_directories[i], split = "v2_out_for_HORs/")[[1]][2], ]
  all_chromosomes <- all_chromosomes[all_chromosomes$is.chr == 1, ]
  chromosomes <- all_chromosomes$chromosome.name
  chromosomes_lengths <- all_chromosomes$size
  
  bp_per_sequence_t = chromosomes_lengths
  bp_per_sequence = c(bp_per_sequence_t, 1000, 1000)
  bp_per_sequence = as.data.frame(bp_per_sequence)
  bp_per_sequence$chr_name <- c(chromosomes, "","")
  bp_per_sequence = bp_per_sequence[order(bp_per_sequence$bp_per_sequence, decreasing = TRUE), ]
  
  chromosomes = bp_per_sequence$chr_name
  print(chromosomes)
  
  repeats = repeats[repeats$seqID %in% chromosomes,]
  arrays = arrays[arrays$seqID %in% chromosomes,]
  if(!no_edta) edta = edta[edta$V1 %in% chromosomes,]
  ### chromosomes GC
  print("chr GC")
  if(i != 178) {
    gc_per_chr = unlist(lapply(seq_along(assembly), function(X) GC(assembly[[X]])))
    gc_per_chr_per_bp = unlist(lapply(seq_along(gc_per_chr), function(X) gc_per_chr[X] * bp_per_sequence_t[X] /  sum(bp_per_sequence_t)))
    genomic_GC = sum(gc_per_chr_per_bp)
  } else  {
    genomic_GC = 0
  }
  
  
  ### genome summary
  print("genome summary")

  genome_summary = data.frame(chromosomes = length(chromosomes),
                              genome_bp = sum(bp_per_sequence_t),
                              genome_GC = genomic_GC,
                              repeat_total_bp = sum(repeats$width),
                              repeat_fraction = sum(repeats$width) / sum(bp_per_sequence_t),
                              edta_total_bp = ifelse(!no_edta, sum(edta$width), 0),
                              edta_fraction = ifelse(!no_edta, (sum(edta$width) / sum(bp_per_sequence_t)), 0))
  
  
  if(!is.na(templates$class.Var1[1])) {
    repeat_1 = repeats[repeats$new_class == templates$class.Var1[1], ]
    genome_summary = cbind(genome_summary, data.frame(repeat_1_seq = templates$template[1],
                                                      repeat_1_bp = sum(repeat_1$width),
                                                      repeat_1_fraction = sum(repeat_1$width)/ sum(bp_per_sequence_t),
                                                      repeat_1_len = templates$length[1],
                                                      repeat_1_width_sd = sd(repeat_1$width),
                                                      repeat_1_mean_edit = mean(repeat_1$score_template),
                                                      repeat_1_edit_sd = sd(repeat_1$score_template),
                                                      repeat_1_HOR_per_rep2 = 0,
                                                      repeat_1_GC = GC(strsplit(paste0(repeat_1$sequence, collapse = ""), split = "")[[1]], forceToLower = TRUE)))
  } else {
    genome_summary = cbind(genome_summary, data.frame(repeat_1_seq = "",
                                                      repeat_1_bp = -1,
                                                      repeat_1_fraction = -1,
                                                      repeat_1_len = -1,
                                                      repeat_1_width_sd = -1,
                                                      repeat_1_mean_edit = -1,
                                                      repeat_1_edit_sd = -1,
                                                      repeat_1_HOR_per_rep2 = -1,
                                                      repeat_1_GC = -1))
  } 
  if(!is.na(templates$class.Var1[2])) {
    repeat_2 = repeats[repeats$new_class == templates$class.Var1[2], ]
    genome_summary = cbind(genome_summary, data.frame(repeat_2_seq = templates$template[2],
                                                      repeat_2_bp = sum(repeat_2$width),
                                                      repeat_2_fraction = sum(repeat_2$width)/ sum(bp_per_sequence_t),
                                                      repeat_2_len = templates$length[1],
                                                      repeat_2_width_sd = sd(repeat_2$width),
                                                      repeat_2_mean_edit = mean(repeat_2$score_template),
                                                      repeat_2_edit_sd = sd(repeat_2$score_template),
                                                      repeat_2_HOR_per_rep2 = 0,
                                                      repeat_2_GC = GC(strsplit(paste0(repeat_2$sequence, collapse = ""), split = "")[[1]], forceToLower = TRUE)))
  } else {
    genome_summary = cbind(genome_summary, data.frame(repeat_2_seq = "",
                                                      repeat_2_bp = -1,
                                                      repeat_2_fraction = -1,
                                                      repeat_2_len = -1,
                                                      repeat_2_width_sd = -1,
                                                      repeat_2_mean_edit = -1,
                                                      repeat_2_edit_sd = -1,
                                                      repeat_2_HOR_per_rep2 = -1,
                                                      repeat_2_GC = -1))
  }  
  if(!is.na(templates$class.Var1[3])) {
    repeat_3 = repeats[repeats$new_class == templates$class.Var1[3], ]
    genome_summary = cbind(genome_summary, data.frame(repeat_3_seq = templates$template[3],
                                                      repeat_3_bp = sum(repeat_3$width),
                                                      repeat_3_fraction = sum(repeat_3$width)/ sum(bp_per_sequence_t),
                                                      repeat_3_len = templates$length[3],
                                                      repeat_3_width_sd = sd(repeat_3$width),
                                                      repeat_3_mean_edit = mean(repeat_3$score_template),
                                                      repeat_3_edit_sd = sd(repeat_3$score_template),
                                                      repeat_3_HOR_per_rep2 = 0,
                                                      repeat_3_GC = GC(strsplit(paste0(repeat_3$sequence, collapse = ""), split = "")[[1]], forceToLower = TRUE)))
  } else {
    genome_summary = cbind(genome_summary, data.frame(repeat_3_seq = "",
                                                      repeat_3_bp = -1,
                                                      repeat_3_fraction = -1,
                                                      repeat_3_len = -1,
                                                      repeat_3_width_sd = -1,
                                                      repeat_3_mean_edit = -1,
                                                      repeat_3_edit_sd = -1,
                                                      repeat_3_HOR_per_rep2 = -1,
                                                      repeat_3_GC = -1))
  }  
  if(!is.na(templates$class.Var1[4])) {
    repeat_4 = repeats[repeats$new_class == templates$class.Var1[4], ]
    genome_summary = cbind(genome_summary, data.frame(repeat_4_seq = templates$template[4],
                                                      repeat_4_bp = sum(repeat_4$width),
                                                      repeat_4_fraction = sum(repeat_4$width)/ sum(bp_per_sequence_t),
                                                      repeat_4_len = templates$length[4],
                                                      repeat_4_width_sd = sd(repeat_4$width),
                                                      repeat_4_mean_edit = mean(repeat_4$score_template),
                                                      repeat_4_edit_sd = sd(repeat_4$score_template),
                                                      repeat_4_HOR_per_rep2 = 0,
                                                      repeat_4_GC = GC(strsplit(paste0(repeat_4$sequence, collapse = ""), split = "")[[1]], forceToLower = TRUE)))
  } else {
    genome_summary = cbind(genome_summary, data.frame(repeat_4_seq = "",
                                                      repeat_4_bp = -1,
                                                      repeat_4_fraction = -1,
                                                      repeat_4_len = -1,
                                                      repeat_4_width_sd = -1,
                                                      repeat_4_mean_edit = -1,
                                                      repeat_4_edit_sd = -1,
                                                      repeat_4_HOR_per_rep2 = -1,
                                                      repeat_4_GC = -1))
  }  
  if(!is.na(templates$class.Var1[5])) {
    repeat_5 = repeats[repeats$new_class == templates$class.Var1[5], ]
    genome_summary = cbind(genome_summary, data.frame(repeat_5_seq = templates$template[5],
                                                      repeat_5_bp = sum(repeat_5$width),
                                                      repeat_5_fraction = sum(repeat_5$width)/ sum(bp_per_sequence_t),
                                                      repeat_5_len = templates$length[5],
                                                      repeat_5_width_sd = sd(repeat_5$width),
                                                      repeat_5_mean_edit = mean(repeat_5$score_template),
                                                      repeat_5_edit_sd = sd(repeat_5$score_template),
                                                      repeat_5_HOR_per_rep2 = 0,
                                                      repeat_5_GC = GC(strsplit(paste0(repeat_5$sequence, collapse = ""), split = "")[[1]], forceToLower = TRUE)))
  } else {
    genome_summary = cbind(genome_summary, data.frame(repeat_5_seq = "",
                                                      repeat_5_bp = -1,
                                                      repeat_5_fraction = -1,
                                                      repeat_5_len = -1,
                                                      repeat_5_width_sd = -1,
                                                      repeat_5_mean_edit = -1,
                                                      repeat_5_edit_sd = -1,
                                                      repeat_5_HOR_per_rep2 = -1,
                                                      repeat_5_GC = -1))
  } 
  
  
  ### chromosome summary
  print("chromosome summary")
  chromosome_summary = NULL
  for(j in seq_along(chromosomes))
  {
    print(j)
    repeats_chr = repeats[repeats$seqID == chromosomes[j],]
    arrays_chr = arrays[arrays$seqID == chromosomes[j],]
    if(!no_edta) edta_chr = edta[edta$V1 == chromosomes[j],]
    
    if(i != 178) {
      which_gc <- which(chromosomes[j] %in% names(assembly))
      chromosome_summary_chr = data.frame(genome = assembly_name,
                                          chromosome = chromosomes[j],
                                          chromosome_bp = bp_per_sequence$bp_per_sequence[j],
                                          chromosome_GC = gc_per_chr[which_gc],
                                          repeat_total_bp = sum(repeats_chr$width),
                                          edta_total_bp = ifelse(!no_edta, sum(edta_chr$width), 0))
    } else {
      which_gc <- 1
      chromosome_summary_chr = data.frame(genome = assembly_name,
                                          chromosome = chromosomes[j],
                                          chromosome_bp = bp_per_sequence$bp_per_sequence[j],
                                          chromosome_GC = 0,
                                          repeat_total_bp = sum(repeats_chr$width),
                                          edta_total_bp = ifelse(!no_edta, sum(edta_chr$width), 0))
    }
    
    
    
    
    if(!is.na(templates$class.Var1[1]) & templates$class.Var1[1] %in% repeats_chr$new_class) {
      repeat_1 = repeats_chr[repeats_chr$new_class == templates$class.Var1[1], ]
      chromosome_summary_chr = cbind(chromosome_summary_chr, data.frame(repeat_1_seq = templates$template[1],
                                                                repeat_1_bp = sum(repeat_1$width),
                                                                repeat_1_fraction = sum(repeat_1$width) / bp_per_sequence$bp_per_sequence[j],
                                                                repeat_1_len = templates$length[1],
                                                                repeat_1_width_sd = sd(repeat_1$width),
                                                                repeat_1_mean_edit = mean(repeat_1$score_template),
                                                                repeat_1_edit_sd = sd(repeat_1$score_template),
                                                                repeat_1_HOR_per_rep2 = 0,
                                                                repeat_1_GC = GC(strsplit(paste0(repeat_1$sequence, collapse = ""), split = "")[[1]], forceToLower = TRUE),
                                                                repeat_1_mean_start_pos = mean(repeat_1$start),
                                                                repeat_1_sd_start_pos = sd(repeat_1$start)))
    } else {
      chromosome_summary_chr = cbind(chromosome_summary_chr, data.frame(repeat_1_seq = "",
                                                                repeat_1_bp = -1,
                                                                repeat_1_fraction = -1,
                                                                repeat_1_len = -1,
                                                                repeat_1_width_sd = -1,
                                                                repeat_1_mean_edit = -1,
                                                                repeat_1_edit_sd = -1,
                                                                repeat_1_HOR_per_rep2 = -1,
                                                                repeat_1_GC = -1,
                                                                repeat_1_mean_start_pos = -1,
                                                                repeat_1_sd_start_pos = -1))
    } 
    if(!is.na(templates$class.Var1[2]) & templates$class.Var1[2] %in% repeats_chr$new_class) {
      repeat_2 = repeats_chr[repeats_chr$new_class == templates$class.Var1[2], ]
      chromosome_summary_chr = cbind(chromosome_summary_chr, data.frame(repeat_2_seq = templates$template[2],
                                                                repeat_2_bp = sum(repeat_2$width),
                                                                repeat_2_fraction = sum(repeat_2$width) / bp_per_sequence$bp_per_sequence[j],
                                                                repeat_2_len = templates$length[2],
                                                                repeat_2_width_sd = sd(repeat_2$width),
                                                                repeat_2_mean_edit = mean(repeat_2$score_template),
                                                                repeat_2_edit_sd = sd(repeat_2$score_template),
                                                                repeat_2_HOR_per_rep2 = 0,
                                                                repeat_2_GC = GC(strsplit(paste0(repeat_2$sequence, collapse = ""), split = "")[[1]], forceToLower = TRUE),
                                                                repeat_2_mean_start_pos = mean(repeat_2$start),
                                                                repeat_2_sd_start_pos = sd(repeat_2$start)))
    } else {
      chromosome_summary_chr = cbind(chromosome_summary_chr, data.frame(repeat_2_seq = "",
                                                                repeat_2_bp = -1,
                                                                repeat_2_fraction = -1,
                                                                repeat_2_len = -1,
                                                                repeat_2_width_sd = -1,
                                                                repeat_2_mean_edit = -1,
                                                                repeat_2_edit_sd = -1,
                                                                repeat_2_HOR_per_rep2 = -1,
                                                                repeat_2_GC = -1,
                                                                repeat_2_mean_start_pos = -1,
                                                                repeat_2_sd_start_pos = -1))
    } 
    if(!is.na(templates$class.Var1[3]) & templates$class.Var1[3] %in% repeats_chr$new_class) {
      repeat_3 = repeats_chr[repeats_chr$new_class == templates$class.Var1[3], ]
      chromosome_summary_chr = cbind(chromosome_summary_chr, data.frame(repeat_3_seq = templates$template[3],
                                                                repeat_3_bp = sum(repeat_3$width),
                                                                repeat_3_fraction = sum(repeat_3$width) / bp_per_sequence$bp_per_sequence[j],
                                                                repeat_3_len = templates$length[3],
                                                                repeat_3_width_sd = sd(repeat_3$width),
                                                                repeat_3_mean_edit = mean(repeat_3$score_template),
                                                                repeat_3_edit_sd = sd(repeat_3$score_template),
                                                                repeat_3_HOR_per_rep2 = 0,
                                                                repeat_3_GC = GC(strsplit(paste0(repeat_3$sequence, collapse = ""), split = "")[[1]], forceToLower = TRUE),
                                                                repeat_3_mean_start_pos = mean(repeat_3$start),
                                                                repeat_3_sd_start_pos = sd(repeat_3$start)))
    } else {
      chromosome_summary_chr = cbind(chromosome_summary_chr, data.frame(repeat_3_seq = "",
                                                                repeat_3_bp = -1,
                                                                repeat_3_fraction = -1,
                                                                repeat_3_len = -1,
                                                                repeat_3_width_sd = -1,
                                                                repeat_3_mean_edit = -1,
                                                                repeat_3_edit_sd = -1,
                                                                repeat_3_HOR_per_rep2 = -1,
                                                                repeat_3_GC = -1,
                                                                repeat_3_mean_start_pos = -1,
                                                                repeat_3_sd_start_pos = -1))
    } 
    if(!is.na(templates$class.Var1[4]) & templates$class.Var1[4] %in% repeats_chr$new_class) {
      repeat_4 = repeats_chr[repeats_chr$new_class == templates$class.Var1[4], ]
      chromosome_summary_chr = cbind(chromosome_summary_chr, data.frame(repeat_4_seq = templates$template[4],
                                                                repeat_4_bp = sum(repeat_4$width),
                                                                repeat_4_fraction = sum(repeat_4$width) / bp_per_sequence$bp_per_sequence[j],
                                                                repeat_4_len = templates$length[4],
                                                                repeat_4_width_sd = sd(repeat_4$width),
                                                                repeat_4_mean_edit = mean(repeat_4$score_template),
                                                                repeat_4_edit_sd = sd(repeat_4$score_template),
                                                                repeat_4_HOR_per_rep2 = 0,
                                                                repeat_4_GC = GC(strsplit(paste0(repeat_4$sequence, collapse = ""), split = "")[[1]], forceToLower = TRUE),
                                                                repeat_4_mean_start_pos = mean(repeat_4$start),
                                                                repeat_4_sd_start_pos = sd(repeat_4$start)))
    } else {
      chromosome_summary_chr = cbind(chromosome_summary_chr, data.frame(repeat_4_seq = "",
                                                                repeat_4_bp = -1,
                                                                repeat_4_fraction = -1,
                                                                repeat_4_len = -1,
                                                                repeat_4_width_sd = -1,
                                                                repeat_4_mean_edit = -1,
                                                                repeat_4_edit_sd = -1,
                                                                repeat_4_HOR_per_rep2 = -1,
                                                                repeat_4_GC = -1,
                                                                repeat_4_mean_start_pos = -1,
                                                                repeat_4_sd_start_pos = -1))
    } 
    if(!is.na(templates$class.Var1[5]) & templates$class.Var1[5] %in% repeats_chr$new_class) {
      repeat_5 = repeats_chr[repeats_chr$new_class == templates$class.Var1[5], ]
      chromosome_summary_chr = cbind(chromosome_summary_chr, data.frame(repeat_5_seq = templates$template[5],
                                                                repeat_5_bp = sum(repeat_5$width),
                                                                repeat_5_fraction = sum(repeat_5$width) / bp_per_sequence$bp_per_sequence[j],
                                                                repeat_5_len = templates$length[5],
                                                                repeat_5_width_sd = sd(repeat_5$width),
                                                                repeat_5_mean_edit = mean(repeat_5$score_template),
                                                                repeat_5_edit_sd = sd(repeat_5$score_template),
                                                                repeat_5_HOR_per_rep2 = 0,
                                                                repeat_5_GC = GC(strsplit(paste0(repeat_5$sequence, collapse = ""), split = "")[[1]], forceToLower = TRUE),
                                                                repeat_5_mean_start_pos = mean(repeat_5$start),
                                                                repeat_5_sd_start_pos = sd(repeat_5$start)))
    } else {
      chromosome_summary_chr = cbind(chromosome_summary_chr, data.frame(repeat_5_seq = "",
                                                                repeat_5_bp = -1,
                                                                repeat_5_fraction = -1,
                                                                repeat_5_len = -1,
                                                                repeat_5_width_sd = -1,
                                                                repeat_5_mean_edit = -1,
                                                                repeat_5_edit_sd = -1,
                                                                repeat_5_HOR_per_rep2 = -1,
                                                                repeat_5_GC = -1,
                                                                repeat_5_mean_start_pos = -1,
                                                                repeat_5_sd_start_pos = -1))
    } 
    
  
    chromosome_summary = rbind(chromosome_summary, chromosome_summary_chr)
  }
  
  write.csv(genome_summary, paste0("genome_summary_", assembly_name, ".csv"), row.names = FALSE)
  write.csv(chromosome_summary, paste0("chromosomes_summary_", assembly_name, ".csv"), row.names = FALSE)
  
  write.csv(genome_summary, paste0("/home/pwlodzimierz/ToL/git_ToL/Summary_tables/genomes/genome_summary_", assembly_name, ".csv"), row.names = FALSE)
  write.csv(chromosome_summary, paste0("/home/pwlodzimierz/ToL/git_ToL/Summary_tables/chromosomes/chromosomes_summary_", assembly_name, ".csv"), row.names = FALSE)
  
  
  
  
  
  
}
