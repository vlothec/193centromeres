

similarity_tables_within <- list.files("/home/pwlodzimierz/ToL/upload_files/21_histograms_similarity_within_between_csvs", 
                                       pattern = "_within_pairs_df.csv", full.names = T)
similarity_tables_between <- list.files("/home/pwlodzimierz/ToL/upload_files/21_histograms_similarity_within_between_csvs", 
                                       pattern = "_between_pairs_df.csv", full.names = T)


### load data
data_directories <- list.dirs(path = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs", recursive = FALSE, full.names = TRUE)

main_table <- data.frame(assembly_short_name = vector(mode = "numeric", length = length(data_directories)),
                         similarity_within = vector(mode = "numeric", length = length(data_directories)),
                         similarity_between = vector(mode = "numeric", length = length(data_directories)),
                         hor_score = vector(mode = "numeric", length = length(data_directories)),
                         hor_score_SD = vector(mode = "numeric", length = length(data_directories)))


for(i in 1 : length(data_directories)) {
  print(i)
  
  setwd(data_directories[i])
  
  assembly_name_long <- strsplit(data_directories[i], split = "/")[[1]][7]
  assembly_name_short <- strsplit(assembly_name_long, split = ".fa")[[1]][1]
  
  main_table$assembly_short_name[i] = assembly_name_short
  
  
  hor_repeats_file <- list.files(path = paste0("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs/self_HOR_out_loose_settings/", assembly_name_long), 
                                 pattern = "HOR_scored_repeats_with_hors", full.names = TRUE)
  
  if(length(hor_repeats_file) == 0) next
  
  
  hor_repeats <- NULL
  for(j in seq_along(hor_repeats_file)) {
    hor_repeats <- rbind(hor_repeats, read.csv(file = hor_repeats_file[j]))
  }
  
  main_table$hor_score[i] = mean(hor_repeats$HOR_score)
  main_table$hor_score_SD[i] = sd(hor_repeats$HOR_score)
  
  which_sim_tables <- grep(assembly_name_short, similarity_tables_within)
  if(length(which_sim_tables) == 0) next
  
  main_table$similarity_within[i] = mean(read.csv(similarity_tables_within[which_sim_tables])$within_similarities)
  
  which_sim_tables <- grep(assembly_name_short, similarity_tables_between)
  if(length(which_sim_tables) == 0) next
  
  main_table$similarity_between[i] = mean(read.csv(similarity_tables_between[which_sim_tables])$between_similarities)
  
  
}


setwd("/home/pwlodzimierz/ToL/upload_files")

write.csv(main_table, file = "all_species_similarities_within_between_hor_scores.csv", row.names = FALSE)

main_table <- main_table[main_table$similarity_within != 0, ]

pdf("HOR_score_vs_similarity_within_minus_between.pdf")
plot(x = main_table$hor_score,
     y = main_table$similarity_within - main_table$similarity_between,
     xlab = "Mean HOR score",
     ylab = "Mean similarity within minus between",
     pch = 16)
dev.off()

pdf("HOR_score_vs_similarity_between.pdf")
plot(x = main_table$hor_score,
     y = main_table$similarity_between,
     xlab = "Mean HOR score",
     ylab = "Mean similarity between",
     pch = 16)
dev.off()

pdf("HOR_score_vs_similarity_within.pdf")
plot(x = main_table$hor_score,
     y = main_table$similarity_within,
     xlab = "Mean HOR score",
     ylab = "Mean similarity within",
     pch = 16)
dev.off()

pdf("hits_similarities_within_minus_between.pdf")
hist(main_table$similarity_within - main_table$similarity_between, breaks = 20,
     xlab = "Mean similarity within minus between",
     main = "Mean similarity within minus between",
     pch = 16)
dev.off()


pdf("HOR_score_SD_vs_HOR_score.pdf")
plot(x = main_table$hor_score,
     y = main_table$hor_score_SD / main_table$hor_score,
     xlab = "Mean HOR score",
     ylab = "Mean HOR score SD in score %",
     pch = 16)
dev.off()













