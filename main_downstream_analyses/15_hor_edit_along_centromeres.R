bins = 30


ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}


####################################
taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 193
print(i)

setwd("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs/self_HOR_out_loose_settings") 

hor_dirs <- list.dirs(path = ".", recursive = FALSE)
hor_dirs_full <- list.dirs(path = ".", recursive = FALSE, full.names = T)


# i = 45 46 47 67 69 70

assembly <- strsplit(hor_dirs[i], split = "/")[[1]][2]
assembly <- strsplit(assembly, split = ".fa")[[1]][1]

setwd(hor_dirs[i])

repeats_files <- list.files(path = ".", pattern = "HOR_scored_repeats_with_hors")

total_HOR_score <- rep(0, bins)
total_repeats_per_bin <- rep(0, bins)

for(rep_file in repeats_files) {
  
  chr_repeats <- read.csv(file = rep_file)
  
  starts <- round(seq(0 , nrow(chr_repeats), length.out = bins + 1))
  ends <- starts[-1]
  starts <- starts[-length(starts)] + 1
  
  total_HOR_score <- unlist(lapply(seq_along(total_HOR_score), function(X) {
    total_HOR_score[X] + sum(chr_repeats$HOR_score[starts[X] : ends[X]])
  }))
  total_repeats_per_bin <- unlist(lapply(seq_along(total_HOR_score), function(X) {
    total_repeats_per_bin[X] + length(chr_repeats$HOR_score[starts[X] : ends[X]])
  }))
    
  
}

normalised_HOR_score <- total_HOR_score / total_repeats_per_bin 

normalised_HOR_score <- ma(c(normalised_HOR_score[1], normalised_HOR_score[1], normalised_HOR_score, normalised_HOR_score[length(normalised_HOR_score)], normalised_HOR_score[length(normalised_HOR_score)]))
normalised_HOR_score <- normalised_HOR_score[!is.na(normalised_HOR_score)]

pdf(file = paste0("/home/pwlodzimierz/ToL/upload_files/15_HOR_score_along/hor_score_along_", assembly, ".pdf"), 
    width = 7, height = 4)

plot(1 : bins, normalised_HOR_score, type = "l", lwd = 2, 
     ylab = "HOR score, %", xlab = paste0(bins, " bins"), main = paste0(" hor score along ", assembly, " chromosomes"))

dev.off()


pdf(file = paste0("/home/pwlodzimierz/ToL/upload_files/15_HOR_score_along/plots_same_ylim/hor_score_along_", assembly, ".pdf"), 
    width = 7, height = 4)

plot(1 : bins, normalised_HOR_score, type = "l", lwd = 2, 
     ylab = "HOR score, %", xlab = paste0(bins, " bins"), main = paste0(" hor score along ", assembly, " chromosomes"),
     ylim = c(0,8))

dev.off()


write.csv(x = normalised_HOR_score, file = paste0("/home/pwlodzimierz/ToL/upload_files/15_HOR_score_along/scores_csvs/scores_csvs_", assembly, ".csv"), row.names = FALSE)


























