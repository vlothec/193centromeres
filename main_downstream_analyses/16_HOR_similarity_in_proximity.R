


samples <- 500000
max_dist_away <- 16

ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}



####################################

setwd("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/HORs/self_HOR_out_loose_settings") 

hor_dirs <- list.dirs(path = ".", recursive = FALSE)
hor_dirs_full <- list.dirs(path = ".", recursive = FALSE, full.names = T)

taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 193
print(i)

assembly <- strsplit(hor_dirs[i], split = "/")[[1]][2]
assembly <- strsplit(assembly, split = ".fa")[[1]][1]

setwd(hor_dirs[i])

repeats_files <- list.files(path = ".", pattern = "HOR_scored_repeats_with_hors")

repeats <- NULL

for(rep_file in repeats_files) repeats <- rbind(repeats, read.csv(file = rep_file))

sample_reps_a <- sample(max_dist_away : nrow(repeats), size = samples, replace = T)

sample_reps_b <- sample_reps_a + sample(c(-max_dist_away : -1, 1 : max_dist_away), 1)

sample_reps_b[sample_reps_b < 1] <- 1

similarity_values <- unlist(lapply(seq_along(sample_reps_a), function(X) {
  100 * adist(repeats$sequence[sample_reps_a[X]], repeats$sequence[sample_reps_b[X]])[[1]] / repeats$width[sample_reps_a[X]]
}))

similarity_values[similarity_values > 100] <- 100
similarity_values <- similarity_values[similarity_values >= 0]
similarity_values <- 100 - similarity_values

pdf(file = paste0("/home/pwlodzimierz/ToL/upload_files/similarity_in_proximity_histograms/hist_similarity_in_proximity_", assembly, ".pdf"), 
    width = 7, height = 4)

hist_data <- hist(similarity_values, breaks = seq(0, 100, length.out = 150), plot = F)
hist(similarity_values, breaks = seq(0, 100, length.out = 150), ylim = c(0, max(hist_data$counts)), 
     main = paste0(assembly, nrow(repeats)))
new_counts <- ma(c(hist_data$counts[1], hist_data$counts[1], hist_data$counts, hist_data$counts[length(hist_data$counts)], hist_data$counts[length(hist_data$counts)]))
new_counts <- new_counts[!is.na(new_counts)]

par(new = TRUE)
plot(1 : length(new_counts), new_counts, type = "l", lwd = 2, 
     xlab = "", ylab = "", ylim = c(0, max(hist_data$counts)), xaxt = "n", yaxt = "n", main = "")

dev.off()

write.csv(x = new_counts, file = paste0("/home/pwlodzimierz/ToL/upload_files/similarity_in_proximity_histograms/histogram_counts/histogram_counts_", assembly, ".csv"), row.names = FALSE, col.names = FALSE)


























