

similarity_files <- list.files(path = "/home/pwlodzimierz/ToL/upload_files/37_similarity_values_within_between_chr", full.names = T)


genome_summary <- data.frame(genome_name = vector(mode = "character"),
                             mean_sim_within = vector(mode = "numeric"),
                             mean_sim_between = vector(mode = "numeric"),
                             p.val = vector(mode = "numeric"),
                             no_of_chr = vector(mode = "numeric"))

for(i in 1 : length(similarity_files)) {
  print(i)
  
  sim_data <- read.csv(file = similarity_files[i])
  
  chromosomes <- sim_data$similarity[sim_data$method == "same_chr"]
  genomes     <- sim_data$similarity[sim_data$method == "diff_chr"]
  
  observed_diff <- mean(chromosomes) - mean(genomes)
  combined <- c(chromosomes, genomes)
  group_labels <- c(rep(1, length(chromosomes)), rep(2, length(genomes)))
  
  perm_diffs <- replicate(1000, {
    shuffled_labels <- sample(group_labels)
    mean(combined[shuffled_labels == 1]) - mean(combined[shuffled_labels == 2])
  })
  
  p_value <- mean(abs(perm_diffs) >= abs(observed_diff))
                    

  
  genome_summary <- rbind(genome_summary, data.frame(genome_name = strsplit(strsplit(similarity_files[i], split = "/similarities_within_between_chrs_")[[1]][2], split = "_")[[1]][1],
                                                     mean_sim_within = mean(chromosomes),
                                                     mean_sim_between = mean(genomes),
                                                     # p.val = wilcox.test(chromosomes, genomes)$p.value,
                                                     p.val = p_value,
                                                     no_of_chr = length(unique(c(sim_data$chr1, sim_data$chr2)))))
  
}
satellite_metadata <- read.csv("/home/pwlodzimierz/ToL/curated_satellites_repDec24_jan2025_may25.csv")

# genome_summary$mean_sim_within - genome_summary$mean_sim_between

genome_summary$Clade <- ""

for(i in 1 : nrow(genome_summary)) {
  
  genome_summary$Clade[i] = satellite_metadata$Clade[satellite_metadata$Genome == genome_summary$genome_name[i]][1]
  
}
genome_summary$Clade[genome_summary$genome_name == "rosCan"] = "Dicot"

genome_summary$group <- "Invertebrate"
genome_summary$group[genome_summary$Clade %in% c("Dicot", "Monocot")] <- "Plant"
genome_summary$group[genome_summary$Clade %in% c("Aves", "Fish", "Chordata", "Mammalia", "Reptilia")] <- "Chordate"

pdf(file = "/home/pwlodzimierz/ToL/upload_files/38_similarity_between_within_summary_scatters/38_scatter_simialrity_within_vs_between_chromosomes.pdf")
# Create a point symbol vector based on p-value thresholds
point_pch <- ifelse(genome_summary$p.val < 0.001, 15,
                    ifelse(genome_summary$p.val < 0.01, 17,
                           ifelse(genome_summary$p.val < 0.05, 16, 1)))

# Create a color vector based on group
point_colors <- ifelse(genome_summary$group == "Plant", "#8ac926",
                       ifelse(genome_summary$group == "Invertebrate", "#3f37c9", "#f72585"))

# Plot with p-value-based symbols and group-based colors
plot(x = genome_summary$mean_sim_within,
     y = genome_summary$mean_sim_between,
     pch = point_pch,
     col = point_colors,
     xlab = "Mean Similarity Within Chromosomes",
     ylab = "Mean Similarity Between Chromosomes",
     main = "Genome Similarity Scatter Plot by Group",
     cex = 1.2,
     xlim = c(40,100), ylim = c(40,100))

# Add a legend for group colors
legend("topleft", 
       legend = c("Plant", "Invertebrate", "Chordate"),
       col = c("#8ac926", "#3f37c9", "#f72585"),
       pch = 16,
       title = "Group",
       cex = 0.8)

# Add a legend for p-value symbols
legend("bottomright", 
       legend = c("p ≥ 0.05", "p < 0.05", "p < 0.01", "p < 0.001"),
       pch = c(1, 16, 17, 15),
       col = "black",
       title = "P-value",
       cex = 0.8)

lines(x = c(40,100), y = c(40,100), lty = "aa")
dev.off()



genome_summary$diff_within_minus_between <- genome_summary$mean_sim_within - genome_summary$mean_sim_between


pdf(file = "/home/pwlodzimierz/ToL/upload_files/38_similarity_between_within_summary_scatters/38_scatter_relative_similarity_vs_chr_number.pdf")
plot(x = genome_summary$no_of_chr, 
     y = genome_summary$diff_within_minus_between,
     xlab = "number of chromosomes",
     ylab = "Similarity within minus similarity between chromosomes",
     col = point_colors,
     pch = 16)
text(x = min(genome_summary$no_of_chr) + max(genome_summary$no_of_chr)/100, y = 30, 
     paste0("p=",round(cor.test(genome_summary$no_of_chr, genome_summary$diff_within_minus_between)$p.value,4)), adj = 0)
dev.off()


genome_organisation_metadata <- read.csv("/home/pwlodzimierz/ToL/genomes_organisation_type.csv")

genome_summary$genome_size <- 0
for(i in 1 : nrow(genome_summary)) {
  
  genome_summary$genome_size[i] = genome_organisation_metadata$Genome.Size[genome_organisation_metadata$fasta == genome_summary$genome_name[i]][1]
}
genome_summary$genome_size[genome_summary$genome_name == "rosCan"] = genome_organisation_metadata$Genome.Size[genome_organisation_metadata$fasta == "rosCan_S27_v1.fasta"][1]

pdf(file = "/home/pwlodzimierz/ToL/upload_files/38_similarity_between_within_summary_scatters/38_scatter_relative_similarity_vs_genome_size.pdf")
plot(x = genome_summary$genome_size, 
     y = genome_summary$diff_within_minus_between,
     xlab = "genome size, bp",
     ylab = "Similarity within minus similarity between chromosomes",
     col = point_colors,
     pch = 16)
text(x = min(genome_summary$genome_size) + max(genome_summary$genome_size)/100, y = 30, 
     paste0("p=",round(cor.test(genome_summary$genome_size, genome_summary$diff_within_minus_between)$p.value,4)), adj = 0)
dev.off()





### mess around






pdf(file = "/home/pwlodzimierz/ToL/upload_files/38_similarity_between_within_summary_scatters/38_scatter_simialrity_within_vs_between_chromosomes_alt.pdf")
# Create a point symbol vector based on p-value thresholds
point_pch <- ifelse(genome_summary$p.val < 0.001, 15,
                    ifelse(genome_summary$p.val < 0.01, 17,
                           ifelse(genome_summary$p.val < 0.05, 16, 1)))

# Create a color vector based on group
point_colors <- ifelse(genome_summary$group == "Plant", "#8ac926",
                       ifelse(genome_summary$group == "Invertebrate", "#3f37c9", "#f72585"))

# Plot with p-value-based symbols and group-based colors
plot(x = genome_summary$mean_sim_within,
     y = genome_summary$mean_sim_between,
     pch = point_pch,
     col = point_colors,
     xlab = "Mean Similarity Within Chromosomes",
     ylab = "Mean Similarity Between Chromosomes",
     main = "Genome Similarity Scatter Plot by Group",
     cex = 1.2,
     xlim = c(40,100), ylim = c(40,100))

# Add a legend for group colors
legend("topleft", 
       legend = c("Plant", "Invertebrate", "Chordate"),
       col = c("#8ac926", "#3f37c9", "#f72585"),
       pch = 16,
       title = "Group",
       cex = 0.8)

# Add a legend for p-value symbols
legend("bottomright", 
       legend = c("p ≥ 0.05", "p < 0.05", "p < 0.01", "p < 0.001"),
       pch = c(1, 16, 17, 15),
       col = "black",
       title = "P-value",
       cex = 0.8)

lines(x = c(40,100), y = c(40,100), lty = "aa")

for(i in 1 : nrow(genome_summary)) {
  text(x = genome_summary$mean_sim_within[i] + 1,
       y = genome_summary$mean_sim_between[i] - 0, 
       genome_summary$genome_name[i],
       col =  ifelse(genome_summary$group[i] == "Plant", "#8ac926",
                     ifelse(genome_summary$group[i] == "Invertebrate", "#3f37c9", "#f72585")),
       cex = 0.15,
       adj = 0)
  
  
}
dev.off()




