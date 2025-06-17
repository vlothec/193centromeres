grand_tables <- list.files(path = "/home/pwlodzimierz/ToL/upload_files/grand_tables", full.names = T)
grand_tables <- grand_tables[!grepl("grandest_table.csv", grand_tables)]

genome_sizes <- NULL
centromeric_repeats_total_bp <- NULL

for(i in seq_along(grand_tables)) {
  cat(i, "\n")
  table <- read.csv(file = grand_tables[i])
  
  genome_sizes <- c(genome_sizes, table$size_chromosomes.bp)
  centromeric_repeats_total_bp <- c(centromeric_repeats_total_bp, table$repeats_centromeric.bp)
  
  
}

table <- read.csv("/home/pwlodzimierz/ToL/upload_files/grand_tables/grandest_table.csv")

setwd("/home/pwlodzimierz/ToL/upload_files")

pdf("grand_table_scatter_genome_bp_vs_repeats_bp.pdf")
plot(table$size_chromosomes.bp[table$repeats_centromeric.bp > 0], table$repeats_centromeric.bp[table$repeats_centromeric.bp > 0])
dev.off()