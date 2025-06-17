#!/usr/bin/env Rscript
.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))

taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 193
print(i)


library(Biostrings)
library(GenomicRanges)
library(Rsamtools)
library(ggplot2)


### load data
data_directories <- list.dirs(path = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs", recursive = FALSE, full.names = TRUE)

assembly_name_long <- strsplit(data_directories[i], split = "/")[[1]][7]
assembly_name_short <- strsplit(assembly_name_long, split = ".fa")[[1]][1]

assembly_files <- list.files(path = "/home/pwlodzimierz/ToL/Assemblies/fastas_2021_Michael", recursive = FALSE, full.names = TRUE)
assembly_files <- assembly_files[-grep("fa.fai", assembly_files)]

assembly_file <- assembly_files[grep(assembly_name_short, assembly_files)]

setwd(data_directories[i])

repeats_file <- list.files(path = ".", pattern = "repeats_filtered.csv", full.names = TRUE)
if(length(repeats_file) != 1) {cat("issue with repeats files: \n", repeats_file); stop()}

edta_file <- list.files(path = ".", pattern = "_edta_filtered.csv.reassigned", full.names = TRUE)
if(length(edta_file) != 1) {cat("issue with edta files: \n", edta_file); stop()}


repeats <- read.csv(repeats_file); cat("repeats in\n")
edta <- read.csv(edta_file, header = F); cat("edta in\n")
indexFa(assembly_file)
fasta_file <- FaFile(assembly_file)
fasta_sequences <- readDNAStringSet(assembly_file); cat("fasta in\n")  # Your FASTA file

satellites_metadata <- read.csv(file = "/home/pwlodzimierz/ToL/curated_satellites_repDec24_jan2025.csv")
chr_sizes_metadata <- read.csv(file = "/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full.Ian.csv")

satellites_metadata <- satellites_metadata[satellites_metadata$Genome == assembly_name_long, ]


chr_sizes_metadata <- chr_sizes_metadata[chr_sizes_metadata$assembly.name == assembly_name_long, ]
chr_sizes_metadata <- chr_sizes_metadata[chr_sizes_metadata$is.chr == 1, ]

TRASH_names <- NULL
if(nrow(satellites_metadata) != 0) {
  for(j in 1 : nrow(satellites_metadata)) {
    TRASH_names <- c(TRASH_names, strsplit(satellites_metadata$TRASH_name_dec2024runs[j], split = ";")[[1]])
  }
}
repeats_full <- repeats

### comment this to make centromeric repeat only files
# repeats <- repeats[repeats$new_class %in% TRASH_names, ]


### sample repeats, use only 50 repeats from each class that has more than 50 repeats in it
table_classes <- sort(table(repeats$new_class), decreasing = T)
table_classes <- table_classes[table_classes >= 50]

repeats_new <- NULL
for(j in seq_along(table_classes)) {
  repeats_new <- rbind(repeats_new, repeats[repeats$new_class == names(table_classes)[j],][sample(1 : nrow(repeats[repeats$new_class == names(table_classes)[j], ]), 50), ])
}

repeats <- repeats_new

### adjust edta start end end by +1

edta$V5[edta$V5 == 0] = 1

edta$width <- edta$V6 - edta$V5
edta <- edta[edta$width >= 500, ]
edta <- edta[edta$V4 != "repeat_region", ]

# Load FASTA and extract EDTA sequences
edta_gr <- GRanges(seqnames = edta$V2, ranges = IRanges(start = edta$V5, end = edta$V6))
edta_sequences <- getSeq(fasta_file, edta_gr)

# Create FASTA headers with seqID:start-end
edta_headers <- paste0(edta$V2, ":", edta$V5, "-", edta$V6)
edta_fasta <- DNAStringSet(edta_sequences)
names(edta_fasta) <- edta_headers

# Write to file
writeXStringSet(edta_fasta, "edta_db.fasta")

# Create BLAST Database Use makeblastdb to build the database:
system("makeblastdb -in edta_db.fasta -dbtype nucl -out edta_blastdb")

# Prepare Repeat Sequences as Query Export repeat sequences to a FASTA file:
repeat_fasta <- DNAStringSet(repeats$sequence)
names(repeat_fasta) <- paste0("repeat_", seq_along(repeats$sequence), "_", repeats$seqID, ":", repeats$start, "-", repeats$end)
writeXStringSet(repeat_fasta, "repeats_query.fasta")

# Run BLAST Search Use blastn to search repeats against the "edta" database:
system("blastn -query repeats_query.fasta -db edta_blastdb -out blast_results_allreps.txt -outfmt 6 -evalue 1e-5")



# Parse BLAST Results in R Load and process the results:
blast_results <- read.table("blast_results_allreps.txt", header = FALSE, sep = "\t",
                            col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", 
                                          "gapopen", "qstart", "qend", "sstart", "send", 
                                          "evalue", "bitscore"))

# Extract repeat and edta metadata
blast_results$repeat_idx <- as.numeric(gsub("repeat_([0-9]+)_.*", "\\1", blast_results$qseqid))
blast_results$repeat_seqID <- gsub("repeat_[0-9]+_([^:]+):.*", "\\1", blast_results$qseqid)
blast_results$repeat_start <- as.numeric(gsub("repeat_[0-9]+_[^:]+:([0-9]+)-.*", "\\1", blast_results$qseqid))
blast_results$repeat_end <- as.numeric(gsub("repeat_[0-9]+_[^:]+:[0-9]+-([0-9]+)", "\\1", blast_results$qseqid))

blast_results$edta_seqID <- gsub(":.*", "", blast_results$sseqid)
blast_results$edta_start <- as.numeric(gsub(".*:([0-9]+)-.*", "\\1", blast_results$sseqid))
blast_results$edta_end <- as.numeric(gsub(".*-([0-9]+)$", "\\1", blast_results$sseqid))
blast_results$edta_width <- blast_results$edta_end - blast_results$edta_start + 1


edta <- read.csv(edta_file, header = F); cat("edta in\n")
edta$V5[edta$V5 == 0] = 1

# Rename edta columns for clarity (if not already done)
colnames(edta)[c(2, 5, 6, 4)] <- c("seqID", "start", "end", "edta_id")

# Merge blast_results with edta to get edta_id
blast_results <- merge(blast_results, 
                       edta[, c("seqID", "start", "end", "edta_id")],
                       by.x = c("edta_seqID", "edta_start", "edta_end"),
                       by.y = c("seqID", "start", "end"),
                       all.x = TRUE)

# Merge with original repeat data
repeats$idx <- seq_len(nrow(repeats))
match_df <- merge(blast_results, 
                  repeats[, c("idx", "sequence", "seqID", "start", "end", "new_class")], 
                  by.x = "repeat_idx", by.y = "idx", all.x = TRUE)



# filter out hits that are annotated by TRASH

# Create GRanges for hits and repeats
hit_gr <- GRanges(match_df$edta_seqID, 
                  IRanges(pmin(match_df$edta_start + match_df$sstart - 1, match_df$edta_start + match_df$send - 1),
                          pmax(match_df$edta_start + match_df$sstart - 1, match_df$edta_start + match_df$send - 1)),
                  hit_idx = seq_len(nrow(match_df)))
repeats_gr <- GRanges(repeats$seqID, IRanges(repeats$start, repeats$end))

# Calculate total overlap percentage per hit
overlaps <- findOverlaps(hit_gr, repeats_gr)
overlap_width <- width(pintersect(hit_gr[queryHits(overlaps)], repeats_gr[subjectHits(overlaps)]))
overlap_pct <- overlap_width / width(hit_gr)[queryHits(overlaps)] * 100
total_overlap_pct <- tapply(overlap_pct, queryHits(overlaps), sum, default = 0)

# Assign overlap percentages to all hits (0 for non-overlapping)
match_df$overlap_pct <- 0
match_df$overlap_pct[as.integer(names(total_overlap_pct))] <- total_overlap_pct

# Filter hits and report
filtered_match_df <- match_df[match_df$overlap_pct <= 50, ]
cat("Discarded", nrow(match_df) - nrow(filtered_match_df), "hits with >50% overlap with repeats\n")


# Apply significance filter
significant_hits <- filtered_match_df[filtered_match_df$pident >= 90 & filtered_match_df$evalue <= 1e-5, ]
print(head(significant_hits))


write.csv(significant_hits, paste0(assembly_name_short, "_blast_hits_allreps.csv"), row.names = FALSE)

match_df <- significant_hits


gff_df <- data.frame(
  seqid = match_df$edta_seqID,
  source = "BLAST",
  type = "repeat_region",
  start = ifelse(match_df$sstart < match_df$send, 
                 match_df$edta_start + match_df$sstart - 1, 
                 match_df$edta_start + match_df$send - 1),
  end = ifelse(match_df$sstart < match_df$send, 
               match_df$edta_start + match_df$send - 1, 
               match_df$edta_start + match_df$sstart - 1),
  score = match_df$bitscore,
  strand = ifelse(match_df$sstart < match_df$send, "+", "-"),
  phase = ".",
  attributes = paste0(
    "ID=repeat_hit_", match_df$repeat_idx, "_", seq_len(nrow(match_df)), ";",
    "repeat_idx=", match_df$repeat_idx, ";",
    "new_class=", match_df$new_class, ";",
    "edta_id=", match_df$edta_id, ";",
    "pident=", match_df$pident, ";",
    "evalue=", match_df$evalue
  )
)

# Sort and write GFF3
gff_df <- gff_df[order(gff_df$seqid, gff_df$start), ]
gff_header <- "##gff-version 3\n"
writeLines(gff_header, paste0(assembly_name_short, "_repeat_hits_allreps.gff3"))
write.table(gff_df, paste0(assembly_name_short, "_repeat_hits_allreps.gff3"), append = TRUE, sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

writeLines(gff_header, paste0("/home/pwlodzimierz/ToL/upload_files/23_gff_repeat_blast_against_tes/", assembly_name_short, "_repeat_hits_allreps.gff3"))
write.table(gff_df, paste0("/home/pwlodzimierz/ToL/upload_files/23_gff_repeat_blast_against_tes/", assembly_name_short, "_repeat_hits_allreps.gff3"), append = TRUE, sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)


significant_hits$label <- significant_hits$new_class

# Calculate hit coordinates (genome-wide)
significant_hits$hit_start <- ifelse(significant_hits$sstart < significant_hits$send, 
                                     significant_hits$edta_start + significant_hits$sstart - 1, 
                                     significant_hits$edta_start + significant_hits$send - 1)
significant_hits$hit_end <- ifelse(significant_hits$sstart < significant_hits$send, 
                                   significant_hits$edta_start + significant_hits$send - 1, 
                                   significant_hits$edta_start + significant_hits$sstart - 1)

# Unique edta annotations with hits (base R replacement for dplyr)
edta_with_hits <- unique(significant_hits[, c("edta_seqID", "edta_start", "edta_end")])
edta_with_hits <- merge(edta_with_hits, 
                        edta[, c("seqID", "start", "end", "edta_id")], 
                        by.x = c("edta_seqID", "edta_start", "edta_end"), 
                        by.y = c("seqID", "start", "end"), 
                        all.x = TRUE)
edta_with_hits$length <- edta_with_hits$edta_end - edta_with_hits$edta_start + 1
edta_with_hits <- edta_with_hits[order(-edta_with_hits$length), ]  # Sort by length, longest to shortest

# Function to prepare data for one edta region 
prepare_region_data <- function(seq_id, edta_start, edta_end, edta_id) {
  buffer <- 10000
  midpoint <- (edta_start + edta_end) / 2
  region_start <- midpoint - buffer
  region_end <- midpoint + buffer
  region_start <- max(1, region_start)
  
  edta_subset <- edta[edta$seqID == seq_id & 
                        edta$start <= region_end & 
                        edta$end >= region_start, ]
  
  edta_subset$type <- "edta"
  edta_subset$label <- edta_subset$edta_id
  n_edta <- nrow(edta_subset)
  edta_subset$y_pos <- seq(2, 2.5, length.out = n_edta)
  
  repeats_subset <- repeats_full[repeats_full$seqID == seq_id & 
                                   repeats_full$start <= region_end & 
                                   repeats_full$end >= region_start, ]
  
  if (nrow(repeats_subset) > 0) {
    repeats_subset$type <- "repeat"
    repeats_subset$label <- repeats_subset$new_class
    repeats_subset$y_pos <- 1
  } else {
    repeats_subset <- data.frame()
  }
  
  hits_subset <- significant_hits[significant_hits$edta_seqID == seq_id & 
                                    significant_hits$edta_start == edta_start & 
                                    significant_hits$edta_end == edta_end, ]
  
  hits_subset_selected <- data.frame()
  if (nrow(hits_subset) > 0) {
    hits_subset$type <- "hit"
    n_hits <- nrow(hits_subset)
    hits_subset$y_pos <- seq(0, 0.5, length.out = n_hits)
    hits_subset_selected <- hits_subset[, c("edta_seqID", "hit_start", "hit_end", "type", "label", "y_pos")]
    colnames(hits_subset_selected) <- c("seqID", "start", "end", "type", "label", "y_pos")
  }
  
  all_annotations <- data.frame()
  if (nrow(edta_subset) > 0) {
    all_annotations <- rbind(all_annotations, 
                             edta_subset[, c("seqID", "start", "end", "type", "label", "y_pos")])
  }
  if (nrow(repeats_subset) > 0) {
    all_annotations <- rbind(all_annotations, 
                             repeats_subset[, c("seqID", "start", "end", "type", "label", "y_pos")])
  }
  if (nrow(hits_subset_selected) > 0) {
    all_annotations <- rbind(all_annotations, hits_subset_selected)
  }
  
  list(region_start = region_start, region_end = region_end, annotations = all_annotations)
}

# Generate plot data
plot_data <- lapply(seq_len(nrow(edta_with_hits)), function(i) {
  row <- edta_with_hits[i, ]
  prepare_region_data(row$edta_seqID, row$edta_start, row$edta_end, row$edta_id)
})

plots_per_page <- 10

# Open PDF device for multi-page output
pdf(paste0(assembly_name_short, "_allreps_repeats_vs_edta_hits.pdf"), width = 15, height = plots_per_page*1,1)

# Loop through each page
for (page in 1 : ceiling(length(plot_data) / plots_per_page)) {
  cat(page, "/", ceiling(length(plot_data) / plots_per_page), "\n")
  
  par(mfrow = c(plots_per_page, 1), mar = c(2,0,0,0))
  for(j in 1 : plots_per_page) {
    if(((page-1) * plots_per_page + j) > length(plot_data)) next
    
    panel_data <- plot_data[[(page-1) * plots_per_page + j]][[3]]
    
    plot(x = NULL, y = NULL, xlim = c(min(panel_data$start), max(panel_data$end)), ylim = c(0,7))
    
    panel_edta <- panel_data[panel_data$type == "edta",]
    panel_repeat <- panel_data[panel_data$type == "repeat",]
    panel_hit <- panel_data[panel_data$type == "hit",]
    
    for(k in seq_len(nrow(panel_edta))) {
      sample.y = 2.5 + sample(0:400, 1)/100
      lines(x = c(panel_edta$start[k], panel_edta$end[k]), y = rep(sample.y, 2), col = "#00ee0080", lwd = 3)
      text(x = panel_edta$start[k] + (panel_edta$end[k] -  panel_edta$start[k]) / 2, y = sample.y+0.2, panel_edta$label[k])
    }
    for(k in seq_len(nrow(panel_repeat))) {
      lines(x = c(panel_repeat$start[k], panel_repeat$end[k]), y = rep(1.5 + sample(0:100, 1)/100, 2), col = "#0000ee80", lwd = 3)
    }
    for(k in seq_len(nrow(panel_hit))) {
      lines(x = c(panel_hit$start[k], panel_hit$end[k]), y = rep(0.5 + sample(0:100, 1)/100, 2), col = "#ee000080", lwd = 3)
    }
    
    text(x = min(panel_data$start), y = 6.9, panel_data$seqID[1], adj = 0)
    
    cex_text <- 1
    col_text <- "#ee000080"
    if(sum(TRASH_names %in% unique(panel_hit$label)) > 0) {
      cex_text <- 1.5
      col_text <- "#ee00ee"
    }
    
    unique_repeat_labels <- paste(unique(panel_hit$label), collapse = "; ")
    text(x = min(panel_data$start), y = 0.1, unique_repeat_labels, adj = 0, col = col_text, cex = cex_text)
    
  }
 
}

# Close PDF device
dev.off()

file.copy(from = paste0(assembly_name_short, "_allreps_repeats_vs_edta_hits.pdf"), to = paste0("/home/pwlodzimierz/ToL/upload_files/23_gff_repeat_blast_against_tes/", assembly_name_short, "_allreps_repeats_vs_edta_hits.pdf"), overwrite = T)

































