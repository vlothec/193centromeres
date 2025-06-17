library(seqinr)

setwd("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs/drGeuUrba1.1.fa")

# Read input files
# edta <- read.csv("./drGeuUrba1.1_edta_filtered.csv.reassigned", header = FALSE)
edta <- read.table("/home/pwlodzimierz/ToL/geum_stitching/drGeuUrba1.1.fa.mod.EDTA.intact.gff3.clean", header = FALSE, sep = "\t")
fasta <- read.fasta("/home/pwlodzimierz/ToL/Assemblies/fastas_2021_Michael/drGeuUrba1.1.fa")

edta$V15 = ""
for(i in 1 : nrow(edta)) {
  edta$V15[i] = strsplit(strsplit(x = edta$V9[i], split = ";Method=")[[1]][2], split = ";")[[1]][1]
}


# Step 1: Filter for Gypsy_LTR_retrotransposon
# gypsy <- edta[edta$V4 == "Gypsy_LTR_retrotransposon", ]

gypsy <- edta

# gypsy_filtered <- gypsy[gypsy$V6 - gypsy$V5 > 9000,] # this just assumes that long elements are intact
gypsy_filtered <- gypsy[gypsy$V15== "structural",] # this means the element is intact


# in case annotations are overlapping, here is a mergiong funtion:
merge_regions <- function(regions) {
  if (nrow(regions) <= 1) return(regions)
  # Sort by start position
  regions <- regions[order(regions$V4), ]
  merged <- data.frame(V4 = regions$V4[1], V5 = regions$V5[1])
  for (i in 2:nrow(regions)) {
    if(i %in% round((1:100) * nrow(regions)/100)) cat("#")
    if (regions$V4[i] <= merged$V5[nrow(merged)] + 1) {
      # Overlapping or adjacent: extend the end position
      merged$V5[nrow(merged)] <- max(merged$V5[nrow(merged)], regions$V5[i])
    } else {
      # Non-overlapping: add new region
      merged <- rbind(merged, data.frame(V4 = regions$V4[i], V5 = regions$V5[i]))
    }
  }; cat("\n")
  return(merged)
}

# Create a data frame to track removed regions
removed_regions <- data.frame(
  seqid = character(),
  original_start = integer(),
  original_end = integer(),
  removed_length = integer(),
  new_start = integer(),
  new_end = integer(),
  stringsAsFactors = FALSE
)

# Process each sequence in the FASTA file
new_fasta <- list()
for (seq_name in names(fasta)) {
  seq <- fasta[[seq_name]]
  seq_regions <- gypsy_filtered[gypsy_filtered$V1 == seq_name, ]
  
  if (nrow(seq_regions) == 0) {
    new_fasta[[seq_name]] <- seq
    next
  }
  
  # Merge overlapping regions
  seq_regions <- merge_regions(seq_regions)
  
  # Initialize variables for stitching
  new_seq <- c()
  current_pos <- 1
  new_pos <- 1
  
  # Iterate through merged regions to remove
  for (i in 1:nrow(seq_regions)) {
    cat(seq_name, i, nrow(seq_regions), "\n")
    start <- seq_regions$V4[i]
    end <- seq_regions$V5[i]
    
    # Add sequence before the region
    if (current_pos < start) {
      new_seq <- c(new_seq, seq[current_pos:(start-1)])
    }
    
    new_pos <- length(new_seq)
    
    # Record removed region
    removed_length <- end - start + 1
    removed_regions <- rbind(removed_regions, data.frame(
      seqid = seq_name,
      original_start = start,
      original_end = end,
      removed_length = removed_length,
      new_start = new_pos,
      new_end = new_pos + 1,
      stringsAsFactors = FALSE
    ))
    
    current_pos <- end + 1
  }
  
  # Add remaining sequence
  if (current_pos <= length(seq)) {
    new_seq <- c(new_seq, seq[current_pos:length(seq)])
  }
  
  new_fasta[[seq_name]] <- new_seq
}

# Write the removed regions data frame
write.csv(removed_regions, "drGeuUrba1.1_removed_regions_3rd.csv", row.names = FALSE)

# Write the new FASTA file
write.fasta(sequences = new_fasta, names = names(new_fasta), 
            file.out = "drGeuUrba1.1_stitched_3rd.fasta")


############
#

# SECOND ROUND

#
############



library(seqinr)

setwd("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs/drGeuUrba1.1.fa")

# Read input files
# edta <- read.csv("./drGeuUrba1.1_edta_filtered.csv.reassigned", header = FALSE)
edta <- read.table("/home/pwlodzimierz/ToL/geum_stitching/drGeuUrba1.1_stitched_3rd.fasta.F2B.mod.EDTA.intact.gff3", header = FALSE, sep = "\t", skip = 6)
fasta <- read.fasta("./drGeuUrba1.1_stitched_3rd.fasta")

edta$V15 = ""
for(i in 1 : nrow(edta)) {
  edta$V15[i] = strsplit(strsplit(x = edta$V9[i], split = ";Method=")[[1]][2], split = ";")[[1]][1]
}


# Step 1: Filter for Gypsy_LTR_retrotransposon
# gypsy <- edta[edta$V4 == "Gypsy_LTR_retrotransposon", ]
  
  gypsy <- edta

# gypsy_filtered <- gypsy[gypsy$V6 - gypsy$V5 > 9000,] # this just assumes that long elements are intact
gypsy_filtered <- gypsy[gypsy$V15== "structural",] # this means the element is intact


# in case annotations are overlapping, here is a mergiong funtion:
merge_regions <- function(regions) {
  if (nrow(regions) <= 1) return(regions)
  # Sort by start position
  regions <- regions[order(regions$V4), ]
  merged <- data.frame(V4 = regions$V4[1], V5 = regions$V5[1])
  for (i in 2:nrow(regions)) {
    if(i %in% round((1:100) * nrow(regions)/100)) cat("#")
    if (regions$V4[i] <= merged$V5[nrow(merged)] + 1) {
      # Overlapping or adjacent: extend the end position
      merged$V5[nrow(merged)] <- max(merged$V5[nrow(merged)], regions$V5[i])
    } else {
      # Non-overlapping: add new region
      merged <- rbind(merged, data.frame(V4 = regions$V4[i], V5 = regions$V5[i]))
    }
  }; cat("\n")
  return(merged)
}

# Create a data frame to track removed regions
removed_regions <- data.frame(
  seqid = character(),
  original_start = integer(),
  original_end = integer(),
  removed_length = integer(),
  new_start = integer(),
  new_end = integer(),
  stringsAsFactors = FALSE
)

# Process each sequence in the FASTA file
new_fasta <- list()
for (seq_name in names(fasta)) {
  seq <- fasta[[seq_name]]
  seq_regions <- gypsy_filtered[gypsy_filtered$V1 == seq_name, ]
  
  if (nrow(seq_regions) == 0) {
    new_fasta[[seq_name]] <- seq
    next
  }
  
  # Merge overlapping regions
  seq_regions <- merge_regions(seq_regions)
  
  # Initialize variables for stitching
  new_seq <- c()
  current_pos <- 1
  new_pos <- 1
  
  # Iterate through merged regions to remove
  for (i in 1:nrow(seq_regions)) {
    cat(seq_name, i, nrow(seq_regions), "\n")
    start <- seq_regions$V4[i]
    end <- seq_regions$V5[i]
    
    # Add sequence before the region
    if (current_pos < start) {
      new_seq <- c(new_seq, seq[current_pos:(start-1)])
    }
    
    new_pos <- length(new_seq)
    
    # Record removed region
    removed_length <- end - start + 1
    removed_regions <- rbind(removed_regions, data.frame(
      seqid = seq_name,
      original_start = start,
      original_end = end,
      removed_length = removed_length,
      new_start = new_pos,
      new_end = new_pos + 1,
      stringsAsFactors = FALSE
    ))
    
    current_pos <- end + 1
  }
  
  # Add remaining sequence
  if (current_pos <= length(seq)) {
    new_seq <- c(new_seq, seq[current_pos:length(seq)])
  }
  
  new_fasta[[seq_name]] <- new_seq
}

# Write the removed regions data frame
write.csv(removed_regions, "drGeuUrba1.1_removed_regions_5th.csv", row.names = FALSE)

# Write the new FASTA file
write.fasta(sequences = new_fasta, names = names(new_fasta), 
            file.out = "drGeuUrba1.1_stitched_5th.fasta")
















