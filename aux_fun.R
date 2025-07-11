

calc.GC.in.a.window = function(start, sequence, end) return(GC(sequence[start:end]))

calculate.GC.in.windows = function(windows.starts, windows.ends, sequence)
{
  if(length(sequence) < windows.ends[length(windows.ends)]) stop("calculate.GC.in.windows: Something's off with he sequence, either to short or not a char vector")
  return(100 * unlist(lapply(seq_along(windows.starts), function(X) calc.GC.in.a.window(windows.starts[X], sequence, windows.ends[X]))))
}

calculate.repeats.sizes.in.windows = function(windows.starts, repeat.starts, repeat.lengths, windows.ends)
{
  score_in_bins <- unlist(lapply(seq_along(windows.starts), function(X) {
    a = median(repeat.lengths[(repeat.starts > windows.starts[X]) & (repeat.starts <= windows.ends[X])])
    if(is.na(a)) return(0)
    return(a)
  }))
  return(score_in_bins)
}

calc_closest_dist_in_a_win = function(sequence, kmer = 10, min_cov = 0.2) {
  kmers = unlist(lapply(seq_len(length(sequence) - kmer - 1), function(X) paste0(sequence[X : (X + kmer + 1)], collapse = "")[[1]]))
  a = table(kmers)
  b = names(a)[which(a > 1)]
  if(length(b) == 0) return(0)
  if((sum(a[which(a > 1)]) / length(kmers)) < min_cov) return(0)
  distances = NULL
  for(i in seq_along(b)) {
    positions = which(kmers == b[1])
    distances = c(distances, positions[2:length(positions)] - positions[1: (length(positions) - 1)])
  }
  distances = table(distances)
  distances = as.numeric(names(distances[order(distances, decreasing = TRUE)][1]))
  return(distances)
}

calculate_closest_dist_in_win = function(windows.starts, windows.ends, sequence, max_repeat) {
  adjusted_window_starts = windows.starts - max_repeat / 2
  adjusted_windows.ends = windows.ends + max_repeat / 2
  adjusted_window_starts[adjusted_window_starts < 1] <- 1
  adjusted_windows.ends[adjusted_windows.ends > length(sequence)] <- length(sequence)
  return(unlist(lapply(seq_along(adjusted_window_starts), function(X) calc_closest_dist_in_a_win(sequence[adjusted_window_starts[X] : adjusted_windows.ends[X]]))))
}

calculate.repeats.percentage.in.windows = function(windows.starts, repeat.starts, repeat.lengths, sequence.length = 0)
{
  # if(sequence.length < windows.starts[length(windows.starts)]) print("calculate.repeats.percentage.in.windows: sequence length is shorter than last window")
  # if(sequence.length < max(repeat.starts)) print("calculate.repeats.percentage.in.windows: sequence length is shorter than last repeat start")
  
  if(length(windows.starts) == 1) return(sum(repeat.lengths) / (sequence.length - windows.starts[1]))
  
  bins.breaks = c(windows.starts, sequence.length)
  
  repeat.lengths = repeat.lengths[repeat.starts > min(windows.starts)]
  repeat.starts = repeat.starts[repeat.starts > min(windows.starts)]
  
  repeat.lengths = repeat.lengths[repeat.starts < sequence.length]
  repeat.starts = repeat.starts[repeat.starts < sequence.length]
  
  hist.data = hist(rep(repeat.starts, repeat.lengths), breaks = bins.breaks, plot = F)
  
  return(100 * hist.data$counts / (bins.breaks[2:length(bins.breaks)] - bins.breaks[1:(length(bins.breaks) - 1)]))
}

genomic_bins_starts <- function(start = 1, end = 0, bin_number = 0, bin_size = 0) {
  if (bin_number > 0 && bin_size > 0) stop("genomic bins starts: Use either bin number or bin size")
  if (bin_number == 0 && bin_size == 0) stop("genomic bins starts: Use either bin number or bin size")
  if (end < start) stop("genomic bins starts: End smaller than start will not work too well...")
  if (bin_size >= (end - start)) return(start)
  
  if (bin_number > 0) {
    seq_per_bin <- (end - start + 1) %/% bin_number
    remaining_seq <- (end - start + 1) %% bin_number
    bin_sizes <- rep(seq_per_bin, bin_number)
    if (remaining_seq > 0) {
      add_remaining_here <- sample(1:bin_number, remaining_seq)
      bin_sizes[add_remaining_here] <- bin_sizes[add_remaining_here] + 1
    }
    start_positions <- bin_sizes
    for (i in 2 : length(start_positions)) {
      start_positions[i] <- start_positions[i] + start_positions[i - 1]
    }
    start_positions <- start_positions - start_positions[1] + start
    remove(seq_per_bin, remaining_seq)
    
    return(start_positions)
  }
  if (bin_size > 0) {
    if ((end - start) < bin_size) return(start)
    start_positions <- seq(start, (end - bin_size), bin_size)
    if ((end - start_positions[length(start_positions)]) < (bin_size / 2)) { #if future last win length is less than half a bin size
      start_positions <- start_positions[-length(start_positions)] # remove the last one
    }
    return(start_positions)
  }
  return(NA)
}

kmer_compare = function(sequence_against, sequence_to_compare, kmer = 8) {
  sequence_against_ext = paste0(sequence_against, sequence_against, collapse = "")
  sequence_to_compare_ext = paste0(sequence_to_compare, sequence_to_compare, collapse = "")
  sequence_to_compare_ext_rv = revCompString(sequence_to_compare_ext)
  sequence_against_kmers = lapply(seq_len(nchar(sequence_against)), function(X) substr(sequence_against_ext, X, (X + kmer - 1)))
  sequence_to_compare_kmers = lapply(seq_len(nchar(sequence_to_compare)), function(X) substr(sequence_to_compare_ext, X, (X + kmer - 1)))
  sequence_to_compare_kmers_rv = lapply(seq_len(nchar(sequence_to_compare)), function(X) substr(sequence_to_compare_ext_rv, X, (X + kmer - 1)))
  score_fw = sum(sequence_against_kmers %in% sequence_to_compare_kmers) / length(sequence_against_kmers)
  score_rv = sum(sequence_against_kmers %in% sequence_to_compare_kmers_rv) / length(sequence_against_kmers)
  if(score_rv > score_fw) score_fw = score_rv
  return(score_fw)
}

revCompString = function(DNAstr) 
{
  return(tolower(toString(Biostrings::reverseComplement(Biostrings::DNAString(DNAstr)))))
}

read.fasta.and.list = function(file = "")
{
  if(!file.exists(file)) {warning(paste0("File ", file, " does not exist")); return()}
  fasta_full = read.dna(file = file, format = "fasta", as.character = T, as.matrix = F) # this ape function is faster than read.fasta()
  if(typeof(fasta_full) == "character") 
  {
    names.save = rownames(fasta_full)
    fasta_full = list(fasta_full)
    names(fasta_full) = names.save
    remove(names.save)
  }
  names(fasta_full) = lapply(names(fasta_full), function(X) strsplit(X, split = " ")[[1]][1])
  return(fasta_full)
}

export_gff = function(annotations.data.frame = "", output = ".", file.name = "gff.table",
                      seqid = ".", source = ".", type = ".", start = "0", end = "0", score = ".", strand = ".", phase = ".", attributes = ".", 
                      attribute.names = ".") 
{
  # print("Export gff function")
  # print("GFF function, either add numeric value (or values) for which column(s) contains data, or text string which will be universal for that gff field")
  # print("attribute.names are names added to each attribute as in attributes columns selection")
  
  if(!is.data.frame(annotations.data.frame)) stop("Provide a data frame with annotations")
  
  if(nrow(annotations.data.frame) == 0) stop("The data frame provided is does contain at least one row")
  
  output = sub("[/\\\\]$", "", output)
  
  if(!is.character(seqid)) 
    if(length(seqid) == 1) seqid = annotations.data.frame[, seqid]
  else seqid = do.call(paste, c(annotations.data.frame[, seqid], sep = "_"))
  else seqid = rep(seqid, nrow(annotations.data.frame))
  
  if(!is.character(source)) 
    if(length(source) == 1) source = annotations.data.frame[, source]
  else source = do.call(paste, c(annotations.data.frame[, source], sep = "_"))
  else source = rep(source, nrow(annotations.data.frame))
  
  if(!is.character(type)) 
    if(length(type) == 1) type = annotations.data.frame[, type]
  else type = do.call(paste, c(annotations.data.frame[, type], sep = "_"))
  else type = rep(type, nrow(annotations.data.frame))
  
  if(!is.character(start)) 
    if(length(start) == 1) start = annotations.data.frame[, start]
  else start = do.call(sum, c(annotations.data.frame[, start]))
  else start = rep(start, nrow(annotations.data.frame))
  
  if(!is.character(end)) 
    if(length(end) == 1) end = annotations.data.frame[, end]
  else end = do.call(sum, c(annotations.data.frame[, end]))
  else end = rep(end, nrow(annotations.data.frame))
  
  if(!is.character(score)) 
    if(length(score) == 1) score = annotations.data.frame[, score]
  else score = do.call(sum, c(annotations.data.frame[, score]))
  else score = rep(score, nrow(annotations.data.frame))
  
  if(!is.character(strand)) 
    if(length(strand) == 1) strand = annotations.data.frame[, strand]
  else strand = do.call(paste, c(annotations.data.frame[, strand], sep = "_"))
  else strand = rep(strand, nrow(annotations.data.frame))
  
  if(!is.character(phase)) 
    if(length(phase) == 1) phase = annotations.data.frame[, phase]
  else phase = do.call(paste, c(annotations.data.frame[, phase], sep = "_"))
  else phase = rep(phase, nrow(annotations.data.frame))
  
  if(!is.character(attributes)) 
    if(length(attributes) == 1) attributes = paste(attribute.names, annotations.data.frame[, attributes], sep = "")
  else attributes = apply(annotations.data.frame[, attributes], 1, function(x) {    paste0(attribute.names, x, collapse = ";")  })
  else attributes = rep(attributes, nrow(annotations.data.frame))
  
  
  gff_format = data.frame(seqid, # sequence ID, like chromosme name
                          source,# annotation source, like TRASH
                          type,  # annotation type, like "gene", "satellite_DNA"
                          start,
                          end,
                          score,
                          strand,# +, - or "." I think
                          phase, # 1 2 or 3 I think
                          attributes) # other things
  
  
  options(scipen=10)
  write.table(x = gff_format, file = paste(output, "/", file.name, ".gff", sep = ""), quote = FALSE, sep = "\t", eol = "\r", row.names = FALSE, col.names = FALSE)
  options(scipen=0)
  
}

consensus_N = function(alignment, N)
{
  alignment.matrix = tolower(as.matrix(Biostrings::unmasked(alignment)))
  # alignment.matrix = alignment
  frequencies = vector(mode = "numeric", length = ncol(alignment.matrix))
  
  for(i in seq_along(frequencies))
  {
    frequencies[i] = frequencies[i] + sum(1 * (alignment.matrix[,i] != "-"))
  }
  
  is.in.consensus = vector(mode = "logical", length = ncol(alignment.matrix))
  
  for(i in 1 : N)
  {
    is.in.consensus[order(frequencies, decreasing = TRUE)[i]] = TRUE
  }
  
  consensus = vector(mode = "character", length = N)
  consensus.ID = 1
  for(i in seq_len(ncol(alignment.matrix)))
  {
    if(is.in.consensus[i])
    {
      consensus[consensus.ID] = c("g","c","t","a")[which.max(c(length(which(alignment.matrix[,i] == "g")),
                                                               length(which(alignment.matrix[,i] == "c")),
                                                               length(which(alignment.matrix[,i] == "t")),
                                                               length(which(alignment.matrix[,i] == "a"))))]
      consensus.ID = consensus.ID + 1
    }
  }
  
  return(tolower(paste(consensus, collapse = "")))
}

