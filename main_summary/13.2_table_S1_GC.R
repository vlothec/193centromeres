

taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 193
print(i)



.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))
suppressMessages(library(msa))
suppressMessages(library(seqinr)) 
suppressMessages(library(GenomicRanges))
suppressMessages(library(ggplot2))

ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}

do_GC <- T
redo <- F

setwd("/home/pwlodzimierz/ToL/git_ToL")
source("./aux_fun.R")

data_directories <- list.dirs(path = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs", recursive = FALSE, full.names = TRUE)
data_directories <- data_directories[!grepl(pattern = "templated_", data_directories)]
data_directories <- data_directories[grepl(pattern = ".fa", data_directories)]

if(i > length(data_directories)) stop()
if(grepl("Rosa_agrestis_DTOL_chrs_5n.fasta", data_directories[i])) stop()
if(grepl("Rosa_canina_DTOL_chrs_5n.fasta", data_directories[i])) stop()
if(grepl("drGeuRiva1.hap1.1.fa", data_directories[i])) stop()
if(grepl("drGeuRiva1.hap2.1.fa", data_directories[i])) stop()
if(grepl("drRosRugo1.1.fa", data_directories[i])) stop()
if(grepl("drRosSpin1.hap1.1.fa", data_directories[i])) stop()
if(grepl("GCA_000001405.29_GRCh38.p14_genomic.fa", data_directories[i])) stop()

assembly_files <- list.files(path = "/home/pwlodzimierz/ToL/Assemblies/fastas_2021_Michael", recursive = FALSE, full.names = TRUE)
assembly_files <- assembly_files[!grepl(".fai", assembly_files)]


print(paste0(i, " / ", length(data_directories)))
### Load data ================================================================
assembly_name = strsplit(strsplit(data_directories[i], split = ".fa")[[1]][1], split = "v2_out_for_HORs/")[[1]][2]
fasta_name = strsplit(data_directories[i], split = "v2_out_for_HORs/")[[1]][2]
assembly_file = grep(assembly_name, assembly_files)
print(assembly_file)
print(assembly_files[assembly_file])

genomes_organisation_data <- read.csv("/home/pwlodzimierz/ToL/genomes_organisation_type.csv")
genomes_organisation_data <- genomes_organisation_data[genomes_organisation_data$fasta == fasta_name,]

setwd(data_directories[i])
print(getwd())

if(!redo) {
  if(file.exists(paste0("/home/pwlodzimierz/ToL/upload_files/13_GC_grand_tables/", fasta_name, "_", genomes_organisation_data$Centromere.architecture[1],"_GC_grand_table.csv"))) {
    stop("File exists, not redoing")
  }
}

print("reading fasta")

fasta <- read.fasta(file = assembly_files[assembly_file])



print("fasta read")

if(fasta_name %in% c("mRhiSin1.1.fa", "iyVesVulg1.1.fa", "rosCan_S27_v1.fasta")) {
  edta_file = list.files(pattern = paste0(assembly_name, "_edta_filtered.csv"), full.names = TRUE)
  if(length(edta_file) != 1) {warning(paste0(i, " no edta!")); setwd(".."); quit(save = "no", status = 1)}
} else {
  edta_file = list.files(pattern = paste0(assembly_name, "_edta_filtered.csv.reassigned"), full.names = TRUE)
  if(length(edta_file) != 1) {warning(paste0(i, " no edta!")); setwd(".."); quit(save = "no", status = 1)}
  
}

repeat_file = list.files(pattern = "_repeats_filtered.csv", full.names = TRUE)
if(length(repeat_file) != 1) {warning(paste0(i, "No repeats!")); setwd(".."); quit(save = "no", status = 1)}

helixer_file = list.files(pattern = "helixer_filtered.gff", full.names = TRUE)
if(length(helixer_file) != 1) {warning(paste0(i, "No genes!")); setwd(".."); quit(save = "no", status = 1)}

array_file = list.files(pattern = "_arrays_filtered.csv", full.names = TRUE)
if(length(array_file) != 1) {warning(paste0(i, " no arrays!")); setwd(".."); quit(save = "no", status = 1)}

classes_file = list.files(pattern = "_classes_merged_filtered", full.names = TRUE)
if(length(classes_file) != 1) {warning(paste0(i, " no classes!")); setwd(".."); quit(save = "no", status = 1)}

gaps_file = list.files(pattern = paste0(fasta_name, "_centromeric_gaps.csv"), full.names = TRUE)
if(length(gaps_file) != 1) {print(paste0(i, " no gaps!"))} else {
  
  gaps <- read.csv(file = gaps_file)
}

cenarray_file = list.files(pattern = paste0(fasta_name, "_centromeric_arrays.csv"), full.names = TRUE)
if(length(cenarray_file) != 1) {print(paste0(i, " no cen arrays!"))} else {
  
  centromeric_arrays <- read.csv(file = cenarray_file)
}

genomemetadata_file = list.files(pattern = paste0(fasta_name, "_genome_metadata.csv"), full.names = TRUE)
if(length(genomemetadata_file) != 1) {warning(paste0(i, " no genome metadata!")); setwd(".."); quit(save = "no", status = 1)}


print("load annotations")
repeats = read.csv(file = repeat_file)
arrays = read.csv(file = array_file)
classes = read.csv(file = classes_file)
classes$num_ID <- 1 : nrow(classes)
edta = read.csv(file = edta_file, header = FALSE)
genes <- read.table(file = helixer_file, header = FALSE, sep = "\t", skip = 4)
genes_full <- genes
genes <- genes[genes$V3 == "CDS", ]
genome_metadata <- read.csv(file = genomemetadata_file)

if(fasta_name %in% c("mRhiSin1.1.fa", "iyVesVulg1.1.fa")) { 
  repeats$new_class <- repeats$class
  
}

if(i == 178) edta <- edta[-1,]

edta$V5 <- as.numeric(edta$V5)
edta$V6 <- as.numeric(edta$V6)

print("annotations loaded")



which_fasta_is_chr <- which(names(fasta) %in% genome_metadata$chromosome.name[genome_metadata$is.chr == 1])

which_fasta_isnt_chr <- seq_along(fasta)[-which(seq_along(fasta) %in% which_fasta_is_chr)]




satellite_metadata <- read.csv("/home/pwlodzimierz/ToL/curated_satellites_metadata_on_chromosomes_only_may.csv")
centromeric_stallite_for_species <- satellite_metadata[satellite_metadata$Genome == fasta_name,]




chromosomes <- genome_metadata$chromosome.name[genome_metadata$is.chr == 1]
chromosomes_lengths <- genome_metadata$size[genome_metadata$is.chr == 1]
non_chromosomes <- genome_metadata$chromosome.name[genome_metadata$is.chr != 1]
non_chr_lengths <- genome_metadata$size[genome_metadata$is.chr != 1]

# unique(genomes_organisation_data$Centromere.architecture)
# "Holocentric no satellite" "Holocentric satellite"    "Satellite"                "Transposon"               "Mixed"                    "Unknown"           


LTR_names <- c("(?i)^(?!.*non[-]?LTR).*LTR|long_terminal_repeat")
nonLTR_names <-c("line|sine|penelope|nonLTR_retrotransposon")
TIR_names <- c("TIR")
nonTIR_names <- c("helitron|MITE|Maverick|Polinton|Class_II_DNA_Transposon|Tyrosine_Recombinase_Elements|Crypton_Tyrosine_Recombinase")


repeats$is_on_chromosomes <- FALSE
repeats$is_on_chromosomes[repeats$seqID %in% chromosomes] <- TRUE

repeats$is_centromeric <- FALSE
repeats$in_centromere <- FALSE
satellite_names <- NULL
if(nrow(centromeric_stallite_for_species) != 0) {
  for(j in 1 : nrow(centromeric_stallite_for_species)) {
    satellite_names <- strsplit(centromeric_stallite_for_species$TRASH_name_dec2024runs[j], split = ";")[[1]]
    repeats$is_centromeric[repeats$new_class %in% satellite_names] <- TRUE
    
    for(j in 1 : nrow(centromeric_arrays)) {
      if(!centromeric_arrays$in_centromere[j]) next
      repeats$in_centromere[repeats$is_centromeric & repeats$seqID == centromeric_arrays$seqID[j] & (repeats$start %in% (centromeric_arrays$start[j] : centromeric_arrays$end[j]))] <- TRUE
      
    }
    
  }
  
}

### Summary ====================================================================



print("1")
fasta_name <- fasta_name

genus <- genomes_organisation_data$Genus

species <- genomes_organisation_data$Species

if(genomes_organisation_data$Centromere.architecture[1] %in% c("Satellite", "Holocentric satellite") ) {
  centromeric_satellites.1_yes_0_no <- 1
} else {
  centromeric_satellites.1_yes_0_no <- 0
}

if(genomes_organisation_data$Centromere.architecture[1] %in% c("Holocentric no satellite", "Holocentric satellite") ) {
  holocentricity.1_yes_0_no <- 1
} else {
  holocentricity.1_yes_0_no <- 0
}

sequences.no <- nrow(genome_metadata)

chromosomes.no <- sum(genome_metadata$is.chr)

chromosomes_with_centromeric_repeats.no <- 0
if(centromeric_satellites.1_yes_0_no == 1) {
  for(j in seq_along(chromosomes)) {
    if(nrow(repeats[repeats$seqID == chromosomes[j] & repeats$is_centromeric,]) > 0) {
      chromosomes_with_centromeric_repeats.no <- chromosomes_with_centromeric_repeats.no + 1
    }
  }
}

print("2")

size_genome.bp <- sum(unlist(lapply(X = seq_along(fasta), function(X) length(fasta[[X]]))))

size_chromosomes.bp <- sum(chromosomes_lengths)

size_non_chromosomes.bp <- sum(non_chr_lengths)

chromosomes_mean_size.bp <- mean(chromosomes_lengths)

chromosomes_mean_size_SD.bp <- sd(chromosomes_lengths)



size_centromeres.bp <- NA
if(centromeric_satellites.1_yes_0_no == 1) {
  size_centromeres.bp <- sum(centromeric_arrays$width[centromeric_arrays$in_centromere])
}

print("2.1")
size_pericentromeres.bp <- 0
if(centromeric_satellites.1_yes_0_no == 1) {
  for(j in which(genome_metadata$is.chr == 1)) {
    gen_meta_temp <- genome_metadata[j,]
    peri_starts <- as.numeric(strsplit(as.character(gen_meta_temp$pericentromere_start[1]), split = ";")[[1]])
    peri_ends <- as.numeric(strsplit(as.character(gen_meta_temp$pericentromere_end[1]), split = ";")[[1]])
    size_pericentromeres.bp <- size_pericentromeres.bp + sum(width(IRanges(peri_starts, peri_ends)))
  }
  size_pericentromeres.bp <-  size_pericentromeres.bp - size_centromeres.bp
}

centromere_mean_size.bp <- NA
if(centromeric_satellites.1_yes_0_no == 1) {
  cen_sizes_vect <- unlist(lapply(chromosomes, function(X) sum(centromeric_arrays$width[centromeric_arrays$seqID == X])))
  centromere_mean_size.bp <- mean(cen_sizes_vect[cen_sizes_vect != 0])
}

print("2.2")
centromere_mean_size_SD.bp <- NA
if(centromeric_satellites.1_yes_0_no == 1) {
  centromere_mean_size_SD.bp <- sd(cen_sizes_vect[cen_sizes_vect != 0])
}

print("3")
GC_genome.perc <- NA
if(do_GC) {
  GC_genome.perc <- 100 * sum(unlist(lapply(X = seq_along(fasta), function(X) GC(fasta[[X]], NA.GC = T) * length(fasta[[X]]) ))) / size_genome.bp
}

print("3.1")
GC_chromosomes.perc <- NA
if(do_GC) {
  GC_chromosomes.perc <- 100 * sum(unlist(lapply(X = which_fasta_is_chr, function(X) GC(fasta[[X]], NA.GC = T) * length(fasta[[X]]) ))) / size_chromosomes.bp
  
}

print("3.2")
GC_non_chromosomes.perc <- NA
if(do_GC) {
  GC_non_chromosomes.perc <- 100 * sum(unlist(lapply(X = which_fasta_isnt_chr, function(X) GC(fasta[[X]], NA.GC = T) * length(fasta[[X]]) ))) / size_non_chromosomes.bp
  
}

print("3.3")
bp_times_gc <- 0
GC_centromeres.perc <- NA
if(do_GC) {
  if(centromeric_satellites.1_yes_0_no == 1) {
    sum_bp <- 0
    for(j in which(genome_metadata$is.chr == 1)) {
      arrays_temp <- centromeric_arrays[centromeric_arrays$seqID == genome_metadata$chromosome.name[j],]
      if(nrow(arrays_temp) == 0) next
      for(k in seq_len(nrow(arrays_temp))) {
        if(!arrays_temp$in_centromere[k]) next
        bp_times_gc <- bp_times_gc + 100*GC(fasta[[j]][arrays_temp$start[k] : arrays_temp$end[k]], NA.GC = T) * arrays_temp$width[k]
        sum_bp <- sum_bp + arrays_temp$width[k]
      }
    }
    GC_centromeres.perc <- bp_times_gc / sum_bp
  }
}


print("3.4")
bp_times_gc_peri <- 0
sum_bp_peri <- 0
GC_pericentromeres.perc <- NA
if(do_GC) {
  if(centromeric_satellites.1_yes_0_no == 1) {
    for(j in which(genome_metadata$is.chr == 1)) {
      arrays_temp <- centromeric_arrays[centromeric_arrays$seqID == genome_metadata$chromosome.name[j],]
      arrays_temp <- arrays_temp[arrays_temp$in_centromere,]
      if(nrow(arrays_temp) == 0) next
      gen_meta_temp <- genome_metadata[genome_metadata$chromosome.name == genome_metadata$chromosome.name[j],]
      
      peri_starts <- as.numeric(strsplit(as.character(gen_meta_temp$pericentromere_start[1]), split = ";")[[1]])
      peri_ends <- as.numeric(strsplit(as.character(gen_meta_temp$pericentromere_end[1]), split = ";")[[1]])
      if(length(peri_ends) == 1) {
        if(peri_ends == 0) next
      } 
      
      pericentromere_only_gr <- setdiff(IRanges(peri_starts, peri_ends), IRanges(arrays_temp$start, arrays_temp$end))
      
      for(k in seq_along(pericentromere_only_gr)) {
        print(c(start(pericentromere_only_gr)[k], end(pericentromere_only_gr)[k]))
        fasta_extract <- unique(fasta[[j]][start(pericentromere_only_gr)[k] : end(pericentromere_only_gr)[k]])
        if(length(fasta_extract) == 1) {
          if(is.na(fasta_extract)) {
            next
          }
        }
        bp_times_gc_peri <- bp_times_gc_peri + 100*GC(fasta[[j]][start(pericentromere_only_gr)[k] : end(pericentromere_only_gr)[k]], NA.GC = T) * width(pericentromere_only_gr)[k]
        sum_bp_peri <- sum_bp_peri + width(pericentromere_only_gr)[k]
      }
    }
  } else {
    for(j in which(genome_metadata$is.chr == 1)) {
      gen_meta_temp <- genome_metadata[genome_metadata$chromosome.name == genome_metadata$chromosome.name[j],]
      if(gen_meta_temp$pericentromere_start[1] != 0) {
        
        peri_starts <- as.numeric(strsplit(as.character(gen_meta_temp$pericentromere_start[1]), split = ";")[[1]])
        peri_ends <- as.numeric(strsplit(as.character(gen_meta_temp$pericentromere_end[1]), split = ";")[[1]])
        if(length(peri_ends) == 1) {
          if(peri_ends == 0) next
        } 
        pericentromere_only_gr <- IRanges(peri_starts, peri_ends)
        
        for(k in seq_along(pericentromere_only_gr)) {
          fasta_extract <- unique(fasta[[j]][start(pericentromere_only_gr)[k] : end(pericentromere_only_gr)[k]])
          if(length(fasta_extract) == 1) {
            if(is.na(fasta_extract)) {
              next
            }
          }
          bp_times_gc_peri <- bp_times_gc_peri + 100*GC(fasta[[j]][start(pericentromere_only_gr)[k] : end(pericentromere_only_gr)[k]], NA.GC = T) * width(pericentromere_only_gr)[k]
          sum_bp_peri <- sum_bp_peri + width(pericentromere_only_gr)[k]
        }
      }
    }
  }
  GC_pericentromeres.perc <- bp_times_gc_peri / sum_bp_peri
}


print("4")
bp_times_gc_arm <- 0
sum_bp_arm <- 0
GC_arms.perc <- NA
if(do_GC) {
  for(j in which(genome_metadata$is.chr == 1)) {
    gen_meta_temp <- genome_metadata[genome_metadata$chromosome.name == genome_metadata$chromosome.name[j],]
    if(gen_meta_temp$pericentromere_start[1] != 0) {
      
      peri_starts <- as.numeric(strsplit(as.character(gen_meta_temp$pericentromere_start[1]), split = ";")[[1]])
      peri_ends <- as.numeric(strsplit(as.character(gen_meta_temp$pericentromere_end[1]), split = ";")[[1]])
      if(length(peri_ends) == 1) {
        if(peri_ends == 0) next
      } 
      
      dontconsider <- (is.na(peri_starts) | is.na(peri_ends))
      
      peri_starts <- peri_starts[!dontconsider]
      peri_ends <- peri_ends[!dontconsider]
      
      print(c(peri_starts, peri_ends))
      
      if(length(peri_starts) != 0) {
        arm_only_gr <- setdiff(IRanges(1, genome_metadata$size[j]), IRanges(peri_starts, peri_ends))
        for(k in seq_along(arm_only_gr)) {
          fasta_extract <- unique(fasta[[j]][start(arm_only_gr)[k] : end(arm_only_gr)[k]])
          if(length(fasta_extract) == 1) {
            if(is.na(fasta_extract)) {
              next
            }
          }
          bp_times_gc_arm <- bp_times_gc_arm + 100*GC(fasta[[j]][start(arm_only_gr)[k] : end(arm_only_gr)[k]], NA.GC = T) * width(arm_only_gr)[k]
          sum_bp_arm <- sum_bp_arm + width(arm_only_gr)[k]
        }
      }
      
      
    }
    GC_arms.perc <- bp_times_gc_arm / sum_bp_arm
  }
}

cat("a ") 
bp_times_gc <- 0
sum_bp <- 0
gr_of_interest <- 0
GC_exons_total.perc <- NA
if(do_GC) {
  for(j in which(genome_metadata$is.chr == 1)) {
    gen_meta_temp <- genome_metadata[genome_metadata$chromosome.name == genome_metadata$chromosome.name[j],]
    exons_temp <- genes[genes$V1 == genome_metadata$chromosome.name[j],]
    gr_of_interest <- overlapsRanges(IRanges(1, genome_metadata$size[j]), IRanges(exons_temp$V4, exons_temp$V5))
    for(k in seq_along(gr_of_interest)) {
      fasta_extract <- unique(fasta[[j]][start(gr_of_interest)[k] : end(gr_of_interest)[k]])
      if(length(fasta_extract) == 1) {
        if(is.na(fasta_extract)) next
      }
      bp_times_gc <- bp_times_gc + 100*GC(fasta[[j]][start(gr_of_interest)[k] : end(gr_of_interest)[k]], NA.GC = T) * width(gr_of_interest)[k]
      sum_bp <- sum_bp + width(gr_of_interest)[k]
    }
  }
  GC_exons_total.perc <- bp_times_gc / sum_bp
}

remove(bp_times_gc, sum_bp, gr_of_interest)
cat("b ")
bp_times_gc <- 0
sum_bp <- 0
gr_of_interest <- NA
GC_TE_total.perc <- NA
if(do_GC) {
  for(j in which(genome_metadata$is.chr == 1)) {
    gen_meta_temp <- genome_metadata[genome_metadata$chromosome.name == genome_metadata$chromosome.name[j],]
    TE_temp <- edta[edta$V2 == genome_metadata$chromosome.name[j],]
    gr_of_interest <- overlapsRanges(IRanges(1, genome_metadata$size[j]), IRanges(TE_temp$V5, TE_temp$V6))
    for(k in seq_along(gr_of_interest)) {
      fasta_extract <- unique(fasta[[j]][start(gr_of_interest)[k] : end(gr_of_interest)[k]])
      if(length(fasta_extract) == 1) {
        if(is.na(fasta_extract)) next
      }
      bp_times_gc <- bp_times_gc + 100*GC(fasta[[j]][start(gr_of_interest)[k] : end(gr_of_interest)[k]], NA.GC = T) * width(gr_of_interest)[k]
      sum_bp <- sum_bp + width(gr_of_interest)[k]
    }
  }
  GC_TE_total.perc <- bp_times_gc / sum_bp
}
remove(bp_times_gc, sum_bp, gr_of_interest)

GC_repeats_total.perc <- NA
if(do_GC) {
  GC_repeats_total.perc <- 100*GC(strsplit(paste(repeats$sequence, collapse = ""), split = "")[[1]])
}


GC_repeats_centromeric.perc = NA
if(do_GC) {
  if(centromeric_satellites.1_yes_0_no == 1) {
    if(sum(repeats$is_centromeric) != 0) {
      GC_repeats_centromeric.perc <- 100*GC(strsplit(paste(repeats$sequence[repeats$is_centromeric], collapse = ""), split = "")[[1]])
    }
  } 
}


print("6")
GC_repeats_noncentromeric.perc <- NA
if(do_GC) {
  if(centromeric_satellites.1_yes_0_no == 1) {
    if(sum(!repeats$is_centromeric) != 0) {
      GC_repeats_noncentromeric.perc <- 100*GC(strsplit(paste(repeats$sequence[!repeats$is_centromeric], collapse = ""), split = "")[[1]])
    }
    
  }
  
}



grand_table <- data.frame(fasta_name,
                          group = genomes_organisation_data$Group[1],
                          genus,
                          species,
                          GC_genome.perc,
                          GC_chromosomes.perc,
                          GC_non_chromosomes.perc,
                          GC_centromeres.perc,
                          GC_pericentromeres.perc,
                          GC_arms.perc,
                          GC_exons_total.perc,
                          GC_TE_total.perc,
                          GC_repeats_total.perc,
                          GC_repeats_centromeric.perc,
                          GC_repeats_noncentromeric.perc
)

write.csv(x = grand_table, file = paste0(fasta_name, "_", genomes_organisation_data$Centromere.architecture[1], "_GC_grand_table.csv"), row.names = FALSE)
write.csv(x = grand_table, file = paste0("/home/pwlodzimierz/ToL/upload_files/13_GC_grand_tables/", fasta_name, "_", genomes_organisation_data$Centromere.architecture[1],"_GC_grand_table.csv"), row.names = FALSE)



