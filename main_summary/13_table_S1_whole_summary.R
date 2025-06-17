

taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 193
print(i)



.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))
suppressMessages(library(msa))
suppressMessages(library(seqinr)) 
suppressMessages(library(GenomicRanges))
suppressMessages(library(ggplot2))

ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}


filter_short_gaps <- 250

max_repeats_to_compare <- 1000

max_repeats_to_align <- 2500 # add this to the 10% of repeats that are used by default

do_GC <- TRUE
redo <- FALSE

setwd("/home/pwlodzimierz/ToL/git_ToL")
source("./aux_fun.R")

data_directories <- list.dirs(path = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs", recursive = FALSE, full.names = TRUE)
data_directories <- data_directories[!grepl(pattern = "templated_", data_directories)]
data_directories <- data_directories[grepl(pattern = ".fa", data_directories)]

assembly_files <- list.files(path = "/home/pwlodzimierz/ToL/Assemblies/fastas_2021_Michael", recursive = FALSE, full.names = TRUE)
assembly_files <- assembly_files[!grepl(".fai", assembly_files)]

if(i > length(data_directories)) stop()
if(grepl("Rosa_agrestis_DTOL_chrs_5n.fasta", data_directories[i])) stop()
if(grepl("Rosa_canina_DTOL_chrs_5n.fasta", data_directories[i])) stop()
if(grepl("drGeuRiva1.hap1.1.fa", data_directories[i])) stop()
if(grepl("drGeuRiva1.hap2.1.fa", data_directories[i])) stop()
if(grepl("drRosRugo1.1.fa", data_directories[i])) stop()
if(grepl("drRosSpin1.hap1.1.fa", data_directories[i])) stop()
if(grepl("GCA_000001405.29_GRCh38.p14_genomic.fa", data_directories[i])) stop()

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
  if(file.exists(paste0("/home/pwlodzimierz/ToL/upload_files/grand_tables/FULL_", fasta_name, "_", genomes_organisation_data$Centromere.architecture[1], "_grand_table.csv"))) {
    stop("File exists, not redoing")
  }
}

# print("reading fasta")
# fasta <- read.fasta(file = assembly_files[assembly_file])
# print("fasta read")

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
  
  max_repeats_to_compare <- 200
  
  max_repeats_to_align <- 50
}

if(i == 178) edta <- edta[-1,]
if(i == 153) edta <- edta[-1,]

edta$V5 <- as.numeric(edta$V5)
edta$V6 <- as.numeric(edta$V6)

print("annotations loaded")

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


# similarity within between files

if(nrow(centromeric_stallite_for_species) > 0) {
  within_between_table_1 <- read.csv(paste0("/home/pwlodzimierz/ToL/upload_files/37_similarity_values_within_between_chr/similarities_within_between_chrs_", fasta_name, "_", centromeric_stallite_for_species$Satellite_name[1], ".csv"))
  if(nrow(centromeric_stallite_for_species) > 1) {
    within_between_table_2 <- read.csv(paste0("/home/pwlodzimierz/ToL/upload_files/37_similarity_values_within_between_chr/similarities_within_between_chrs_", fasta_name, "_", centromeric_stallite_for_species$Satellite_name[2], ".csv"))
  }
}


### Data present:
# > str(repeats)
# 'data.frame':   99230 obs. of  16 variables:
#   $ seqID            : chr  "chr_1" "chr_1" "chr_1" "chr_1" ...
# $ arrayID          : int  1 1 1 1 1 1 1 1 1 1 ...
# $ start            : int  11 18 25 32 38 45 52 59 66 74 ...
# $ end              : int  17 24 31 37 44 51 58 65 73 80 ...
# $ strand           : chr  "-" "-" "-" "-" ...
# $ score            : num  1.00e-08 1.43e+01 1.43e+01 1.43e+01 1.43e+01 ...
# $ eval             : num  -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 ...
# $ width            : int  7 7 7 6 7 7 7 7 8 7 ...
# $ class            : chr  "7_4" "7_4" "7_4" "7_4" ...
# $ score_template   : num  -1 -1 -1 0.143 -1 ...
# $ sequence         : chr  "gtttagg" "gttttgg" "gttttgg" "gttagg" ...
# $ new_class        : chr  "7_4" "7_4" "7_4" "7_4" ...
# $ new_class_num_ID : int  5 5 5 5 5 5 5 5 5 5 ...
# $ is_on_chromosomes: logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
# $ is_centromeric   : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
# $ in_centromere    : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
# 
# > str(edta)
# 'data.frame':   31117 obs. of  24 variables:
#   $ V1 : int  1 2 3 4 5 6 7 8 9 10 ...
# $ V2 : chr  "chr_1" "chr_1" "chr_1" "chr_1" ...
# $ V3 : chr  "EDTA" "EDTA" "EDTA" "EDTA" ...
# $ V4 : chr  "helitron" "helitron" "helitron" "TE_unclass" ...
# $ V5 : num  20088 20334 41813 45458 47767 ...
# $ V6 : num  20333 21704 41906 45693 48025 ...
# $ V7 : chr  "." "." "." "." ...
# $ V8 : chr  "+" "+" "-" "+" ...
# $ V9 : chr  "." "." "." "." ...
# $ V10: chr  "TE_homo_1" "TE_homo_2" "TE_homo_3" "TE_homo_4" ...
# $ V11: chr  "ATREP3" "ATREP20" "TE_00000237" "TE_00000153" ...
# $ V12: chr  "Unspecified" "Unspecified" "DNA/Helitron" "Unknown" ...
# $ V13: chr  "SO:0000657" "SO:0000657" "SO:0000544" "SO:0000657" ...
# $ V14: num  0.918 0.952 0.769 0.915 0.957 0.941 0.763 0.833 0.962 0.807 ...
# $ V15: chr  "homology" "homology" "homology" "homology" ...
# $ V16: chr  NA NA NA NA ...
# $ V17: chr  NA NA NA NA ...
# $ V18: chr  "." "." "." "." ...
# $ V19: chr  "." "." "." "." ...
# $ V20: chr  "repeat_region" "repeat_region" "helitron" "repeat_region" ...
# $ V21: int  0 0 0 0 0 0 0 0 0 0 ...
# $ V22: int  246 1371 94 236 259 901 239 164 960 130 ...
# $ V23: num  0 0 0 0 0 0 0 0 0 0 ...
# $ V24: chr  "repeat_region/helitron/atrep3&unspecified/atrep3/step2" "repeat_region/helitron/atrep20&unspecified/atrep20/step2" "unaffected" "repeat_region/TE_unclass/te_00000153&unknown/te_00/step1" ...
# 
# 
# > str(classes)
# 'data.frame':   98 obs. of  8 variables:
#   $ class        : chr  "178_1" "500_2" "159_9" "497_6" ...
# $ count        : int  79791 3493 1099 342 12490 162 558 82 25 246 ...
# $ consensus    : chr  "caaccttcttcttgcttctcaaagctttcatggtgtagccaaagtccatatgagtctttggctttgtgtcttctaacaaggaaacactacttaggcttttaagatcgggtt"| __truncated__ "gttaagcgtgcttgggcgagagtagtactaggatgggtgacctcccgggaagtcctcgtgttgcatccctctttttttttttttttttttttttggttaaaactttatgac"| __truncated__ "acaaagctattgactgcttctaagcaattttttgttggttttagcctcttttgggagaaaatgggtataagtgttgtctaaacactcctaatccatctctaactcttataa"| __truncated__ "gtaggatttaggttttaaagttctatggtttagggtttagggtttagggttaagggtttagggttaagagtttatttagggttaggggttacggtttagggtttaggattt"| __truncated__ ...
# $ median_length: int  178 500 159 487 7 447 97 355 868 75 ...
# $ length_SD    : num  4.585 9.44 26.63 84.045 0.297 ...
# $ score        : num  0.0545 0.0336 0.0928 0.2158 0.0906 ...
# $ sum_coverage : int  14202798 1748300 174741 166554 88522 131382 54126 29110 34524 18450 ...
# $ num_ID       : int  1 2 3 4 5 6 7 8 9 10 ...
# 
# 
# > str(genes)
# 'data.frame':   144062 obs. of  9 variables:
#   $ V1: chr  "chr_1" "chr_1" "chr_1" "chr_1" ...
# $ V2: chr  "Helixer" "Helixer" "Helixer" "Helixer" ...
# $ V3: chr  "CDS" "CDS" "CDS" "CDS" ...
# $ V4: int  6834 7070 7546 7766 8234 8500 26657 27680 27890 28179 ...
# $ V5: int  6987 7350 7665 8155 8386 8691 27589 27793 28100 28573 ...
# $ V6: chr  "." "." "." "." ...
# $ V7: chr  "+" "+" "+" "+" ...
# $ V8: chr  "." "." "." "." ...
# $ V9: chr  "ID=ddAraThal4_chr_1_000001.1.CDS.1;Parent=ddAraThal4_chr_1_000001.1" "ID=ddAraThal4_chr_1_000001.1.CDS.2;Parent=ddAraThal4_chr_1_000001.1" "ID=ddAraThal4_chr_1_000001.1.CDS.3;Parent=ddAraThal4_chr_1_000001.1" "ID=ddAraThal4_chr_1_000001.1.CDS.4;Parent=ddAraThal4_chr_1_000001.1" ...
# 
# 
# > str(genome_metadata)
# 'data.frame':   23 obs. of  57 variables:
#   $ X                                                              : int  3753 3754 3755 3756 3757 3758 3759 3760 3761 3762 ...
# $ assembly.name                                                  : chr  "ddAraThal4.1.fa" "ddAraThal4.1.fa" "ddAraThal4.1.fa" "ddAraThal4.1.fa" ...
# $ chromosome.name                                                : chr  "chr_1" "chr_2" "chr_3" "chr_4" ...
# $ size                                                           : int  33257523 21973073 27693526 23271492 30562185 193407 168487 108111 71548 67749 ...
# $ is.chr                                                         : int  1 1 1 1 1 0 0 0 0 0 ...
# $ Ians                                                           : int  1 1 1 1 1 0 0 0 0 0 ...
# $ analysed_centromeric_Satellite_name                            : chr  "ddAraThal.178" "ddAraThal.178" "ddAraThal.178" "ddAraThal.178" ...
# $ analysed_centromeric_Trash_name                                : chr  "178_1" "178_1" "178_1" "178_1" ...
# $ Architecture                                                   : chr  "metacentric" "metacentric" "metacentric" "metacentric" ...
# $ centromere_start                                               : int  14310876 3947504 13391015 3294666 12351160 0 0 0 0 0 ...
# $ centromere_end                                                 : int  18681071 5818620 18125527 9039249 15911252 0 0 0 0 0 ...
# $ pericentromere_start                                           : int  13423218 2809344 11438631 2608013 10521409 0 0 0 0 0 ...
# $ pericentromere_end                                             : int  19433612 8127027 19566077 9930507 17235068 0 0 0 0 0 ...
# $ centromere_arrays_starts                                       : chr  "14310876;18569369" "3947504" "13391015;17131930;17573657;18109587" "3294666;5646676;8437528" ...
# $ centromere_arrays_ends                                         : chr  "17218293;18681071" "5818620" "17027242;17256990;17970155;18125527" "3864472;7659868;9039249" ...
# $ centromere_mean_TE_coverage                                    : num  74.6 57.8 96.8 82.7 98.3 ...
# $ pericentromere_mean_TE_coverage                                : num  54.5 71.1 70.3 58.2 52.7 ...
# $ arm_mean_TE_coverage                                           : num  9.54 13.74 9.89 12.22 10.62 ...
# $ centromere_mean_Gyspy_coverage                                 : num  61 46.5 91.9 72.8 96.6 ...
# $ pericentromere_mean_Gyspy_coverage                             : num  26.1 32.1 36.2 32.3 29.9 ...
# $ arm_mean_Gyspy_coverage                                        : num  0.521 1.16 0.532 0.861 0.658 ...
# $ centromere_mean_Copia_coverage                                 : num  1.8 3.91 0 1.76 0 ...
# $ pericentromere_mean_Copia_coverage                             : num  2.72 3.68 3.9 2.26 3.32 ...
# $ arm_mean_Copia_coverage                                        : num  0.962 1.19 1.308 1.249 0.935 ...
# $ centromere_mean_Class1_LTR_coverage                            : num  63 50.4 94.8 75.7 97.2 ...
# $ pericentromere_mean_Class1_LTR_coverage                        : num  30.6 37.6 41.9 36.4 34.8 ...
# $ arm_mean_Class1_LTR_coverage                                   : num  1.82 3.89 2.2 4.35 1.95 ...
# $ centromere_mean_Class1_nonLTR_coverage                         : num  2.65 0 0 1.03 0 ...
# $ pericentromere_mean_Class1_nonLTR_coverage                     : num  1.2 4.58 1.94 1.18 0.72 ...
# $ arm_mean_Class1_nonLTR_coverage                                : num  0.726 1.175 0.555 0.814 0.807 ...
# $ centromere_mean_Class2_TIR_coverage                            : num  3.6711 7.0108 0.7205 4.9889 0.0632 ...
# $ pericentromere_mean_Class2_TIR_coverage                        : num  16.1 20 19.1 12.3 11 ...
# $ arm_mean_Class2_TIR_coverage                                   : num  3.27 3.94 3.28 3.41 3.57 ...
# $ centromere_mean_Class2_nonTIR_coverage                         : num  4.877 0 0 0.331 0 ...
# $ pericentromere_mean_Class2_nonTIR_coverage                     : num  3.69 6.56 6.1 4.29 4.79 ...
# $ arm_mean_Class2_nonTIR_coverage                                : num  2.99 3.98 2.84 3.06 3.52 ...
# $ centromere_mean_unspecified_coverage                           : num  0 0 0.0272 0.0338 0 ...
# $ pericentromere_mean_unspecified_coverage                       : num  0.278 0.35 0.349 0.196 0.214 ...
# $ arm_mean_unspecified_coverage                                  : num  0.552 0.43 0.656 0.42 0.443 ...
# $ centromere_GC                                                  : num  0.376 0.378 0.386 0.379 0.381 ...
# $ pericentromere_GC                                              : num  0.367 0.373 0.374 0.388 0.387 ...
# $ arm_GC                                                         : num  0.358 0.358 0.362 0.362 0.357 ...
# $ centromeric_arrays_number                                      : int  2 1 4 3 1 0 0 0 0 0 ...
# $ mean_centromeric_arrays_size                                   : num  1509560 1871116 1043431 1061573 3560092 ...
# $ mean_centromeric_arrays_size_sd                                : num  1976869 NA 1735920 824281 NA ...
# $ gaps_in_centromeric_arrays                                     : int  20 6 35 35 38 0 0 0 0 0 ...
# $ gaps_in_centromeric_arrays_mean_size                           : num  16541 21576 11746 10234 10859 ...
# $ mean_centromeric_repeat_position_along_chromosome_normalised   : num  47.3 22 55.5 29.1 46.2 ...
# $ mean_centromeric_repeat_position_along_chromosome_normalised_sd: num  2.51 2.32 4.17 6.18 3.16 ...
# $ fraction_of_plus_strand_centromeric_repeats                    : num  0.106 99.98 96.528 73.632 2.845 ...
# $ mean_centromeric_repeat_width                                  : num  178 178 177 178 177 ...
# $ mean_centromeric_repeat_width_sd                               : num  4.4 3.24 5.32 3.6 5.15 ...
# $ mean_centromeric_repeat_pairwise_similarity                    : num  94.4 92.4 94.6 90.9 92.5 ...
# $ mean_centromeric_repeat_pairwise_similarity_sd                 : num  3.4 4.21 3.51 5.38 4.45 ...
# $ percent_centromeric_repeat_in_all_class_repeats                : num  100 99.9 100 100 100 ...
# $ relative_centromere_size_vs_chr                                : num  9.08 8.52 15.07 13.69 11.65 ...
# $ relative_pericentromere_size_vs_chr                            : num  8.99 15.69 14.28 17.78 10.32 ...
# 
# 
# > str(gaps)
# 'data.frame':   135 obs. of  5 variables:
#   $ start          : int  1916494 3348219 3359615 3407381 3413902 3417498 3433871 3444868 3467212 3492567 ...
# $ end            : int  1917021 3356767 3402383 3411174 3415904 3428150 3440299 3464468 3484677 3497406 ...
# $ repeat_array_id: chr  "array_chr_2_artha.178_1" "array_chr_4_artha.178_1" "array_chr_4_artha.178_1" "array_chr_4_artha.178_1" ...
# $ species        : chr  "ddAraThal4.1.fa" "ddAraThal4.1.fa" "ddAraThal4.1.fa" "ddAraThal4.1.fa" ...
# $ chromosome     : chr  "chr_2" "chr_4" "chr_4" "chr_4" ...
# 
# 
# > str(centromeric_stallite_for_species)
# 'data.frame':   1 obs. of  22 variables:
#   $ Kolumna1                                   : int 34
# $ Clade                                      : chr "Dicot"
# $ group2                                     : chr "plant"
# $ Species                                    : chr "Arabidopsis_thaliana"
# $ Genome                                     : chr "ddAraThal4.1.fa"
# $ Trash_name_aug2024runs                     : chr "178_1"
# $ Satellite_name                             : chr "ddAraThal.178"
# $ Architecture                               : chr "Satellite"
# $ Typic                                      : chr "Monotypic"
# $ same_species_as_above                      : int 0
# $ TRASH_name_dec2024runs                     : chr "178_1"
# $ Satellite_name_current                     : chr "ddAraThal.178"
# $ mean_repeat_size                           : num 178
# $ repeat_no                                  : int 79791
# $ repeat_bp                                  : int 14168775
# $ repeat_size_SD_in_perc                     : num 2.58
# $ consensus                                  : logi NA
# $ mean_similairty_within_chromosomes         : logi NA
# $ mean_similairty_between_chromosomes        : logi NA
# $ holocentrics.mean_similairty_between_arrays: logi NA
# $ holocentrics.mean_array_size               : logi NA
# $ repeat_no.including_non_chr_seq            : int 79791

# > str(genomes_organisation_data)
# 'data.frame':   1 obs. of  7 variables:
#   $ fasta                  : chr "ddAraThal4.1.fa"
# $ Group                  : chr "Dicots"
# $ Genus                  : chr "Arabidopsis"
# $ Species                : chr "thaliana"
# $ Chromosomes            : int 5
# $ Genome.Size            : num 1.37e+08
# $ Centromere.architecture: chr "Satellite"


### genomes_organisation_data
# one column with the organisation type

### centromeric_stallite_for_species
# 1 or 2 columns for satellite species with the family/families details, nothing worth copying over

### genome_metadata
# A column per sequence in the fasta file, can be used for centromere and percnetromere starts and ends

### gaps
# gaps in the centromeric arrays, $repeat_array_id can identify which repeat the gap is found in

### classes
# useless for this

### genes
# exons

### edta
# transposons

### repeats
# all repeats



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

chromosomes_with_centromeric_repeats.no <- NA
if(centromeric_satellites.1_yes_0_no == 1) {
  chromosomes_with_centromeric_repeats.no <- 0
  for(j in seq_along(chromosomes)) {
    if(nrow(repeats[repeats$seqID == chromosomes[j] & repeats$is_centromeric,]) > 0) {
      chromosomes_with_centromeric_repeats.no <- chromosomes_with_centromeric_repeats.no + 1
    }
  }
}

print("2")

size_genome.bp <- sum(chromosomes_lengths) + sum(non_chr_lengths)

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
for(j in which(genome_metadata$is.chr == 1)) {
  gen_meta_temp <- genome_metadata[j,]
  peri_starts <- as.numeric(strsplit(as.character(gen_meta_temp$pericentromere_start[1]), split = ";")[[1]])
  peri_ends <- as.numeric(strsplit(as.character(gen_meta_temp$pericentromere_end[1]), split = ";")[[1]])
  size_pericentromeres.bp <- size_pericentromeres.bp + sum(width(IRanges(peri_starts, peri_ends)))
}
if(centromeric_satellites.1_yes_0_no == 1) {
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


genes_chr.no <- nrow(genes_full[genes_full$V1 %in% chromosomes & genes_full$V3 == "gene",])

genes_chr.bp <- sum(genes_full$V5[genes_full$V1 %in% chromosomes & genes_full$V3 == "gene"] - genes_full$V4[genes_full$V1 %in% chromosomes & genes_full$V3 == "gene"])

genes_chr.perc <- 100 * genes_chr.bp / size_chromosomes.bp

exons_chr.no <- nrow(genes[genes$V1 %in% chromosomes,])

exons_chr.bp <- sum(genes$V5[genes$V1 %in% chromosomes] - genes$V4[genes$V1 %in% chromosomes])

exons_chr.perc <- 100 * exons_chr.bp / size_chromosomes.bp

genes$in_centromere <- FALSE
if(centromeric_satellites.1_yes_0_no == 1) {
  for(j in which(genome_metadata$is.chr == 1)) {
    
    arrays_temp <- centromeric_arrays[centromeric_arrays$seqID == genome_metadata$chromosome.name[j],]
    arrays_temp <- arrays_temp[arrays_temp$in_centromere,]
    
    genes_temp <- genes[genes$V1 == genome_metadata$chromosome.name[j],]
    
    genes[genes$V1 == genome_metadata$chromosome.name[j],]$in_centromere[subjectHits(findOverlaps(query = IRanges(arrays_temp$start, arrays_temp$end),
                                                                                                  subject = IRanges(genes_temp$V4, genes_temp$V5)))] <- TRUE
  }
}
exons_centromeres.no <- NA
exons_centromeres.bp <- NA
exons_centromeres.perc <- NA
if(centromeric_satellites.1_yes_0_no == 1) {
  exons_centromeres.no <- nrow(genes[genes$in_centromere,])
  
  exons_centromeres.bp <- sum(genes$V5[genes$in_centromere] - genes$V4[genes$in_centromere])
  
  exons_centromeres.perc <- 100 * exons_centromeres.bp / size_centromeres.bp
}

print("7")

genes$in_pericentromere <- FALSE
for(j in which(genome_metadata$is.chr == 1)) {
  
  genes_temp <- genes[genes$V1 == genome_metadata$chromosome.name[j],]
  
  gen_meta_temp <- genome_metadata[genome_metadata$chromosome.name == genome_metadata$chromosome.name[j],]
  peri_starts <- as.numeric(strsplit(as.character(gen_meta_temp$pericentromere_start[1]), split = ";")[[1]])
  peri_ends <- as.numeric(strsplit(as.character(gen_meta_temp$pericentromere_end[1]), split = ";")[[1]])
  if(length(peri_ends) == 1) {
    if(peri_ends == 0) next
  } 
  
  genes[genes$V1 == genome_metadata$chromosome.name[j],]$in_pericentromere[subjectHits(findOverlaps(query = IRanges(peri_starts, peri_ends),
                                                                                                    subject = IRanges(genes_temp$V4, genes_temp$V5)))] <- TRUE
}
genes$in_pericentromere[genes$in_centromere] = FALSE

exons_pericentromeres.no <- nrow(genes[genes$in_pericentromere,])

exons_pericentromeres.bp <- sum(genes$V5[genes$in_pericentromere] - genes$V4[genes$in_pericentromere])

exons_pericentromeres.perc <- 100 * exons_pericentromeres.bp / size_pericentromeres.bp

exons_arms.no <- nrow(genes[!genes$in_pericentromere & !genes$in_centromere,])

exons_arms.bp <- sum(genes$V5[!genes$in_pericentromere & !genes$in_centromere] - genes$V4[!genes$in_pericentromere & !genes$in_centromere])

exons_arms.perc <- 100 * exons_arms.bp / (size_chromosomes.bp - size_pericentromeres.bp - size_centromeres.bp)

TE_total.no <- nrow(edta)

TE_total.bp <- sum(edta$V6 - edta$V5)

print("8")
TE_total.perc <- 100 * TE_total.bp / size_chromosomes.bp

edta$extra_class = ""
edta$extra_class[grep(LTR_names, edta$V4, ignore.case = TRUE, perl = TRUE)] = "C1_LTR"
edta$extra_class[grep(nonLTR_names, edta$V4, ignore.case = TRUE, perl = TRUE)] = "C1_nonLTR"
edta$extra_class[grep(TIR_names, edta$V4, ignore.case = TRUE, perl = TRUE)] = "C2_TIR"
edta$extra_class[grep(nonTIR_names, edta$V4, ignore.case = TRUE, perl = TRUE)] = "C2_nonTIR"

TE_C1_LTR.no <- nrow(edta[edta$extra_class == "C1_LTR",])

TE_C1_LTR.bp <- sum(edta$V6[edta$extra_class == "C1_LTR"] - edta$V5[edta$extra_class == "C1_LTR"])

TE_C1_LTR.perc <- 100 * TE_C1_LTR.bp / size_chromosomes.bp

TE_C1_nonLTR.no <- nrow(edta[edta$extra_class == "C1_nonLTR",])

TE_C1_nonLTR.bp <- sum(edta$V6[edta$extra_class == "C1_nonLTR"] - edta$V5[edta$extra_class == "C1_nonLTR"])

TE_C1_nonLTR.perc <- 100 * TE_C1_nonLTR.bp / size_chromosomes.bp

TE_C2_TIR.no <- nrow(edta[edta$extra_class == "C2_TIR",])

TE_C2_TIR.bp <- sum(edta$V6[edta$extra_class == "C2_TIR"] - edta$V5[edta$extra_class == "C2_TIR"])

TE_C2_TIR.perc <- 100 * TE_C2_TIR.bp / size_chromosomes.bp

TE_C2_nonTIR.no <- nrow(edta[edta$extra_class == "C2_nonTIR",])

TE_C2_nonTIR.bp <- sum(edta$V6[edta$extra_class == "C2_nonTIR"] - edta$V5[edta$extra_class == "C2_nonTIR"])

TE_C2_nonTIR.perc <- 100 * TE_C2_nonTIR.bp / size_chromosomes.bp


LINE_names <- "line"
edta$extra_class[grep(LINE_names, edta$V4, ignore.case = TRUE, perl = TRUE)] = "LINE"

TE_LINE.no <- nrow(edta[edta$extra_class == "LINE",])

TE_LINE.bp <- sum(edta$V6[edta$extra_class == "LINE"] - edta$V5[edta$extra_class == "LINE"])

TE_LINE.perc <- 100 * TE_LINE.bp / size_chromosomes.bp


HELITRON_names <- "helitron"
edta$extra_class[grep(HELITRON_names, edta$V4, ignore.case = TRUE, perl = TRUE)] = "HELITRON"

TE_HELITRON.no <- nrow(edta[edta$extra_class == "HELITRON",])

TE_HELITRON.bp <- sum(edta$V6[edta$extra_class == "HELITRON"] - edta$V5[edta$extra_class == "HELITRON"])

TE_HELITRON.perc <- 100 * TE_HELITRON.bp / size_chromosomes.bp



edta$in_centromere <- FALSE
if(centromeric_satellites.1_yes_0_no == 1) {
  for(j in which(genome_metadata$is.chr == 1)) {
    
    arrays_temp <- centromeric_arrays[centromeric_arrays$seqID == genome_metadata$chromosome.name[j],]
    arrays_temp <- arrays_temp[arrays_temp$in_centromere,]
    
    edta_temp <- edta[edta$V2 == genome_metadata$chromosome.name[j],]
    
    edta[edta$V2 == genome_metadata$chromosome.name[j],]$in_centromere[subjectHits(findOverlaps(query = IRanges(arrays_temp$start, arrays_temp$end),
                                                                                                subject = IRanges(edta_temp$V5, edta_temp$V6)))] <- TRUE
  }
}
edta$in_pericentromere <- FALSE
for(j in which(genome_metadata$is.chr == 1)) {
  
  edta_temp <- edta[edta$V2 == genome_metadata$chromosome.name[j],]
  
  gen_meta_temp <- genome_metadata[genome_metadata$chromosome.name == genome_metadata$chromosome.name[j],]
  peri_starts <- as.numeric(strsplit(as.character(gen_meta_temp$pericentromere_start[1]), split = ";")[[1]])
  peri_ends <- as.numeric(strsplit(as.character(gen_meta_temp$pericentromere_end[1]), split = ";")[[1]])
  if(length(peri_ends) == 1) {
    if(peri_ends == 0) next
  } 
  
  edta[edta$V2 == genome_metadata$chromosome.name[j],]$in_pericentromere[subjectHits(findOverlaps(query = IRanges(peri_starts, peri_ends),
                                                                                                  subject = IRanges(edta_temp$V5, edta_temp$V6)))] <- TRUE
}
edta$in_pericentromere[edta$in_centromere] = FALSE


TE_total_centromere.perc <- 100 * sum(edta$V6[edta$in_centromere] - edta$V5[edta$in_centromere]) / size_centromeres.bp

print("9")
TE_total_pericentromere.perc <- 100 * sum(edta$V6[edta$in_pericentromere] - edta$V5[edta$in_pericentromere]) / size_pericentromeres.bp

TE_total_arms.perc <- 100 * sum(edta$V6[!edta$in_pericentromere & !edta$in_pericentromere] - edta$V5[!edta$in_pericentromere & !edta$in_pericentromere]) / (size_chromosomes.bp - size_pericentromeres.bp - size_centromeres.bp)

TE_COPIA_centromere.perc <- 100 * sum(edta$V6[edta$in_centromere & grepl("copia", edta$V4, ignore.case = TRUE)] -
                                        edta$V5[edta$in_centromere & grepl("copia", edta$V4, ignore.case = TRUE)]) / size_centromeres.bp

TE_COPIA_pericentromere.perc <- 100 * sum(edta$V6[edta$in_pericentromere & grepl("copia", edta$V4, ignore.case = TRUE)] -
                                            edta$V5[edta$in_pericentromere & grepl("copia", edta$V4, ignore.case = TRUE)]) / size_pericentromeres.bp

TE_COPIA_arms.perc <- 100 * sum(edta$V6[!edta$in_pericentromere & !edta$in_pericentromere & grepl("copia", edta$V4, ignore.case = TRUE)] -
                                  edta$V5[!edta$in_pericentromere & !edta$in_pericentromere & grepl("copia", edta$V4, ignore.case = TRUE)]) / (size_chromosomes.bp - size_pericentromeres.bp - size_centromeres.bp)

TE_GYPSY_centromere.perc <- 100 * sum(edta$V6[edta$in_centromere & grepl("gypsy", edta$V4, ignore.case = TRUE)] -
                                        edta$V5[edta$in_centromere & grepl("gypsy", edta$V4, ignore.case = TRUE)]) / size_centromeres.bp

TE_GYPSY_pericentromere.perc <- 100 * sum(edta$V6[edta$in_pericentromere & grepl("gypsy", edta$V4, ignore.case = TRUE)] -
                                            edta$V5[edta$in_pericentromere & grepl("gypsy", edta$V4, ignore.case = TRUE)]) / size_pericentromeres.bp

TE_GYPSY_arms.perc <- 100 * sum(edta$V6[!edta$in_pericentromere & !edta$in_pericentromere & grepl("gypsy", edta$V4, ignore.case = TRUE)] -
                                  edta$V5[!edta$in_pericentromere & !edta$in_pericentromere & grepl("gypsy", edta$V4, ignore.case = TRUE)]) / (size_chromosomes.bp - size_pericentromeres.bp - size_centromeres.bp)

TE_HELITRON_centromere.perc <- 100 * sum(edta$V6[edta$in_centromere & grepl("helitron", edta$V4, ignore.case = TRUE)] -
                                           edta$V5[edta$in_centromere & grepl("helitron", edta$V4, ignore.case = TRUE)]) / size_centromeres.bp

TE_HELITRON_pericentromere.perc <- 100 * sum(edta$V6[edta$in_pericentromere & grepl("helitron", edta$V4, ignore.case = TRUE)] -
                                               edta$V5[edta$in_pericentromere & grepl("helitron", edta$V4, ignore.case = TRUE)]) / size_pericentromeres.bp

print("10")
TE_HELITRON_arms.perc <- 100 * sum(edta$V6[!edta$in_pericentromere & !edta$in_pericentromere & grepl("helitron", edta$V4, ignore.case = TRUE)] -
                                     edta$V5[!edta$in_pericentromere & !edta$in_pericentromere & grepl("helitron", edta$V4, ignore.case = TRUE)]) / (size_chromosomes.bp - size_pericentromeres.bp - size_centromeres.bp)

TE_C1_LTR_centromere.perc <- 100 * sum(edta$V6[edta$in_centromere & edta$extra_class == "C1_LTR"] -
                                         edta$V5[edta$in_centromere & edta$extra_class == "C1_LTR"]) / size_centromeres.bp

TE_C1_LTR_pericentromere.perc <- 100 * sum(edta$V6[edta$in_pericentromere & edta$extra_class == "C1_LTR"] -
                                             edta$V5[edta$in_pericentromere & edta$extra_class == "C1_LTR"]) / size_pericentromeres.bp

TE_C1_LTR_arms.perc <- 100 * sum(edta$V6[!edta$in_pericentromere & !edta$in_pericentromere & edta$extra_class == "C1_LTR"] -
                                   edta$V5[!edta$in_pericentromere & !edta$in_pericentromere & edta$extra_class == "C1_LTR"]) / (size_chromosomes.bp - size_pericentromeres.bp - size_centromeres.bp)

TE_C1_nonLTR_centromere.perc <- 100 * sum(edta$V6[edta$in_centromere & edta$extra_class == "C1_nonLTR"] -
                                            edta$V5[edta$in_centromere & edta$extra_class == "C1_nonLTR"]) / size_centromeres.bp

TE_C1_nonLTR_pericentromere.perc <- 100 * sum(edta$V6[edta$in_pericentromere & edta$extra_class == "C1_nonLTR"] -
                                                edta$V5[edta$in_pericentromere & edta$extra_class == "C1_nonLTR"]) / size_pericentromeres.bp

TE_C1_nonLTR_arms.perc <- 100 * sum(edta$V6[!edta$in_pericentromere & !edta$in_pericentromere & edta$extra_class == "C1_nonLTR"] -
                                      edta$V5[!edta$in_pericentromere & !edta$in_pericentromere & edta$extra_class == "C1_nonLTR"]) / (size_chromosomes.bp - size_pericentromeres.bp - size_centromeres.bp)

TE_C2_TIR_centromere.perc <- 100 * sum(edta$V6[edta$in_centromere & edta$extra_class == "C2_TIR"] -
                                         edta$V5[edta$in_centromere & edta$extra_class == "C2_TIR"]) / size_centromeres.bp

TE_C2_TIR_pericentromere.perc <- 100 * sum(edta$V6[edta$in_pericentromere & edta$extra_class == "C2_TIR"] -
                                             edta$V5[edta$in_pericentromere & edta$extra_class == "C2_TIR"]) / size_pericentromeres.bp

TE_C2_TIR_arms.perc <- 100 * sum(edta$V6[!edta$in_pericentromere & !edta$in_pericentromere & edta$extra_class == "C2_TIR"] -
                                   edta$V5[!edta$in_pericentromere & !edta$in_pericentromere & edta$extra_class == "C2_TIR"]) / (size_chromosomes.bp - size_pericentromeres.bp - size_centromeres.bp)

TE_C2_nonTIR_centromere.perc <- 100 * sum(edta$V6[edta$in_centromere & edta$extra_class == "C2_nonTIR"] -
                                            edta$V5[edta$in_centromere & edta$extra_class == "C2_nonTIR"]) / size_centromeres.bp

TE_C2_nonTIR_pericentromere.perc <- 100 * sum(edta$V6[edta$in_pericentromere & edta$extra_class == "C2_nonTIR"] -
                                                edta$V5[edta$in_pericentromere & edta$extra_class == "C2_nonTIR"]) / size_pericentromeres.bp

TE_C2_nonTIR_arms.perc <- 100 * sum(edta$V6[!edta$in_pericentromere & !edta$in_pericentromere & edta$extra_class == "C2_nonTIR"] -
                                      edta$V5[!edta$in_pericentromere & !edta$in_pericentromere & edta$extra_class == "C2_nonTIR"]) / (size_chromosomes.bp - size_pericentromeres.bp - size_centromeres.bp)

print("11")
repeats_total.no <- nrow(repeats)

repeats_total.bp <- sum(repeats$width)

repeats_total.perc <- 100 * repeats_total.bp / sum(genome_metadata$size)

centromeric_families.no <- nrow(centromeric_stallite_for_species)

repeats_centromeric.no <- NA
repeats_centromeric.bp <- NA
repeats_centromeric.perc <- NA

if(centromeric_satellites.1_yes_0_no == 1) {
  repeats_centromeric.no <- sum(repeats$is_centromeric)
  
  repeats_centromeric.bp <- sum(repeats$width[repeats$is_centromeric])
  
  repeats_centromeric.perc <- 100 * repeats_centromeric.bp / sum(genome_metadata$size)
}


repeats_noncentromeric.no <- repeats_total.no - repeats_centromeric.no

repeats_noncentromeric.bp <- repeats_total.bp - repeats_centromeric.bp

repeats_noncentromeric.perc <- repeats_total.perc - repeats_centromeric.perc

print("12")
repeats_centromeric_on_chromosomes.bp <- NA
if(centromeric_satellites.1_yes_0_no == 1) {
  repeats_centromeric_on_chromosomes.bp <- sum(repeats$width[repeats$is_centromeric & repeats$seqID %in% chromosomes])
}

repeats_centromeric_not_on_chromosomes.bp <- NA
if(centromeric_satellites.1_yes_0_no == 1) {
  repeats_centromeric_not_on_chromosomes.bp <- sum(repeats$width[!repeats$is_centromeric & repeats$seqID %in% chromosomes])
}

repeats_centromeric_not_on_chromosomes_vs_all_centromeric.perc <- NA
if(centromeric_satellites.1_yes_0_no == 1) {
  repeats_centromeric_not_on_chromosomes_vs_all_centromeric.perc <- 100 * repeats_centromeric_on_chromosomes.bp / repeats_centromeric.bp
}

centromeric_repeats_fraction_on_dominant_strand.perc <- NA

if(centromeric_satellites.1_yes_0_no == 1) {
  reps_on_dom_temp <- 0
  reps_bp_temp <- 0
  for(j in which(genome_metadata$is.chr == 1)) {
    reps_temp <- repeats[repeats$is_centromeric & repeats$seqID == genome_metadata$chromosome.name[j],]
    if(nrow(reps_temp) == 0) next
    perc <- 100 * sum(reps_temp$strand == "+") / nrow(reps_temp)
    if(sum(reps_temp$strand == "+") == 0) perc <- 0
    if(length(perc) == 0) next
    if(perc < 50) perc = 100 - perc
    reps_on_dom_temp <- reps_on_dom_temp + perc * sum(repeats$width)
    reps_bp_temp <- reps_bp_temp + sum(repeats$width)
  }
  centromeric_repeats_fraction_on_dominant_strand.perc <- reps_on_dom_temp / reps_bp_temp
  
}


centromeric_repeats_in_centromere.perc <- NA
if(centromeric_satellites.1_yes_0_no == 1) {
  centromeric_repeats_in_centromere.perc <- 100 * sum(repeats$width[repeats$in_centromere]) / size_centromeres.bp
}

centromeric_repeat_1_name <- NA
if(centromeric_families.no > 0) {
  centromeric_repeat_1_name <- centromeric_stallite_for_species$Satellite_name[1]
}

centromeric_repeat_1_TRASH_identifier <- NA
if(centromeric_families.no > 0) {
  centromeric_repeat_1_TRASH_identifier <- strsplit(centromeric_stallite_for_species$TRASH_name_dec2024runs[1], split = ";")[[1]]
}

centromeric_repeat_1_mean_size.bp <- NA
if(centromeric_families.no > 0) {
  repeats_1 <- repeats[repeats$is_centromeric & repeats$is_on_chromosomes & repeats$new_class %in% centromeric_repeat_1_TRASH_identifier,]
  centromeric_repeat_1_mean_size.bp <- mean(repeats_1$width)
}

centromeric_repeat_1_consensus_unambiguous <- NA
centromeric_repeat_1_consensus_ambiguous <- NA

if(centromeric_families.no > 0) {
  if((nrow(repeats_1) %/% 100 + max_repeats_to_align) > nrow(repeats_1)) {
    repeats_to_align = seq_len(nrow(repeats_1))
  } else if((nrow(repeats_1) %/% 100 + max_repeats_to_align) < max_repeats_to_align) {
    repeats_to_align = seq_len(nrow(repeats_1))
  } else {
    repeats_to_align = sample(seq_len(nrow(repeats_1)), (nrow(repeats_1) %/% 100 + max_repeats_to_align), replace = FALSE)
  }
  sequences_to_align <- repeats_1$sequence[repeats_to_align]
  if(length(sequences_to_align) == 0) {
    cat("\n\n\n\n\n\n\n\n\n\n")
    stop(paste0(assembly_file, " did not find repeats in one of the classes: ", centromeric_repeat_1_TRASH_identifier, ", investigate"))
  } else if(length(sequences_to_align) == 1) {
    centromeric_repeat_1_consensus_unambiguous <- sequences_to_align
    centromeric_repeat_1_consensus_ambiguous <- sequences_to_align
  } else {
    a <- capture.output({alignment_matrix = msa(sequences_to_align, method = "ClustalOmega", type = "dna")})
    centromeric_repeat_1_consensus_unambiguous <- consensus_N(alignment_matrix, centromeric_repeat_1_mean_size.bp)
    
    consensus_ambi_1 <- consensusString(alignment_matrix)
    
    cm <- consensusMatrix(alignment_matrix@unmasked)
    coverage <- colSums(cm[setdiff(rownames(cm), "-"), ])
    
    # Keep positions with max coverage up to specified length
    max_cov_positions <- order(coverage, decreasing = TRUE)[1:centromeric_repeat_1_mean_size.bp]
    centromeric_repeat_1_consensus_ambiguous <- tolower(paste(strsplit(consensus_ambi_1, "")[[1]][sort(max_cov_positions)], collapse = ""))
    
    
  }
}

print("13")


centromeric_repeat_1_mean_size_SD.bp <- NA
if(centromeric_families.no > 0) {
  centromeric_repeat_1_mean_size_SD.bp <- sd(repeats_1$width)
}

centromeric_repeat_1_mean_size_SD.perc <- NA
if(centromeric_families.no > 0) {
  centromeric_repeat_1_mean_size_SD.perc <- 100 * sd(repeats_1$width) / centromeric_repeat_1_mean_size.bp
}

centromeric_repeat_1_mean_pairwise_similarity.perc <- NA
if(centromeric_families.no > 0) {
  
  sample_repeats <- sample(1 : nrow(repeats_1), max_repeats_to_compare, replace = T)
  
  similarity_1_matrix <- adist(repeats_1$sequence[sample_repeats])
  centromeric_repeat_1_mean_pairwise_similarity.perc <- 100 * (1 - mean(similarity_1_matrix) / centromeric_repeat_1_mean_size.bp)
  
}


centromeric_repeat_1_mean_pairwise_similarity_SD.perc <- NA
if(centromeric_families.no > 0) {
  centromeric_repeat_1_mean_pairwise_similarity_SD.perc <- 100 * (sd(similarity_1_matrix) / centromeric_repeat_1_mean_size.bp)
}



centromeric_repeat_1_similarity_between_chromosomes_mean_pairwise.perc <- NA
if(centromeric_families.no > 0) {
  centromeric_repeat_1_similarity_between_chromosomes_mean_pairwise.perc <- mean(within_between_table_1$similarity[within_between_table_1$method == "diff_chr"])
}

centromeric_repeat_1_similarity_between_chromosomes_mean_pairwise_SD.perc_point <- NA
if(centromeric_families.no > 0) {
  centromeric_repeat_1_similarity_between_chromosomes_mean_pairwise_SD.perc_point <- sd(within_between_table_1$similarity[within_between_table_1$method == "diff_chr"])
}

centromeric_repeat_1_similarity_within_chromosomes_mean_pairwise.perc <- NA
if(centromeric_families.no > 0) {
  centromeric_repeat_1_similarity_within_chromosomes_mean_pairwise.perc <- mean(within_between_table_1$similarity[within_between_table_1$method == "same_chr"])
}

centromeric_repeat_1_similarity_within_chromosomes_mean_pairwise_SD.perc_point <- NA
if(centromeric_families.no > 0) {
  centromeric_repeat_1_similarity_within_chromosomes_mean_pairwise_SD.perc_point <- sd(within_between_table_1$similarity[within_between_table_1$method == "same_chr"])
}

centromeric_repeat_1_GC.perc <- NA
if(do_GC) {
  if(centromeric_families.no > 0) {
    if(nrow(repeats_1) != 0) {
      if("sequence" %in% names(repeats_1)) {
        if(repeats_1$sequence[1] != "") {
          centromeric_repeat_1_GC.perc <- 100 * GC(strsplit(paste(repeats_1$sequence, collapse = ""), split = "")[[1]])
        }
      }
    }
  }
}






centromeric_repeat_2_name <- NA
if(centromeric_families.no > 1) {
  centromeric_repeat_2_name <- centromeric_stallite_for_species$Satellite_name[2]
}

centromeric_repeat_2_TRASH_identifier <- NA
if(centromeric_families.no > 1) {
  centromeric_repeat_2_TRASH_identifier <- strsplit(centromeric_stallite_for_species$TRASH_name_dec2024runs[2], split = ";")[[1]]
}

centromeric_repeat_2_mean_size.bp <- NA
if(centromeric_families.no > 1) {
  repeats_2 <- repeats[repeats$is_centromeric & repeats$is_on_chromosomes & repeats$new_class %in% centromeric_repeat_2_TRASH_identifier,]
  centromeric_repeat_2_mean_size.bp <- mean(repeats_2$width)
}

centromeric_repeat_2_consensus_unambiguous <- NA
centromeric_repeat_2_consensus_ambiguous <- NA
if(centromeric_families.no > 1) {
  if((nrow(repeats_2) %/% 100 + max_repeats_to_align) > nrow(repeats_2)) {
    repeats_to_align = seq_len(nrow(repeats_2))
  } else if(nrow(repeats_2) < max_repeats_to_align) {
    repeats_to_align = seq_len(nrow(repeats_2))
  } else {
    repeats_to_align = sample(seq_len(nrow(repeats_2)), (nrow(repeats_2) %/% 100 + max_repeats_to_align), replace = FALSE)
  }
  sequences_to_align <- repeats_2$sequence[repeats_to_align]
  if(length(sequences_to_align) == 0) {
    cat("\n\n\n\n\n\n\n\n\n\n")
    stop(paste0(assembly_file, " did not find repeats in one of the classes: ", centromeric_repeat_2_TRASH_identifier, ", investigate"))
  } else if(length(sequences_to_align) == 1) {
    centromeric_repeat_2_consensus_unambiguous <- sequences_to_align
    centromeric_repeat_2_consensus_ambiguous <- sequences_to_align
  } else {
    a <- capture.output({alignment_matrix = msa(sequences_to_align, method = "ClustalOmega", type = "dna")})
    centromeric_repeat_2_consensus_unambiguous <- consensus_N(alignment_matrix, centromeric_repeat_2_mean_size.bp)
    
    consensus_ambi_2 <- consensusString(alignment_matrix)
    
    cm <- consensusMatrix(alignment_matrix@unmasked)
    coverage <- colSums(cm[setdiff(rownames(cm), "-"), ])
    
    # Keep positions with max coverage up to specified length
    max_cov_positions <- order(coverage, decreasing = TRUE)[1:centromeric_repeat_2_mean_size.bp]
    centromeric_repeat_2_consensus_ambiguous <- tolower(paste(strsplit(consensus_ambi_2, "")[[1]][sort(max_cov_positions)], collapse = ""))
    
    
  }
}

print("16")

centromeric_repeat_2_mean_size_SD.bp <- NA
if(centromeric_families.no > 1) {
  centromeric_repeat_2_mean_size_SD.bp <- sd(repeats_2$width)
}

centromeric_repeat_2_mean_size_SD.perc <- NA
if(centromeric_families.no > 1) {
  centromeric_repeat_2_mean_size_SD.perc <- 100 * sd(repeats_2$width) / centromeric_repeat_2_mean_size.bp
}

centromeric_repeat_2_mean_pairwise_similarity.perc <- NA
if(centromeric_families.no > 1) {
  sample_repeats <- sample(1 : nrow(repeats_2), max_repeats_to_compare, replace = TRUE)
  similarity_2_matrix <- adist(repeats_1$sequence[sample_repeats])
  centromeric_repeat_2_mean_pairwise_similarity.perc <- 100 * (1 - mean(similarity_2_matrix) / centromeric_repeat_2_mean_size.bp)
}


centromeric_repeat_2_mean_pairwise_similarity_SD.perc_point <- NA
if(centromeric_families.no > 1) {
  centromeric_repeat_2_mean_pairwise_similarity_SD.perc_point <- sd(similarity_2_matrix)
}



centromeric_repeat_2_similarity_between_chromosomes_mean_pairwise.perc <- NA
if(centromeric_families.no > 1) {
  centromeric_repeat_2_similarity_between_chromosomes_mean_pairwise.perc <- mean(within_between_table_2$similarity[within_between_table_2$method == "diff_chr"])
}

centromeric_repeat_2_similarity_between_chromosomes_mean_pairwise_SD.perc_point <- NA
if(centromeric_families.no > 1) {
  centromeric_repeat_2_similarity_between_chromosomes_mean_pairwise_SD.perc_point <- sd(within_between_table_2$similarity[within_between_table_2$method == "diff_chr"])
}

centromeric_repeat_2_similarity_within_chromosomes_mean_pairwise.perc <- NA
if(centromeric_families.no > 1) {
  centromeric_repeat_2_similarity_within_chromosomes_mean_pairwise.perc <- mean(within_between_table_2$similarity[within_between_table_2$method == "same_chr"])
}

centromeric_repeat_2_similarity_within_chromosomes_mean_pairwise_SD.perc_point <- NA
if(centromeric_families.no > 1) {
  centromeric_repeat_2_similarity_within_chromosomes_mean_pairwise_SD.perc_point <- sd(within_between_table_2$similarity[within_between_table_2$method == "same_chr"])
}
print("17")

centromeric_repeat_2_GC.perc <- NA
if(do_GC) {
  if(centromeric_families.no > 1) {
    if(nrow(repeats_2) != 0) {
      centromeric_repeat_2_GC.perc <- 100 * GC(strsplit(paste(repeats_2$sequence, collapse = ""), split = "")[[1]])
    }
    
  }
}

print("18")

gaps_in_centromeric_arrays.no <- NA
if(centromeric_families.no > 0) {
  if(nrow(gaps) != 0) {
    gaps$width <- gaps$end - gaps$start
    gaps <- gaps[gaps$width >= filter_short_gaps,]
    gaps$in_centromere <- FALSE
    gaps$on_edge_array <- FALSE
    
    if(length(centromeric_arrays$arrayID[centromeric_arrays$in_centromere]) != 0) {
      gaps$in_centromere[gaps$repeat_array_id %in% centromeric_arrays$arrayID[centromeric_arrays$in_centromere]] <- TRUE
      
      gaps_in_centromeric_arrays.no <- nrow(gaps[gaps$in_centromere,])
    }
  } else {
    gaps_in_centromeric_arrays.no <- 0
  }
  
}

gaps_in_centromeric_arrays.total_bp <- NA
if(centromeric_families.no > 0) {
  if(nrow(gaps) == 0) {
    gaps_in_centromeric_arrays.total_bp<- 0
  } else {
    gaps_in_centromeric_arrays.total_bp <- sum(gaps$width[gaps$in_centromere])
  }
  
  
}

gaps_in_centromeric_arrays_size_against_arrays.perc <- NA
if(centromeric_families.no > 0) {
  if(nrow(gaps) == 0) {
    gaps_in_centromeric_arrays_size_against_arrays.perc<- 0
  } else {
    gaps_in_centromeric_arrays_size_against_arrays.perc <- 100 * gaps_in_centromeric_arrays.total_bp / sum(centromeric_arrays$width[centromeric_arrays$in_centromere])
  }
  
  
}

gaps_in_centromeric_arrays_size_against_arrays_no_edge_arrays.perc <- NA
if(centromeric_families.no > 0) {
  if(nrow(gaps) == 0) {
    gaps_in_centromeric_arrays_size_against_arrays_no_edge_arrays.perc <- 0
  } else {
    if("edge_filter" %in% names(centromeric_arrays)) {
      gaps_in_centromeric_arrays_size_against_arrays_no_edge_arrays.perc <- 
        100 * sum(gaps$width[gaps$repeat_array_id %in% centromeric_arrays$arrayID[centromeric_arrays$in_centromere & !centromeric_arrays$edge_filter]]) / 
        sum(centromeric_arrays$width[centromeric_arrays$in_centromere & !centromeric_arrays$edge_filter])
    } else {
      gaps_in_centromeric_arrays_size_against_arrays_no_edge_arrays.perc <- NA
    }
  }
  
  
  
}


gaps_in_centromeric_arrays.mean_bp <- NA
if(centromeric_families.no > 0) {
  if(nrow(gaps) == 0) {
    gaps_in_centromeric_arrays.mean_bp <- 0
  } else {
    gaps_in_centromeric_arrays.mean_bp <- mean(gaps$width[gaps$in_centromere])
  }
  
}

gaps_in_centromeric_arrays.median_bp <- NA
if(centromeric_families.no > 0) {
  if(nrow(gaps) == 0) {
    gaps_in_centromeric_arrays.median_bp<- 0
  } else {
    gaps_in_centromeric_arrays.median_bp <- median(gaps$width[gaps$in_centromere])
  }
  
}

gaps_in_centromeric_arrays_TE_coverage.perc <- NA
if(centromeric_families.no > 0) {
  if(nrow(gaps) != 0) {
    overlap_bp <- 0
    
    for(j in chromosomes) {
      edta_temp <- edta[edta$V2 == j,]
      if(nrow(edta_temp) == 0) next
      
      overlap_bp <- overlap_bp + sum(width(reduce(overlapsRanges(IRanges(gaps$start[gaps$in_centromere & gaps$chromosome == j], gaps$end[gaps$in_centromere & gaps$chromosome == j]),
                                                                 IRanges(edta_temp$V5[edta_temp$V2 == j], edta_temp$V6[edta_temp$V2 == j]))  )))
    }
    
    gaps_in_centromeric_arrays_TE_coverage.perc <- 100 * overlap_bp / sum(gaps$width[gaps$in_centromere])
  } else {
    gaps_in_centromeric_arrays_TE_coverage.perc <- 0
  }
  
}

gaps_in_centromeric_arrays_LTR_coverage.perc <- NA

if(centromeric_families.no > 0) {
  if(nrow(gaps) != 0) {
    overlap_bp <- 0
    for(j in chromosomes) {
      edta_temp <- edta[edta$V2 == j & grepl(LTR_names, edta$V4, perl = T, ignore.case = T),]
      if(nrow(edta_temp) == 0) next
      
      overlap_bp <- overlap_bp + sum(width(reduce(overlapsRanges(IRanges(gaps$start[gaps$in_centromere & gaps$chromosome == j], gaps$end[gaps$in_centromere & gaps$chromosome == j]),
                                                                 IRanges(edta_temp$V5[edta_temp$V2 == j], edta_temp$V6[edta_temp$V2 == j]))  )))
    }
    
    gaps_in_centromeric_arrays_LTR_coverage.perc <- 100 * overlap_bp / sum(gaps$width[gaps$in_centromere & !gaps$on_edge_array])
    
  } else {
    gaps_in_centromeric_arrays_LTR_coverage.perc <- 0
  }
  
  
}


print("19")
gaps_in_centromeric_arrays_no_edge_arrays.no <- NA
if(centromeric_families.no > 0) {
  if(nrow(gaps) != 0) {
    gaps_in_centromeric_arrays_no_edge_arrays.no <- nrow(gaps[!gaps$on_edge_array,])
  } else {
    gaps_in_centromeric_arrays_no_edge_arrays.no <- 0
  }
  
}


gaps_in_centromeric_arrays_no_edge_arrays.total_bp <- NA
if(centromeric_families.no > 0) {
  if(nrow(gaps) == 0) {
    gaps_in_centromeric_arrays_no_edge_arrays.total_bp <- 0
  } else {
    gaps_in_centromeric_arrays_no_edge_arrays.total_bp <- sum(gaps$width[gaps$in_centromere & !gaps$on_edge_array])
  }
}

gaps_in_centromeric_arrays_no_edge_arrays.mean_bp <- NA
if(centromeric_families.no > 0) {
  
  if(nrow(gaps) == 0) {
    gaps_in_centromeric_arrays_no_edge_arrays.mean_bp<- 0
  } else {
    gaps_in_centromeric_arrays_no_edge_arrays.mean_bp <- mean(gaps$width[gaps$in_centromere])
    
  }
  
  
  
}

gaps_in_centromeric_arrays_no_edge_arrays.median_bp <- NA
if(centromeric_families.no > 0) {
  if(nrow(gaps) == 0) {
    gaps_in_centromeric_arrays_no_edge_arrays.median_bp<- 0
  } else {
    gaps_in_centromeric_arrays_no_edge_arrays.median_bp <- mean(gaps$width[gaps$in_centromere & !gaps$on_edge_array])
  }
  
  
}

print("20")
gaps_in_centromeric_arrays_no_edge_arrays_TE_coverage.perc <- NA

if(centromeric_families.no > 0) {
  if(nrow(gaps) == 0) {
    gaps_in_centromeric_arrays_no_edge_arrays_TE_coverage.perc <- 0
  } else {
    overlap_bp <- 0
    for(j in chromosomes) {
      edta_temp <- edta[edta$V2 == j,]
      if(nrow(edta_temp) == 0) next
      
      overlap_bp <- overlap_bp + sum(width(reduce(overlapsRanges(IRanges(gaps$start[gaps$in_centromere & !gaps$on_edge_array & gaps$chromosome == j], gaps$end[gaps$in_centromere & !gaps$on_edge_array & gaps$chromosome == j]),
                                                                 IRanges(edta_temp$V5[edta_temp$V2 == j], edta_temp$V6[edta_temp$V2 == j]))  )))
    }
    
    gaps_in_centromeric_arrays_no_edge_arrays_TE_coverage.perc <- 100 * overlap_bp / sum(gaps$width[gaps$in_centromere & !gaps$on_edge_array])
    
  }
  
}


gaps_in_centromeric_arrays_no_edge_arrays_LTR_coverage.perc <- NA

if(centromeric_families.no > 0) {
  if(nrow(gaps) == 0) {
    gaps_in_centromeric_arrays_no_edge_arrays_LTR_coverage.perc<- 0
  } else {
    
    overlap_bp <- 0
    for(j in chromosomes) {
      edta_temp <- edta[edta$V2 == j & grepl(LTR_names, edta$V4, perl = T, ignore.case = T),]
      if(nrow(edta_temp) == 0) next
      
      overlap_bp <- overlap_bp + sum(width(reduce(overlapsRanges(IRanges(gaps$start[gaps$in_centromere & !gaps$on_edge_array & gaps$chromosome == j], gaps$end[gaps$in_centromere & !gaps$on_edge_array & gaps$chromosome == j]),
                                                                 IRanges(edta_temp$V5[edta_temp$V2 == j], edta_temp$V6[edta_temp$V2 == j]))  )))
    }
    
    gaps_in_centromeric_arrays_no_edge_arrays_LTR_coverage.perc <- 100 * overlap_bp / sum(gaps$width[gaps$in_centromere & !gaps$on_edge_array])
    
  }
  
  
}

CENPA_presence.1_yes_0_no <- NA



architecture = genomes_organisation_data$Centromere.architecture[1]

centromeric_repeat_1_mean_HOR_score = NA
centromeric_repeat_1_hors_no = NA
centromeric_repeat_1_hors_mean_block_size = NA
centromeric_repeat_1_hors_mean_block_distance = NA
centromeric_repeat_2_mean_HOR_score = NA
centromeric_repeat_2_hors_no = NA
centromeric_repeat_2_hors_mean_block_size = NA
centromeric_repeat_2_hors_mean_block_distance = NA

grand_table <- data.frame(fasta_name,
                          group = genomes_organisation_data$Group[1],
                          genus,
                          species,
                          architecture,
                          typic = centromeric_stallite_for_species$Typic[1],
                          centromeric_satellites.1_yes_0_no,
                          holocentricity.1_yes_0_no,
                          sequences.no,
                          chromosomes.no,
                          chromosomes_with_centromeric_repeats.no,
                          size_genome.bp,              
                          size_chromosomes.bp,
                          size_non_chromosomes.bp,
                          chromosomes_mean_size.bp,
                          chromosomes_mean_size_SD.bp,
                          size_centromeres.bp,
                          size_pericentromeres.bp,
                          centromere_mean_size.bp,
                          centromere_mean_size_SD.bp,
                          GC_genome.perc = 0,
                          GC_chromosomes.perc = 0,
                          GC_non_chromosomes.perc = 0,
                          GC_centromeres.perc = 0,
                          GC_pericentromeres.perc = 0,
                          GC_arms.perc = 0,
                          GC_exons_total.perc = 0,
                          GC_TE_total.perc = 0,
                          GC_repeats_total.perc = 0,
                          GC_repeats_centromeric.perc = 0,
                          GC_repeats_noncentromeric.perc = 0,
                          genes_chr.no,
                          genes_chr.bp,
                          genes_chr.perc,
                          exons_chr.no,
                          exons_chr.bp,
                          exons_chr.perc,
                          exons_centromeres.no,
                          exons_centromeres.bp,
                          exons_centromeres.perc,
                          exons_pericentromeres.no,
                          exons_pericentromeres.bp,
                          exons_pericentromeres.perc,
                          exons_arms.no,
                          exons_arms.bp,
                          exons_arms.perc,
                          TE_total.no,
                          TE_total.bp,
                          TE_total.perc,
                          TE_C1_LTR.no,
                          TE_C1_LTR.bp,
                          TE_C1_LTR.perc,
                          TE_C1_nonLTR.no,
                          TE_C1_nonLTR.bp,
                          TE_C1_nonLTR.perc,
                          TE_C2_TIR.no,
                          TE_C2_TIR.bp,
                          TE_C2_TIR.perc,
                          TE_C2_nonTIR.no,
                          TE_C2_nonTIR.bp,
                          TE_C2_nonTIR.perc,
                          TE_LINE.no,
                          TE_LINE.bp,
                          TE_LINE.perc,
                          TE_HELITRON.no,
                          TE_HELITRON.bp,
                          TE_HELITRON.perc,
                          TE_total_centromere.perc,
                          TE_total_pericentromere.perc,
                          TE_total_arms.perc,
                          TE_COPIA_centromere.perc,
                          TE_COPIA_pericentromere.perc,
                          TE_COPIA_arms.perc,
                          TE_GYPSY_centromere.perc,
                          TE_GYPSY_pericentromere.perc,
                          TE_GYPSY_arms.perc,
                          TE_HELITRON_centromere.perc,
                          TE_HELITRON_pericentromere.perc,
                          TE_HELITRON_arms.perc,
                          TE_C1_LTR_centromere.perc,
                          TE_C1_LTR_pericentromere.perc,
                          TE_C1_LTR_arms.perc,
                          TE_C1_nonLTR_centromere.perc,
                          TE_C1_nonLTR_pericentromere.perc,
                          TE_C1_nonLTR_arms.perc,
                          TE_C2_TIR_centromere.perc,
                          TE_C2_TIR_pericentromere.perc,
                          TE_C2_TIR_arms.perc,
                          TE_C2_nonTIR_centromere.perc,
                          TE_C2_nonTIR_pericentromere.perc,
                          TE_C2_nonTIR_arms.perc,
                          repeats_total.no,
                          repeats_total.bp,
                          repeats_total.perc,
                          centromeric_families.no,
                          repeats_centromeric.no,
                          repeats_centromeric.bp,
                          repeats_centromeric.perc,
                          repeats_noncentromeric.no,
                          repeats_noncentromeric.bp,
                          repeats_noncentromeric.perc,
                          repeats_centromeric_on_chromosomes.bp,
                          repeats_centromeric_not_on_chromosomes.bp,
                          repeats_centromeric_not_on_chromosomes_vs_all_centromeric.perc,
                          centromeric_repeats_fraction_on_dominant_strand.perc,
                          centromeric_repeats_in_centromere.perc,
                          centromeric_repeat_1_name,
                          centromeric_repeat_1_TRASH_identifier = centromeric_stallite_for_species$TRASH_name_dec2024runs[1],
                          centromeric_repeat_1_consensus_unambiguous,
                          centromeric_repeat_1_consensus_ambiguous,
                          centromeric_repeat_1_mean_size.bp,
                          centromeric_repeat_1_mean_size_SD.bp,
                          centromeric_repeat_1_mean_size_SD.perc,
                          centromeric_repeat_1_mean_pairwise_similarity.perc,
                          centromeric_repeat_1_mean_pairwise_similarity_SD.perc,
                          centromeric_repeat_1_similarity_between_chromosomes_mean_pairwise.perc,
                          centromeric_repeat_1_similarity_between_chromosomes_mean_pairwise_SD.perc_point,
                          centromeric_repeat_1_similarity_within_chromosomes_mean_pairwise.perc,
                          centromeric_repeat_1_similarity_within_chromosomes_mean_pairwise_SD.perc_point,
                          centromeric_repeat_1_GC.perc,
                          centromeric_repeat_1_mean_HOR_score,
                          centromeric_repeat_1_hors_no,
                          centromeric_repeat_1_hors_mean_block_size,
                          centromeric_repeat_1_hors_mean_block_distance,
                          centromeric_repeat_2_name,
                          centromeric_repeat_2_TRASH_identifier = centromeric_stallite_for_species$TRASH_name_dec2024runs[2],
                          centromeric_repeat_2_consensus_unambiguous,
                          centromeric_repeat_2_consensus_ambiguous,
                          centromeric_repeat_2_mean_size.bp,
                          centromeric_repeat_2_mean_size_SD.bp,
                          centromeric_repeat_2_mean_size_SD.perc,
                          centromeric_repeat_2_mean_pairwise_similarity.perc,
                          centromeric_repeat_2_mean_pairwise_similarity_SD.perc_point,
                          centromeric_repeat_2_similarity_between_chromosomes_mean_pairwise.perc,
                          centromeric_repeat_2_similarity_between_chromosomes_mean_pairwise_SD.perc_point,
                          centromeric_repeat_2_similarity_within_chromosomes_mean_pairwise.perc,
                          centromeric_repeat_2_similarity_within_chromosomes_mean_pairwise_SD.perc_point,
                          centromeric_repeat_2_GC.perc,
                          centromeric_repeat_2_mean_HOR_score,
                          centromeric_repeat_2_hors_no,
                          centromeric_repeat_2_hors_mean_block_size,
                          centromeric_repeat_2_hors_mean_block_distance,
                          gaps_in_centromeric_arrays.no,
                          gaps_in_centromeric_arrays.total_bp,
                          gaps_in_centromeric_arrays_size_against_arrays.perc,
                          gaps_in_centromeric_arrays_size_against_arrays_no_edge_arrays.perc,
                          gaps_in_centromeric_arrays.mean_bp,
                          gaps_in_centromeric_arrays.median_bp,
                          gaps_in_centromeric_arrays_TE_coverage.perc,
                          gaps_in_centromeric_arrays_LTR_coverage.perc,
                          gaps_in_centromeric_arrays_no_edge_arrays.no,
                          gaps_in_centromeric_arrays_no_edge_arrays.total_bp,
                          gaps_in_centromeric_arrays_no_edge_arrays.mean_bp,
                          gaps_in_centromeric_arrays_no_edge_arrays.median_bp,
                          gaps_in_centromeric_arrays_no_edge_arrays_TE_coverage.perc,
                          gaps_in_centromeric_arrays_no_edge_arrays_LTR_coverage.perc,
                          CENPA_presence.1_yes_0_no
)

gc_exists <- FALSE
# add GC
if(file.exists(paste0("./", fasta_name, "_", genomes_organisation_data$Centromere.architecture[1], "_GC_grand_table.csv"))) {
  gc <- read.csv(paste0(fasta_name, "_", genomes_organisation_data$Centromere.architecture[1], "_GC_grand_table.csv"))
  grand_table$GC_genome.perc <- gc$GC_genome.perc
  grand_table$GC_chromosomes.perc <- gc$GC_chromosomes.perc
  grand_table$GC_non_chromosomes.perc <- gc$GC_non_chromosomes.perc
  grand_table$GC_centromeres.perc <- gc$GC_centromeres.perc
  grand_table$GC_pericentromeres.perc <- gc$GC_pericentromeres.perc
  grand_table$GC_arms.perc <- gc$GC_arms.perc
  grand_table$GC_exons_total.perc <- gc$GC_exons_total.perc
  grand_table$GC_TE_total.perc <- gc$GC_TE_total.perc
  grand_table$GC_repeats_total.perc <- gc$GC_repeats_total.perc
  grand_table$GC_repeats_centromeric.perc <- gc$GC_repeats_centromeric.perc
  grand_table$GC_repeats_noncentromeric.perc <- gc$GC_repeats_noncentromeric.perc
  gc_exists <- TRUE
}


hor_exists <- FALSE
# add HORs for repeat 1  and 2
if(file.exists(paste0("./", fasta_name, "_", genomes_organisation_data$Centromere.architecture[1], "_", centromeric_stallite_for_species$Satellite_name_current[1], "_HORs_grand_table.csv"))) {
  hors_repeat_1 <- read.csv(paste0(fasta_name, "_", genomes_organisation_data$Centromere.architecture[1], "_", centromeric_stallite_for_species$Satellite_name_current[1], "_HORs_grand_table.csv"))
  
  grand_table$centromeric_repeat_1_mean_HOR_score <- hors_repeat_1$repeats_mean_HOR_score
  grand_table$centromeric_repeat_1_hors_no <- hors_repeat_1$hors_no
  grand_table$centromeric_repeat_1_hors_mean_block_size <- hors_repeat_1$hors_mean_block_size
  grand_table$centromeric_repeat_1_hors_mean_block_distance <- hors_repeat_1$hors_mean_block_distance
  
  hor_exists <- TRUE
  
}

if(file.exists(paste0("./", fasta_name, "_", genomes_organisation_data$Centromere.architecture[2], "_", centromeric_stallite_for_species$Satellite_name_current[2], "_HORs_grand_table.csv"))) {
  hors_repeat_2 <- read.csv(paste0(fasta_name, "_", genomes_organisation_data$Centromere.architecture[2], "_", centromeric_stallite_for_species$Satellite_name_current[2], "_HORs_grand_table.csv"))
  
  grand_table$centromeric_repeat_2_mean_HOR_score <- hors_repeat_2$repeats_mean_HOR_score
  grand_table$centromeric_repeat_2_hors_no <- hors_repeat_2$hors_no
  grand_table$centromeric_repeat_2_hors_mean_block_size <- hors_repeat_2$hors_mean_block_size
  grand_table$centromeric_repeat_2_hors_mean_block_distance <- hors_repeat_2$hors_mean_block_distance
}

if(centromeric_satellites.1_yes_0_no == 0) {
  hor_exists <- TRUE
}



if(gc_exists & hor_exists) {
  write.csv(x = grand_table, file = paste0("FULL_", fasta_name, "_", genomes_organisation_data$Centromere.architecture[1], "_grand_table.csv"), row.names = FALSE)
  write.csv(x = grand_table, file = paste0("/home/pwlodzimierz/ToL/upload_files/grand_tables/FULL_", fasta_name, "_", genomes_organisation_data$Centromere.architecture[1],"_grand_table.csv"), row.names = FALSE)
  
} else if(gc_exists) {
  write.csv(x = grand_table, file = paste0("No_HOR_", fasta_name, "_", genomes_organisation_data$Centromere.architecture[1], "_grand_table.csv"), row.names = FALSE)
  write.csv(x = grand_table, file = paste0("/home/pwlodzimierz/ToL/upload_files/grand_tables/No_HOR_", fasta_name, "_", genomes_organisation_data$Centromere.architecture[1],"_grand_table.csv"), row.names = FALSE)
  
} else if(hor_exists) {
  write.csv(x = grand_table, file = paste0("No_GC_", fasta_name, "_", genomes_organisation_data$Centromere.architecture[1], "_grand_table.csv"), row.names = FALSE)
  write.csv(x = grand_table, file = paste0("/home/pwlodzimierz/ToL/upload_files/grand_tables/No_GC_", fasta_name, "_", genomes_organisation_data$Centromere.architecture[1],"_grand_table.csv"), row.names = FALSE)
  
} else {
  write.csv(x = grand_table, file = paste0("No_HOR_No_GC_", fasta_name, "_", genomes_organisation_data$Centromere.architecture[1], "_grand_table.csv"), row.names = FALSE)
  write.csv(x = grand_table, file = paste0("/home/pwlodzimierz/ToL/upload_files/grand_tables/No_HOR_No_GC_", fasta_name, "_", genomes_organisation_data$Centromere.architecture[1],"_grand_table.csv"), row.names = FALSE)
  
}







