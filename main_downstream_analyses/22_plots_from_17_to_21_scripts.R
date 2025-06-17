
# for each genome, make a set of plots:
# a. repeat sequence based PCA for a 10 000 repeat sample from each genome, coloured by chromosome
# --> 20_sample_ali_...R two-column PCA df csv
# b. repeat sequence based tree for a 10 000 repeat sample from each genome, coloured by chromosome
# --> 20_sample_ali_...R newick format tree
# c. sample 10 000 pairs of repeats, calculate their similarity, plot simialrity vs physical distance (within chromosomes)
#!/usr/bin/env Rscript
.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))

# --> 21_repeat_similarity_prox_far_inter_intra_chr.R two-column df csv 1
# d. sample 10 000 repeats, calculate the similarity with:
#     10 000 repeats in a 100-repeat distance 
#     10 000 repeats in a > 1000-repeat distance 
#     10 000 repeats on different chromosomes
#    and make a 3-bar bar-plot for those (try to keep the y axis on the same range as the c plot for easier comparison)
# --> 21_repeat_similarity_prox_far_inter_intra_chr.R two-column df csv 2 3 and 4
# e. sequence length conservation: sample 10 000 repeat pairs, align them, check how many insertions, deletions and substitutions there are from
#    the perspective of one of the repeats. Make bar plots for each insertion:deletion combination with y axis being substitutions number: 0-0, 0-1, 
#    1-0, 1-1, ect. The idea is to show that it's more common to have the 1-1, 2-2, 3-3 etc. arrangement, that preserve the repeat length, let's call
#    it "compensating mutations"
# --> 19_sequence_...R occurances of each mutation csv
# f. HOR score in a 100 kbp distance from gaps. For each repeat, check which TEs are in 100 kbp distance, for each case save into a df the two
#    parameters: HOR score and physical distance
# --> 18_hor_score...R two-column df csv
# g. make a matrix of gap similarity scores, make a tree with gap lengths and TE annotation ploted alongside
# --> 17_TE_...R matrix csv 1
# h. add to the values calculated in g: merge up and downstream repeats and make matrix of those repeat regions similarity (direct alignment). 
#    Plot a scatter of gap (g) vs repeat similarity scores. 
# --> 17_TE_...R matrix csv 2
# i. measure similarity between 1 000 randomly picked repeat region pairs using the method as in g. plot the histogram of those similarities, marking 
#    where the values from (g) fall along the distribution
# --> 17_TE_...R two-column df csv
# 
# ----> 22_...R this file plots
taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 193
print(i)

suppressMessages(library(seqinr))




### load data
data_directories <- list.dirs(path = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs", recursive = FALSE, full.names = TRUE)

setwd(data_directories[i])
repeats_file <- list.files(path = ".", pattern = "repeats_filtered.csv", full.names = TRUE)
if(length(repeats_file) != 1) {cat("issue with repeats files: \n", repeats_file); stop()}

edta_file <- list.files(path = ".", pattern = "_edta_filtered.csv.reassigned", full.names = TRUE)
if(length(edta_file) != 1) {cat("issue with edta files: \n", edta_file); stop()}

arrays_file <- list.files(path = ".", pattern = "_centromeric_arrays.csv", full.names = TRUE)
if(length(arrays_file) != 1) {cat("issue with gaps files: \n", arrays_file); stop()}

gaps_file <- list.files(path = ".", pattern = "_centromeric_gaps.csv", full.names = TRUE)
if(length(gaps_file) != 1) {cat("issue with gaps files: \n", gaps_file); stop()}

repeats <- read.csv(repeats_file)
edta <- read.csv(edta_file)
arrays <- read.csv(arrays_file)
gaps <- read.csv(gaps_file)

satellites_metadata <- read.csv(file = "/home/pwlodzimierz/ToL/curated_satellites_repDec24_jan2025.csv")
chr_sizes_metadata <- read.csv(file = "/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full.Ian.csv")

### filter the data

assembly_name_long <- strsplit(data_directories[i], split = "/")[[1]][7]
assembly_name_short <- strsplit(assembly_name_long, split = ".fa")[[1]][1]

satellites_metadata <- satellites_metadata[satellites_metadata$Genome == assembly_name_long, ]

satellite_TRASH_names <- strsplit(satellites_metadata$TRASH_name_dec2024runs[1], split = ";")[[1]]

chr_sizes_metadata <- chr_sizes_metadata[chr_sizes_metadata$assembly.name == assembly_name_long, ]

repeats <- repeats[repeats$new_class %in% satellite_TRASH_names, ]

###








