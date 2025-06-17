#!/bin/bash

###
#
# This will also rescore the scores of 6_find_centromeric_repeats.R
#
###

#
# There is a hardcoded limit to 200 plot rows, so gotta limit number of plotted chrs to 30...
#


### EDTA unique classes in the new V3 edta filtered cls tables:
edta_classes <- list(
  # class I (retrotransposons)
  ## LTR retrotransposons
  c("Gypsy_LTR_retrotransposon"),
  c("Copia_LTR_retrotransposon"),
  c("Bel_Pao_LTR_retrotransposon"),
  c("TRIM_LTR_retrotransposon"),
  c("Caulimoviridae"),
  c("Retrovirus", "LTR_retrotransposon", "long_terminal_repeat"),
  ## Non-LTR retrotransposons
  c("LINE_element"),
  c("SINE_element"),
  c("Penelope_retrotransposon"),
  c("DIRS_YR_retrotransposon"),
  c("non_LTR_retrotransposon"),
  # class II (DNA transposons)
  ## TIRs
  c("Kolobok_TIR_transposon" , "Ginger_TIR_transposon", "Academ_TIR_transposon", "Novosib_TIR_transposon", "Sola_TIR_transposon", "Merlin_TIR_transposon", "IS3EU_TIR_transposon", "PiggyBac_TIR_transposon", "hAT_TIR_transposon", "Mutator_TIR_transposon", "Tc1_Mariner_TIR_transposon", "Dada_TIR_transposon", "CACTA_TIR_transposon", "Zisupton_TIR_transposon", "PIF_Harbinger_TIR_transposon"),
  ## other class II
  c("DNA_transposon"),
  c("helitron"),
  c("MITE"),
  c("Maverick_Polinton", "polinton"),
  # other, recombinase element based
  c("Tyrosine_Recombinase_Elements", "Crypton_Tyrosine_Recombinase"),
  # others
  c("TE", "TE_unclass"),
  # likely not TEs, remove for plotting?
  c("repeat_region", "SUPER", "Sequence_Ontology", "rRNA_gene", "target_site_duplication", "chr"))
edta_classes_colours <-  c(
  "#E31A1C",  # red
  "#D55E00",  # reddish-orange
  "#F5793A",  # bright orange
  "#FF7F00",  # orange
  "#FDBF6F",  # peach
  "#F0E442",  # yellow
  "#6A3D9A",  # dark purple
  "#A95AA1",  # purple
  "#CC79A7",  # pink
  "#DDA0DD",  # light pinkish purple (plum)
  "#CAB2D6",  # lavender
  "#1F78B4",  # blue
  "#B2DF8A",  # light green
  "#33A02C",  # green
  "#009E73",  # teal green
  "#56B4E9",  # light blue
  "#7F7F7F",  # grey
  "#000000",   # black
  "#000000"   # black
)


.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))

args = commandArgs(trailingOnly = TRUE)
print(args)


library(seqinr)
library(stringr)
# library(gridExtra)
library(Biostrings)
library(GenomicRanges)
library(scales)
mafft.bat.file = "/home/pwlodzimierz/TRASH/src/mafft-linux64/mafft.bat"
source(file = "/home/pwlodzimierz/useful.scripts.R")

setwd("/home/pwlodzimierz/ToL/git_ToL")
source("./aux_fun.R")
ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}
replace_existing_analysis = TRUE

add_centromeric_files_to_plots <- F


data_directories <- list.dirs(path = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs", recursive = FALSE, full.names = TRUE)
data_directories <- data_directories[!grepl(pattern = "templated_", data_directories)]
data_directories <- data_directories[grepl(pattern = ".fa", data_directories)]
assembly_files <- list.files(path = "/home/pwlodzimierz/ToL/Assemblies/fastas_2021_Michael", recursive = FALSE, full.names = TRUE)
assembly_files <- assembly_files[!grepl(".fai", assembly_files)]


taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 15
print(i)

print(paste0(i, " / ", length(data_directories)))
### Load data ================================================================
setwd(data_directories[i])
print(getwd())

assembly_name = strsplit(strsplit(data_directories[i], split = ".fa")[[1]][1], split = "v2_out_for_HORs/")[[1]][2]
assembly_file = grep(assembly_name, assembly_files)
print(assembly_file)
print(assembly_files[assembly_file])



if(!replace_existing_analysis) {
  if(file.exists(paste0("global_plot_plus_edta_", assembly_name, ".png"))) {
    if(file.size(paste0("global_plot_plus_edta_", assembly_name, ".png")) != 0) {
      print("Analysis finished")
      quit(save = "no", status = 1)
    }
  }
} else {
  if(file.exists(paste0("global_plot_plus_edta_", assembly_name, ".png"))) {
    file.remove(paste0("global_plot_plus_edta_", assembly_name, ".png"))
  }
}

if(add_centromeric_files_to_plots) {
  centromere_files <- list.files(path = "/home/pwlodzimierz/ToL/final_summary_analyses/chromosome_data_csvs", full.names = TRUE)
  centromere_file <- centromere_files[grep(assembly_name, centromere_files)]
  if(length(centromere_file) != 1) {warning(paste0(i, "No centromere_file!")); setwd(".."); quit(save = "no", status = 1)}
  centromere_data = read.csv(file = centromere_file, fileEncoding = "UTF-8")
}



helixer_file = list.files(pattern = "helixer_filtered.gff", full.names = TRUE)
if(length(helixer_file) != 1) {warning(paste0(i, "No genes!")); setwd(".."); quit(save = "no", status = 1)}

repeat_file = list.files(pattern = "_repeats_filtered.csv", full.names = TRUE)
if(length(repeat_file) != 1) {warning(paste0(i, "No repeats!")); setwd(".."); quit(save = "no", status = 1)}

classes_file = list.files(pattern = "_classes_merged_filtered", full.names = TRUE)
if(length(classes_file) != 1) {warning(paste0(i, " no classes!")); setwd(".."); quit(save = "no", status = 1)}

array_file = list.files(pattern = "_arrays_filtered.csv", full.names = TRUE)
if(length(array_file) != 1) {warning(paste0(i, " no arrays!")); setwd(".."); quit(save = "no", status = 1)}

edta_file = list.files(pattern = paste0(assembly_name, "_edta_filtered.csv"), full.names = TRUE)
if(length(edta_file) != 1) {warning(paste0(i, " no edta!")); setwd(".."); quit(save = "no", status = 1)}



print("load annotations")
repeats = read.csv(file = repeat_file)
arrays = read.csv(file = array_file)
classes = read.csv(file = classes_file)
classes$num_ID <- 1 : nrow(classes)
edta = read.csv(file = edta_file)
genes <- read.table(file = helixer_file, header = FALSE, sep = "\t", skip = 4)
genes <- genes[genes$V3 == "CDS", ]
scores <- read.csv(file = paste0("genome_classes_", assembly_name, ".csv"))
scores <- scores[order(scores$chromosome, decreasing = FALSE), ]
assembly.seq = read.fasta(assembly_files[assembly_file], forceDNAtolower = T, seqtype = "DNA")

edta_cls <- read.csv(file = paste0("/home/pwlodzimierz/ToL/mnt/RAID/data/estela_temp/getfasta_d/solo_ltr/ltr_fine_tuning/map_app/v2/edta_filtered_cls/", 
                                   assembly_name, "_edta_filtered_cls.csv"))


print("annotations loaded")


### filter data to include only chromosomes ====================================

print("filter chromosomes")

all_chromosomes <- read.csv("/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full.Ian.csv")
all_chromosomes <- all_chromosomes[all_chromosomes$assembly.name == strsplit(data_directories[i], split = "v2_out_for_HORs/")[[1]][2], ]
all_chromosomes <- all_chromosomes[all_chromosomes$is.chr == 1, ]
chromosomes <- all_chromosomes$chromosome.name
chromosomes_lengths <- all_chromosomes$size

# ==============================================================================


# filter top 30 longest chromosomes ============
if(length(chromosomes_lengths) > 30) {
  chromosomes <- chromosomes[chromosomes_lengths %in% sort(chromosomes_lengths, decreasing = TRUE)[1:30]]
  chromosomes_lengths <- chromosomes_lengths[chromosomes_lengths %in% sort(chromosomes_lengths, decreasing = TRUE)[1:30]]
  
}
# ==============================================

# rescore the classes in "scores" df =======================================
### copy over to 1_visualise_scores on the laptop 
for(j in seq_along(chromosomes)) {
  
  scores_chr <- scores[scores$chromosome == chromosomes[j], ]
  
  if(nrow(scores_chr) == 0) next
  
  {
    # soft normalise by total bp to avoid very short arrays scoring high just because they are in a good place
    size_normalisation <- scores_chr$total_bp / max(scores_chr$total_bp) # multiplier x0 to x1 as a fraction of the biggest repeat class
    size_normalisation <- size_normalisation + (1 - size_normalisation) / 1.5 # multiplier x0.(6) to x1 as a fraction of the biggest repeat class
    # 1+0/2=1   0.5+0.5/2=0.75    0+1/2=0.5         0.5:1
    # 1+0/1.5=1   0.5+0.5/1.5=0.83    0+1/1.5=0.66  0.(6):1
    scores_chr$score_total <- 0
    
    # total_bp_norm_chr, 2 at 10%
    values <- scores_chr$total_bp_norm_chr
    values[values > 0.1] = 0.1
    scores_chr$score_total <- scores_chr$score_total + (values / 0.1) * 2
    
    # total_bp_norm_rep, 1 at 50%
    values <- scores_chr$total_bp_norm_rep
    values[values > 0.5] = 0.5
    scores_chr$score_total <- scores_chr$score_total + (values / 0.5) * 1
    
    # start_sd_norm_chr, 0 at 15%, 1 at 0%
    values <- scores_chr$start_sd_norm_chr
    values[values > 0.15] = 0.15
    scores_chr$score_total <- scores_chr$score_total + (1 - (values / 0.15)) * 1 * size_normalisation
    
    # start_norm_chr_0_50, 1 at 0%, 0 at 25% and 1 at 50%
    values <- scores_chr$start_norm_chr_0_50
    values[values > 0.25] = abs(values[values > 0.25] - 0.5)
    scores_chr$score_total <- scores_chr$score_total + (1-(values / 0.25)) * 1 * size_normalisation
    
    # gaps_with_TEs_fraction, 0.5 at 75%
    values <- scores_chr$gaps_with_TEs_fraction
    values[values == -1] = 0
    values[values > 0.75] = 0.75
    scores_chr$score_total <- scores_chr$score_total + (values / 0.75) * 0.5 * size_normalisation
    
    # centre_array_edit, 0 at 20% 1 at 0%
    values <- scores_chr$centre_array_edit
    values[values == -1] = 20
    values[values > 20] = 20
    scores_chr$score_total <- scores_chr$score_total + (1 - (values / 20)) * 2 * size_normalisation
    
    # centre_array_width_sd, 0 at 20 1 at 0%
    values <- scores_chr$centre_array_width_sd
    values[values == -1] = 20
    values[values > 20] = 20
    scores_chr$score_total <- scores_chr$score_total + (1 - (values / 20)) * 2 * size_normalisation
    
    # centre_chromosome_edit, 0 at 15% 0.5 at 0%
    values <- scores_chr$centre_chromosome_edit
    values[values == -1] = 15
    values[values > 15] = 15
    scores_chr$score_total <- scores_chr$score_total + (1 - (values / 15)) * size_normalisation
    
    # centre_chromosome_width_sd, 0.5 at 15%
    values <- scores_chr$centre_chromosome_width_sd
    values[values == -1] = 15
    values[values > 15] = 15
    scores_chr$score_total <- scores_chr$score_total + (1 - (values / 15)) * size_normalisation
    
    # rescore the TE scores to not reach super big values by adding 1 to the denominator
    scores_chr$TE_prox_score <- abs(scores_chr$TE_lm_coef) / (1 + scores_chr$TE_prox_dist + scores_chr$TE_prox_SD) * 100
    
    # TE_prox_score, 5 at 5
    values <- scores_chr$TE_prox_score
    values[values == -1] = 0
    values[values > 5] = 5
    scores_chr$score_total <- scores_chr$score_total + (values / 100) * size_normalisation * 5
    
    # rescore the TE scores to not reach super big values by adding 1 to the denominator
    scores_chr$gene_prox_score <- abs(scores_chr$gene_lm_coef) / (1 + scores_chr$gene_prox_dist + scores_chr$gene_prox_SD) * 100
    
    # gene_prox_score, 5 at 5
    values <- scores_chr$gene_prox_score
    values[values == -1] = 0
    values[values > 5] = 5
    scores_chr$score_total <- scores_chr$score_total + (values / 100) * size_normalisation * 5
    
    # # t_test_p_val, 0 at 0.001; t_test_t_val must be negative
    # values <- scores_chr$t_test_p_val
    # values[values == -1] = 0.001
    # values[values > 0.001] = 0.001
    # values[scores_chr$t_test_t_val > 0] = 0.001
    # scores_chr$score_total <- scores_chr$score_total + (1 - (values / 0.001)) * size_normalisation
  }
  scores$score_total[scores$chromosome == chromosomes[j]] <- scores_chr$score_total
}

scores <- scores[, -1]





# ===========================================================

# choose classes from "scores" df to plot ==============

# 1 is one of 3 top scores of at least 4 on at least one chromosome
# 2 at least 5 kbp in total size
# not more than 12 classes total

scores_filtered <- NULL

# get all classes that are not top 3 or scored 5 or 5 kbp total out
for(j in seq_along(chromosomes)) {
  scores_chromosome <- scores[scores$chromosome == chromosomes[j],]
  scores_chromosome <- scores_chromosome[order(scores_chromosome$score_total, decreasing = TRUE), ]
  scores_chromosome <- scores_chromosome[scores_chromosome$score_total >= 4,]
  if(nrow(scores_chromosome) > 3) scores_chromosome <- scores_chromosome[(1:3), ]
  if(nrow(scores_chromosome) == 0) next
  for(k in seq_len(nrow(scores_chromosome))) {
    scores_chromosome$total_genome_bp[k] = classes$sum_coverage[classes$class == scores_chromosome$class[k]][1]
  }
  scores_chromosome <- scores_chromosome[scores_chromosome$total_genome_bp >= 5000,]
  if(nrow(scores_chromosome) == 0) next
  scores_filtered <- rbind(scores_filtered, scores_chromosome)
}

classes$to_plot <- FALSE
classes$to_plot[classes$class %in% scores_filtered$class] <- TRUE
classes <- classes[classes$to_plot,]
if(nrow(classes) > 12) classes <- classes[(1:12), ]

scores$to_plot <- FALSE
scores$to_plot[scores$class %in% classes$class] <- TRUE
scores <- scores[scores$to_plot, ]

# =================================================



# start the plot ====================================

# Top half: chromosomes one under the other:
#  table 1: repeat scores
#  plot A: GC, repeats and repeats families
#  plot B: TE+repeats peak vs gene valley plots

bin.size.plots = 10000
bin.size.plots.for.GC = 2000
bin.size.plots.for.EDTA = 50000
plot.x.ticks = 1000000
min.seq.size.to.plot = plot.x.ticks
safe_colorblind_palette = rep(c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                                "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888"),100)
print(paste0("Genome ", i, ", ", assembly_name, " chromosomes count: ", length(chromosomes)))

png(filename = paste0("global_plot_plus_edta_", assembly_name, ".png"), width = 3000, 
    height = (700 + 780*length(chromosomes) + 3000 + 200 + 80 * nrow(classes)), pointsize = 32)

layout(matrix(1:(4*length(chromosomes)+3), nrow = (3 + 4*length(chromosomes)), ncol = 1), 
       widths = 1, heights = c(700, rep(c(80,300,200,200),length(chromosomes)), 3000, (200 + 80 * nrow(classes))))

par(mar = c(0.5,2,0.2,0.2), mgp = c(0.5, 0.5, 0), oma = c(2, 1, 3, 2))


plot(NA,NA, xlim = c(1,100), ylim = c(1,100), xlab = "", ylab = "", axes = F)
text(x = 50, y = 95, labels = assembly_name, pos = 1, cex = 5)
text(x = 1, y = 75, 
     labels = paste0("Total size: ", formatC(sum(chromosomes_lengths)/1000000, format = "f", big.mark = " ", digits = 3), " Mbp"), 
     pos = 4, cex = 3)
text(x = 40, y = 75, 
     labels = paste0("Chromosomes count: ", formatC(length(chromosomes_lengths), format = "f", big.mark = " ", digits = 0), " "), 
     pos = 4, cex = 3)
text(x = 1, y = 60, 
     labels = "Top table: significant repeat families (up to 12 with highest count plotted per genome)", 
     pos = 4, cex = 2)
text(x = 1, y = 50, 
     labels = "First plot: GC content per 2 kbp (black line); repeat per 10 kbp (grey bars); repeat classes per 10 kbp (colours as in tables)", 
     pos = 4, cex = 2)
text(x = c(1,56,77), y = 40, 
     labels = c("Second plot: TEs per 50 kbp (total: grey bars; classes: colours), ",
                "averaged genes content, ", 
                "averaged TE and repeats content"), 
     col = c("black", "#00bb33", "#0066aa"),
     pos = 4, cex = 2)
text(x = c(1,15,30,45,60,75,90), y = 32, 
     labels =  c("class I LTR: ",
                 "Gypsy",
                 "Copia",
                 "Bel Pao",
                 "TRIM",
                 "Caulimoviridae", 
                 "unspecified"), 
     col = c("black", edta_classes_colours[1:6]),
     pos = 4, cex = 1.5)
text(x = c(1,15,30,45,60,75), y = 26, 
     labels = c("class I non-LTR: ", 
                "LINE",
                "SINE",
                "Penelope",
                "DIRS YR",
                "unspecified"),
     col = c("black", edta_classes_colours[7:11]),
     pos = 4, cex = 1.5)
text(x = c(1,15), y = 20, 
     labels = c("class II TIRs: ", "Kolobok;Ginger;Academ;Novosib;Sola;Merlin;IS3EU;PiggyBac;hAT;Mutator;Tc1 Mariner;Dada;CACTA;Zisupton;PIF Harbinger"),
     col = c("black", edta_classes_colours[12]), 
     pos = 4, cex = 1.5)
text(x = c(1,15,30,45,60), y = 14, 
     labels = c("class II others: ",  
                "DNA_transposon",
                "helitron",
                "MITE",
                "Maverick Polinton"),
     col = c("black", edta_classes_colours[13:16]), 
     pos = 4, cex = 1.5)
text(x = c(1,15,30), y = 8, 
     labels = c("other: ", "Tyrosine Recombinase",
     "unspecified"),
     col = c("black", edta_classes_colours[17:18]), 
     pos = 4, cex = 1.5)



# ===================================================

# Chromosome tracks, each chromosome one under another, GC, repeats, families
for(j in 1 : length(chromosomes))
{
  print(paste0("Genome ", i, ", ", assembly_name, " | Chromosome ", j, "/", length(chromosomes)))
  
  plot(NA,NA, xlim = c(1,100), ylim = c(1,100), xlab = "", ylab = "", axes = F)
  text(x = 1, y = 15, labels = assembly_name, pos = 4, cex = 1.4)
  text(x = 15, y = 15, labels = chromosomes[j], pos = 4, cex = 1.2)
  text(x = 30, y = 15, 
       labels = paste0(formatC(chromosomes_lengths[j]/1000000, format = "f", big.mark = " ", digits = 3), " Mbp"), 
       pos = 4, cex = 1.2)
  abline(h = 50, lwd = 3)
  abline(h = 60, lwd = 3)
  
  cat(" plot started;")
  window.starts = genomic.bins.starts(start = 1, end = chromosomes_lengths[j], bin.size = bin.size.plots)
  
  ### repeat scores table ===============================
  
  scores_chromosome <- scores[scores$chromosome == chromosomes[j], ]
  scores_chromosome_plot <- scores_chromosome[,c(1,2,3,4,12,13,29)]
  scores_chromosome_plot[,c(3,5,6,7)] = round(scores_chromosome_plot[,c(3,5,6,7)], 2)
  
  if(nrow(scores_chromosome) != 0) {
    
    scores_chromosome_plot[scores_chromosome_plot == -1] = "NA"
    
    scores_chromosome_plot$colours <- ""
    for(k in seq_len(nrow(scores_chromosome_plot))) {
      scores_chromosome_plot$colours[k] = safe_colorblind_palette[which(classes$class == scores_chromosome_plot$class[k])]
    }
    
  } 
  
  cat(" table;")
  create_table(scores_chromosome_plot[, c(1:7)], 
               c("Class", "Repeats no", "Mean length", "Total length", "Centre array ED",
                 "Centre array width SD", "Centromeric score"), 
               colours = scores_chromosome_plot$colours)
  
 
  
  # =====================================================
  
  
  ### plot A  ===========================================
  
  cat(" plot A;")
  # plot Repeats
  repeats.chromosome = repeats[repeats$seqID == chromosomes[j],]
  repeat.frequencies = calculate.repeats.percentage.in.windows(windows.starts = window.starts, 
                                                               repeat.starts = repeats.chromosome$start, 
                                                               repeat.lengths = repeats.chromosome$width, 
                                                               sequence.length = chromosomes_lengths[j])
  repeat.frequencies[repeat.frequencies == 0] = NA
  repeat.frequencies[repeat.frequencies > 100] = 100
  plot(window.starts, repeat.frequencies, type = "h", frame.plot = F, xlim = c(1, chromosomes_lengths[j]), col = "#CCCCCC", 
       ylab = "", xlab = "", xaxt = "n", ylim = c(0,100))
  
  if(add_centromeric_files_to_plots) {
    # add centromere background colour
    starts <- as.numeric(strsplit(centromere_data$centromere_arrays_starts[j], split = ";")[[1]])
    ends <- as.numeric(strsplit(centromere_data$centromere_arrays_ends[j], split = ";")[[1]])
    for(as in seq_along(starts)) {
      rect(starts[as], 0, ends[as], 100, col = "#00000030", border = NULL)
    }
    
    
  }
  
  
  
  # plot GC
  cat(" plot GC;")
  window.starts.gc = genomic.bins.starts(start = 1, end = chromosomes_lengths[j], bin.size = bin.size.plots.for.GC)
  gc.values = calculate.GC.in.windows.2(windows.starts = window.starts.gc, sequence = assembly.seq[[which(names(assembly.seq) == chromosomes[j])]], bin.size = bin.size.plots.for.GC)
  lines(window.starts.gc, gc.values, type = "l", xlim = c(1, chromosomes_lengths[j]), #col = "#FFA500",
        xaxt = "n", ylab = "", yaxt = "n", xlab = "", lwd = 1, ylim = c(0,100))
  
  cat(" plot repeats family;")
  # plot Family repeats
  if(nrow(classes) > 0)
  {
    for(k in 1 : nrow(classes))
    {
      repeats.family = repeats.chromosome[repeats.chromosome$new_class == classes$class[k],]
      if(nrow(repeats.family) == 0) next
      repeat.class.frequencies = calculate.repeats.percentage.in.windows(windows.starts = window.starts, 
                                                                   repeat.starts = repeats.family$start, 
                                                                   repeat.lengths = repeats.family$width, 
                                                                   sequence.length = chromosomes_lengths[j])
      repeat.class.frequencies[repeat.class.frequencies == 0] = NA
      repeat.class.frequencies[repeat.class.frequencies > 100] = 100
      lines(window.starts, repeat.class.frequencies, pch = 16, type = "o", col = safe_colorblind_palette[k])
    }
  }
  
  
  #add bottom axis
  axis(1, at = seq(1, chromosomes_lengths[j], by = plot.x.ticks), labels = F, col.ticks = "black", lwd = 1)
  
  #add sequence names
  text(1, bin.size.plots, labels = chromosomes[j])
  # ====================================================
  
  ### plot B =========================================== 
  cat(" plot B;")
  sequence_arrays = arrays[arrays$seqID == chromosomes[j], ]
  sequence_edta = edta_cls[edta_cls$V1 == chromosomes[j], ]
  sequence_repeats = repeats[repeats$seqID == chromosomes[j], ]
  sequence_genes <- genes[genes$V1 == chromosomes[j],]
  # start with EDTA
  
  window.starts.edta = genomic.bins.starts(start = 1, end = chromosomes_lengths[j], bin.size = bin.size.plots.for.EDTA)
  
  cat(" edta freq;")
  edta.frequencies = calculate.repeats.percentage.in.windows(windows.starts = window.starts.edta, 
                                                               repeat.starts = sequence_edta$V4, 
                                                               repeat.lengths = sequence_edta$width, 
                                                               sequence.length = chromosomes_lengths[j])
  edta.frequencies[edta.frequencies == 0] = NA
  edta.frequencies[edta.frequencies > 100] = 100
  plot(window.starts.edta, edta.frequencies, type = "h", frame.plot = F, 
       xlim = c(1, chromosomes_lengths[j]), col = "#CCCCCC", 
       ylab = "", xlab = "", xaxt = "n", yaxt = "n", ylim = c(0,100))
  
  if(add_centromeric_files_to_plots) {
    # add centromere background colour
    starts <- centromere_data$pericentromere_start[j]
    ends <- centromere_data$pericentromere_end[j]
    for(as in seq_along(starts)) {
      rect(starts[as], 0, ends[as], 100, col = "#00000030", border = NULL)
    }
  }
    
    
  
  for(k in length(edta_classes) : 1) {
    class_edta <- sequence_edta[sequence_edta$V3 %in% edta_classes[[k]], ]
    if(nrow(class_edta) == 0) next
    edta.class.frequencies = calculate.repeats.percentage.in.windows(windows.starts = window.starts.edta, 
                                                                 repeat.starts = class_edta$V4, 
                                                                 repeat.lengths = class_edta$width, 
                                                                 sequence.length = chromosomes_lengths[j])
    edta.class.frequencies[edta.class.frequencies == 0] = NA
    edta.class.frequencies[edta.class.frequencies > 100] = 100
    lines(window.starts.edta, edta.class.frequencies, pch = 16, type = "o", col = edta_classes_colours[k])
    
  }
  
  
  
  
  # gene and TE peaks/valleys
  
  gr2 <- with(sequence_edta, GRanges(chromosomes[j], IRanges(V4, V5)))
  
  ### Calculate TE and Gene landscape =====================================================
  print("calculate TE repeat and gene landscapes")

  # Make a TE plus repeats all bp position map
  TE_coordinates <- list()
  for(edta_id in seq_len(nrow(sequence_edta))) {
    TE_coordinates <- append(TE_coordinates, list(sequence_edta$V4[edta_id] : sequence_edta$V5[edta_id]))
  }
  
  for(repeat_id in seq_len(nrow(sequence_repeats))) {
    TE_coordinates <- append(TE_coordinates, list(sequence_repeats$start[repeat_id] : sequence_repeats$end[repeat_id]))
  }
  
  
  TE_coordinates <- unlist(TE_coordinates)
  TE_coordinates <- c(0, TE_coordinates, chromosomes_lengths[j])
  length(TE_coordinates)
  length(unique(TE_coordinates))
  
  # make a bp histogram to plot
  hist_EDTA <- hist(TE_coordinates, breaks = seq(min(TE_coordinates), max(TE_coordinates), length.out = 50), plot = FALSE)
  counts <- c(hist_EDTA$counts[1], hist_EDTA$counts[1], hist_EDTA$counts, hist_EDTA$counts[length(hist_EDTA$counts)], hist_EDTA$counts[length(hist_EDTA$counts)])
  ma_values_edta <- ma(counts)[3 : (length(counts) - 2)]
  edta_peak <- hist_EDTA$mids[which.max(ma_values_edta)]
  
  par(new=TRUE)
  plot(y = ma_values_edta, x = hist_EDTA$mids, type = "b", lwd = 4, ylim = c(0, max(ma_values_edta)),
       main = "", col = "#0066aa", yaxt = "n", ylab = "")
  axis(side = 2, col = "#0066aa")
  
  # Now calculations for GENES
  
  gene_coordinates <- list()
  for(edta_id in seq_len(nrow(sequence_genes))) {
    gene_coordinates <- append(gene_coordinates, list(sequence_genes$V4[edta_id] : sequence_genes$V5[edta_id]))
  }
  gene_coordinates <- unlist(gene_coordinates)
  gene_coordinates <- c(0, gene_coordinates, chromosomes_lengths[j])
  length(gene_coordinates)
  length(unique(gene_coordinates))
  
  hist_gene <- hist(gene_coordinates, breaks = seq(min(gene_coordinates), max(gene_coordinates), length.out = 50), plot = FALSE)
  counts <- c(hist_gene$counts[1], hist_gene$counts[1], hist_gene$counts, hist_gene$counts[length(hist_gene$counts)], hist_gene$counts[length(hist_gene$counts)])
  ma_values_genes <- ma(counts)[3 : (length(counts) - 2)]
  gene_valley <- hist_gene$mids[which.min(ma_values_genes)]
  par(new=TRUE)
  plot(y = ma_values_genes, x = hist_gene$mids, type = "b", lwd = 4, ylim = c(0, max(ma_values_genes)),
       ylab = "", yaxt = "n", main = "", col = "#00bb33")
  axis(side = 4, col = "#00bb33")
  
  # ====================================================
  
}
# Add dotplot of consensus sequences of all families and reverse complementary too

plot.consensus.divisions = nchar(classes$consensus[1])
if(nrow(classes) > 1) for(j in 2 : nrow(classes)) 
  plot.consensus.divisions = c(plot.consensus.divisions, 
                               (plot.consensus.divisions[length(plot.consensus.divisions)] + nchar(classes$consensus[j])))
midline = plot.consensus.divisions[length(plot.consensus.divisions)]

if(nrow(classes) > 1) for(j in nrow(classes) : 1) 
  plot.consensus.divisions = c(plot.consensus.divisions, 
                               (plot.consensus.divisions[length(plot.consensus.divisions)] + nchar(classes$consensus[j])))

sequence = strsplit(paste(classes$consensus[1:nrow(classes)], collapse = ""), split = "")[[1]]
sequence = c(sequence, strsplit(revCompString(paste(sequence, collapse = "")), split = "")[[1]])

dotPlot(sequence, sequence,
        wsize = 4, wstep = 1, nmatch = 4, col = c("white", "black"),
        xlab = "n", ylab = "n", cex = 10, yaxt="n", xaxt="n")

abline(v = c(1, plot.consensus.divisions), col = safe_colorblind_palette[1 : length(plot.consensus.divisions)], lwd = 12)
abline(h = c(1, plot.consensus.divisions), col = safe_colorblind_palette[1 : length(plot.consensus.divisions)], lwd = 12)

abline(v = plot.consensus.divisions[order(plot.consensus.divisions, decreasing = T)][1:nrow(classes)], 
       col = safe_colorblind_palette[1 : nrow(classes)], lwd = 12)
abline(h = plot.consensus.divisions[order(plot.consensus.divisions, decreasing = T)][1:nrow(classes)], 
       col = safe_colorblind_palette[1 : nrow(classes)], lwd = 12)
abline(h = midline, v = midline, col = "black", lwd = 24)

create_table(classes[,c(1,2,4,7)], c("Class", "Repeats no", "Median width", "Total width"), 
             colours = safe_colorblind_palette[1 : nrow(classes)], font_size = 4)

dev.off()

system2(command = "cp", args = c(paste0("./global_plot_plus_edta_", assembly_name, ".png"), 
                                 paste0("../../../plots_repeats_classes_gc_edta_gene_landscape/global_plot_plus_edta_", assembly_name, ".png")))

print(paste0("Assembly ", assembly_name, " done"))


