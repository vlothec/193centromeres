
.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))
suppressMessages(library(msa))
suppressMessages(library(seqinr))
suppressMessages(library(GenomicRanges))
suppressMessages(library(ggplot2))
library(ape)
library(nlme)
library(splines)


ma <- function(x, n = 5){filter(x, rep(1 / n, n), sides = 2)}


setwd("/home/pwlodzimierz/ToL/git_ToL")
source("./aux_fun.R")
setwd("/home/pwlodzimierz/ToL/upload_files/grand_tables")

grand_tables <- list.files(path = "./", pattern = "FULL_")

grand_table <- NULL
for(i in seq_along(grand_tables)) {
  grand_table <- rbind(grand_table, read.csv(file = grand_tables[i]))
}
# Read the tree
tree_full <- read.tree("../../IIkjj8ySVMo1lh6NE5YX8A_newick.txt")

grand_table$tree_name_species <- ""
for(i in 1 : nrow(grand_table)) {
  grand_table$tree_name_species[i] <- paste0(strsplit(grand_table$genus[i], split = "")[[1]][1], "._", grand_table$species[i])
}


for(i in seq_along(grand_tables)) {
  identifier <- strsplit(grand_tables[i], split = "_")[[1]][2]
  print(paste0(i, " ", sum(grepl(identifier, grand_table$fasta_name))))
}

grand_table$col <- "#00000099"
grand_table$col[grand_table$architecture == "Satellite"] <- "#ff00ff99"
grand_table$col[grand_table$architecture == "Transposon"] <- "#ffff0099"
grand_table$col[grand_table$architecture == "Holocentric no satellite"] <- "#00ffff99"
grand_table$col[grand_table$architecture == "Holocentric satellite"] <- "#00ffff99"


write.csv(grand_table, file = "/home/pwlodzimierz/ToL/upload_files/grand_tables/grandest_table.csv", row.names = FALSE)
# grand_table <- read.csv(file = "/home/pwlodzimierz/ToL/upload_files/grand_tables/grandest_table.csv")

for(i in 1 : nrow(grand_table)) {
  if(grepl("centromeric_repeat_", names(grand_table)[i])) {
    names(grand_table)[i] <- paste0("cen_rep_", strsplit(names(grand_table)[i], split = "centromeric_repeat_")[[1]][2], collapse = "")
  }
}


# > colnames(grand_table)
# [1] "fasta_name"
# [2] "group"
# [3] "genus"
# [4] "species"
# [5] "architecture"
# [6] "typic"
# [7] "centromeric_satellites.1_yes_0_no"
# [8] "holocentricity.1_yes_0_no"
# [9] "sequences.no"
# [10] "chromosomes.no"
# [11] "chromosomes_with_centromeric_repeats.no"
# [12] "size_genome.bp"
# [13] "size_chromosomes.bp"
# [14] "size_non_chromosomes.bp"
# [15] "chromosomes_mean_size.bp"
# [16] "chromosomes_mean_size_SD.bp"
# [17] "size_centromeres.bp"
# [18] "size_pericentromeres.bp"
# [19] "centromere_mean_size.bp"
# [20] "centromere_mean_size_SD.bp"
# [21] "GC_genome.perc"
# [22] "GC_chromosomes.perc"
# [23] "GC_non_chromosomes.perc"
# [24] "GC_centromeres.perc"
# [25] "GC_pericentromeres.perc"
# [26] "GC_arms.perc"
# [27] "GC_exons_total.perc"
# [28] "GC_TE_total.perc"
# [29] "GC_repeats_total.perc"
# [30] "GC_repeats_centromeric.perc"
# [31] "GC_repeats_noncentromeric.perc"
# [32] "genes_chr.no"
# [33] "genes_chr.bp"
# [34] "genes_chr.perc"
# [35] "exons_chr.no"
# [36] "exons_chr.bp"
# [37] "exons_chr.perc"
# [38] "exons_centromeres.no"
# [39] "exons_centromeres.bp"
# [40] "exons_centromeres.perc"
# [41] "exons_pericentromeres.no"
# [42] "exons_pericentromeres.bp"
# [43] "exons_pericentromeres.perc"
# [44] "exons_arms.no"
# [45] "exons_arms.bp"
# [46] "exons_arms.perc"
# [47] "TE_total.no"
# [48] "TE_total.bp"
# [49] "TE_total.perc"
# [50] "TE_C1_LTR.no"
# [51] "TE_C1_LTR.bp"
# [52] "TE_C1_LTR.perc"
# [53] "TE_C1_nonLTR.no"
# [54] "TE_C1_nonLTR.bp"
# [55] "TE_C1_nonLTR.perc"
# [56] "TE_C2_TIR.no"
# [57] "TE_C2_TIR.bp"
# [58] "TE_C2_TIR.perc"
# [59] "TE_C2_nonTIR.no"
# [60] "TE_C2_nonTIR.bp"
# [61] "TE_C2_nonTIR.perc"
# [62] "TE_LINE.no"
# [63] "TE_LINE.bp"
# [64] "TE_LINE.perc"
# [65] "TE_HELITRON.no"
# [66] "TE_HELITRON.bp"
# [67] "TE_HELITRON.perc"
# [68] "TE_total_centromere.perc"
# [69] "TE_total_pericentromere.perc"
# [70] "TE_total_arms.perc"
# [71] "TE_COPIA_centromere.perc"
# [72] "TE_COPIA_pericentromere.perc"
# [73] "TE_COPIA_arms.perc"
# [74] "TE_GYPSY_centromere.perc"
# [75] "TE_GYPSY_pericentromere.perc"
# [76] "TE_GYPSY_arms.perc"
# [77] "TE_HELITRON_centromere.perc"
# [78] "TE_HELITRON_pericentromere.perc"
# [79] "TE_HELITRON_arms.perc"
# [80] "TE_C1_LTR_centromere.perc"
# [81] "TE_C1_LTR_pericentromere.perc"
# [82] "TE_C1_LTR_arms.perc"
# [83] "TE_C1_nonLTR_centromere.perc"
# [84] "TE_C1_nonLTR_pericentromere.perc"
# [85] "TE_C1_nonLTR_arms.perc"
# [86] "TE_C2_TIR_centromere.perc"
# [87] "TE_C2_TIR_pericentromere.perc"
# [88] "TE_C2_TIR_arms.perc"
# [89] "TE_C2_nonTIR_centromere.perc"
# [90] "TE_C2_nonTIR_pericentromere.perc"
# [91] "TE_C2_nonTIR_arms.perc"
# [92] "repeats_total.no"
# [93] "repeats_total.bp"
# [94] "repeats_total.perc"
# [95] "centromeric_families.no"
# [96] "repeats_centromeric.no"
# [97] "repeats_centromeric.bp"
# [98] "repeats_centromeric.perc"
# [99] "repeats_noncentromeric.no"
# [100] "repeats_noncentromeric.bp"
# [101] "repeats_noncentromeric.perc"
# [102] "repeats_centromeric_on_chromosomes.bp"
# [103] "repeats_centromeric_not_on_chromosomes.bp"
# [104] "repeats_centromeric_not_on_chromosomes_vs_all_centromeric.perc"
# [105] "centromeric_repeats_fraction_on_dominant_strand.perc"
# [106] "centromeric_repeats_in_centromere.perc"
# [107] "cen_rep_1_name"
# [108] "cen_rep_1_TRASH_identifier"
# [109] "cen_rep_1_consensus_unambiguous"
# [110] "cen_rep_1_consensus_ambiguous"
# [111] "cen_rep_1_mean_size.bp"
# [112] "cen_rep_1_mean_size_SD.bp"
# [113] "cen_rep_1_mean_size_SD.perc"
# [114] "cen_rep_1_mean_pairwise_similarity.perc"
# [115] "cen_rep_1_mean_pairwise_similarity_SD.perc"
# [116] "cen_rep_1_similarity_between_chromosomes_mean_pairwise.perc"
# [117] "cen_rep_1_similarity_between_chromosomes_mean_pairwise_SD.perc_point"
# [118] "cen_rep_1_similarity_within_chromosomes_mean_pairwise.perc"
# [119] "cen_rep_1_similarity_within_chromosomes_mean_pairwise_SD.perc_point"
# [120] "cen_rep_1_GC.perc"
# [121] "cen_rep_1_mean_HOR_score"
# [122] "cen_rep_1_hors_no"
# [123] "cen_rep_1_hors_mean_block_size"
# [124] "cen_rep_1_hors_mean_block_distance"
# [125] "cen_rep_2_name"
# [126] "cen_rep_2_TRASH_identifier"
# [127] "cen_rep_2_consensus_unambiguous"
# [128] "cen_rep_2_consensus_ambiguous"
# [129] "cen_rep_2_mean_size.bp"
# [130] "cen_rep_2_mean_size_SD.bp"
# [131] "cen_rep_2_mean_size_SD.perc"
# [132] "cen_rep_2_mean_pairwise_similarity.perc"
# [133] "cen_rep_2_mean_pairwise_similarity_SD.perc_point"
# [134] "cen_rep_2_similarity_between_chromosomes_mean_pairwise.perc"
# [135] "cen_rep_2_similarity_between_chromosomes_mean_pairwise_SD.perc_point"
# [136] "cen_rep_2_similarity_within_chromosomes_mean_pairwise.perc"
# [137] "cen_rep_2_similarity_within_chromosomes_mean_pairwise_SD.perc_point"
# [138] "cen_rep_2_GC.perc"
# [139] "cen_rep_2_mean_HOR_score"
# [140] "cen_rep_2_hors_no"
# [141] "cen_rep_2_hors_mean_block_size"
# [142] "cen_rep_2_hors_mean_block_distance"
# [143] "gaps_in_centromeric_arrays.no"
# [144] "gaps_in_centromeric_arrays.total_bp"
# [145] "gaps_in_centromeric_arrays_size_against_arrays.perc"
# [146] "gaps_in_centromeric_arrays_size_against_arrays_no_edge_arrays.perc"
# [147] "gaps_in_centromeric_arrays.mean_bp"
# [148] "gaps_in_centromeric_arrays.median_bp"
# [149] "gaps_in_centromeric_arrays_TE_coverage.perc"
# [150] "gaps_in_centromeric_arrays_LTR_coverage.perc"
# [151] "gaps_in_centromeric_arrays_no_edge_arrays.no"
# [152] "gaps_in_centromeric_arrays_no_edge_arrays.total_bp"
# [153] "gaps_in_centromeric_arrays_no_edge_arrays.mean_bp"
# [154] "gaps_in_centromeric_arrays_no_edge_arrays.median_bp"
# [155] "gaps_in_centromeric_arrays_no_edge_arrays_TE_coverage.perc"
# [156] "gaps_in_centromeric_arrays_no_edge_arrays_LTR_coverage.perc"
# [157] "CENPA_presence.1_yes_0_no"


# clean_data <- clean_data[clean_data$genus != "Trichoderma",]
# clean_data <- clean_data[clean_data$genus != "Bibio",]


# model <- gls(response_var = "TE_total.perc" (dependent variable, Y-axis) ~ predictor_var = "size_genome.bp" (independent variable, X-axis), ...)

# 
# fitted(model) returns values on the scale of TE_total.perc — the predicted % TE content.

model_and_plot = function(response_var, predictor_var, data_frame, tree_full, remove_outliers = TRUE, plot_model = FALSE, significance = 0.05,
                          transform = "non", add_xisy_line = FALSE, plot_dir = "/home/pwlodzimierz/ToL/upload_files/plots_from_grand_table") {
  
  which.response_var <- which(response_var == names(data_frame))
  which.predictor_var <- which(predictor_var == names(data_frame))
  
  clean_data <- data_frame[!is.na(data_frame[, which.response_var]) & 
                             !is.na(data_frame[, which.predictor_var]), ]
  
  common_species <- intersect(tree_full$tip.label, clean_data$tree_name_species)
  tree <- suppressMessages(drop.tip(tree_full, setdiff(tree_full$tip.label, common_species)))
  if(is.null(tree)) {
    return(list(NA,NA))
  }
  data_sub <- clean_data[clean_data$tree_name_species %in% common_species, ]
  
  data_sub <- data_sub[match(tree$tip.label, data_sub$tree_name_species), ]
  
  
  if(transform != "non") {
    if(transform == "Quadratic") {
      data_sub[[paste0(predictor_var, "_test")]] <- data_sub[[predictor_var]]^2
      clean_data[[paste0(predictor_var, "_test")]] <- clean_data[[predictor_var]]^2
    }
    if(transform == "Cubic") {
      data_sub[[paste0(predictor_var, "_test")]] <- data_sub[[predictor_var]]^3
      clean_data[[paste0(predictor_var, "_test")]] <- clean_data[[predictor_var]]^3
    }
    if(transform == "Logarithmic") {
      data_sub[[paste0(predictor_var, "_test")]] <- log10(data_sub[[predictor_var]])
      clean_data[[paste0(predictor_var, "_test")]] <- log10(clean_data[[predictor_var]])
    }
    if(transform == "Square root") {
      data_sub[[paste0(predictor_var, "_test")]] <- sqrt(data_sub[[predictor_var]])
      clean_data[[paste0(predictor_var, "_test")]] <- sqrt(clean_data[[predictor_var]])
    }
    if(transform == "Inverse") {
      data_sub[[paste0(predictor_var, "_test")]] <- 1/data_sub[[predictor_var]]
      clean_data[[paste0(predictor_var, "_test")]] <- 1/clean_data[[predictor_var]]
    }
    if(transform == "Exponential") {
      data_sub[[paste0(predictor_var, "_test")]] <- exp(data_sub[[predictor_var]])
      clean_data[[paste0(predictor_var, "_test")]] <- exp(clean_data[[predictor_var]])
      if(sum(!is.infinite(data_sub[[paste0(predictor_var, "_test")]])) == 0) return("Do not use exponential with such high values")
    }
    if(transform == "Polynomial") {
      data_sub[[paste0(predictor_var, "_test")]] <- data_sub[[predictor_var]] + data_sub[[predictor_var]]^2
      clean_data[[paste0(predictor_var, "_test")]] <- clean_data[[predictor_var]] + clean_data[[predictor_var]]^2
    }
  } else {
    data_sub[[paste0(predictor_var, "_test")]] <- data_sub[[predictor_var]]
    clean_data[[paste0(predictor_var, "_test")]] <- clean_data[[predictor_var]]
  }
  
  cor_struct <- corBrownian(phy = tree, form = ~tree_name_species)
  
  if(sum(clean_data[,which.response_var] != clean_data[,which.predictor_var]) == 0) {
    return(list(1,0))
  }
  
  
  model <- nlme::gls(as.formula(paste(response_var, "~", paste0(predictor_var, "_test"))),
                     data = data_sub,
                     correlation = cor_struct,
                     method = "ML")
  
  if(remove_outliers) {
    std_resid <- resid(model, type = "pearson")
    outlier_threshold <- 2  # or 3 for stricter
    outliers <- which(abs(std_resid) > outlier_threshold)
    outlier_info <- data_sub[outliers, ]
    outlier_info$std_resid <- std_resid[outliers]
    cat("Outliers: ", data_sub$tree_name_species[outliers], "\n")
    
    if(length(outliers) != 0) {
      clean_data$col[clean_data$tree_name_species %in% data_sub$tree_name_species[outliers]] <- "#eeeeee99"
      
      data_sub <- data_sub[-outliers, ]
      
      model <- nlme::gls(as.formula(paste(response_var, "~", paste0(predictor_var, "_test"))),
                         data = data_sub,
                         correlation = cor_struct,
                         method = "ML")
      
    }
    
  }
  
  significant <- FALSE
  if(coefficients(summary(model))[8] < significance) significant <- TRUE
  
  pdf(file = paste0(plot_dir, "/", response_var, "_vs_", predictor_var, 
                    "_outl:", !remove_outliers, "_sig:", significant, "_trans:", transform, ".pdf"))
  plot(x = clean_data[[paste0(predictor_var, "_test")]], 
       y = clean_data[, which.response_var], 
       pch = 16,
       col = clean_data$col, 
       main = paste0(transform, " transformed ", response_var, " vs ", predictor_var), cex.main = 0.7,
       xlab = paste(transform, " transformed ", predictor_var),
       ylab = response_var)
  points(x = clean_data[[paste0(predictor_var, "_test")]], 
         y = clean_data[, which.response_var], pch = 1)
  if(remove_outliers) {
    legend("bottomright", legend = c("Satellite", "Transposon", "Holocentric", "unknown", "outlier"), 
           fill = c("#ff00ff", "#ffff00", "#00ffff", "#00000099", "#eeeeee99"),
           border = c("black", "black", "black", "black", "black"))
  } else {
    legend("bottomright", legend = c("Satellite", "Transposon", "Holocentric", "unknown"), 
           fill = c("#ff00ff", "#ffff00", "#00ffff", "#00000099"),
           border = c("black", "black", "black", "black"))
  }
  if(add_xisy_line) {
    lines(x = c(0,100), y = c(0,100), lty = 2)
  }
  
  mtext(paste0("Phylogeny corrected coefficient: ", format(coefficients(summary(model))[2], scientific = TRUE, digits = 3), " p-value: ", format(coefficients(summary(model))[8], scientific = TRUE, digits = 3)), 
        side = 3, line = 0.7, cex = 0.8)
  mtext(paste0("Non corrected Pearsons test p-value: ", format(cor.test(clean_data[[paste0(predictor_var, "_test")]], 
                                                                        clean_data[, which.response_var], method = 'pearson')$p.value, scientific = TRUE, digits = 3)), 
        side = 3, line = 0, cex = 0.8)
  dev.off()
  
  if(plot_model) {
    pdf(paste0(plot_dir, "/", response_var, "_vs_", predictor_var, 
               "_model_outl:", !remove_outliers, "_sig:", significant, "_trans:", transform, ".pdf"), onefile = T)
    par(mfrow = c(2, 2))
    plot(fitted(model), resid(model, type = "pearson"),
         xlab = "Fitted values",
         ylab = "Standardized residuals",
         main = paste0(transform, " transformed Residuals vs Fitted"))
    abline(h = 0, lty = 2, col = "gray")
    qqnorm(resid(model, type = "pearson"), main = "Normal Q-Q")
    qqline(resid(model, type = "pearson"), col = "gray")
    plot(fitted(model), abs(resid(model, type = "pearson")),
         xlab = "Fitted values",
         ylab = "Absolute residuals",
         main = "Scale-Location")
    abline(h = 0, lty = 2, col = "gray")
    hist(resid(model, type = "pearson"),
         breaks = 20,
         main = "Histogram of Residuals",
         xlab = "Standardized residuals",
         col = "lightgray")
    dev.off()
  }
  
  return(list(coefficients(summary(model))[2], coefficients(summary(model))[8]))
  
}

## available transforms:
# Quadratic
# Cubic
# Logarithmic
# Square root
# Exponential
# Polynomial
# Inverse



################################################################################
### gene, transposon and repeat content vs genome size
################################################################################
# model_and_plot(response_var = "exons_chr.perc", 
#                predictor_var = "size_genome.bp", 
#                data_frame = grand_table, 
#                tree_full = tree_full, 
#                remove_outliers = T,
#                plot_model = T,
#                transform = "non")

grand_table$log10_exons_chr.perc <- log10(grand_table$exons_chr.perc) 

model_and_plot(response_var = "log10_exons_chr.perc", 
               predictor_var = "size_genome.bp", 
               data_frame = grand_table, 
               tree_full = tree_full, 
               remove_outliers = T,
               plot_model = T,
               transform = "Logarithmic")

# model_and_plot(response_var = "repeats_total.perc", 
#                predictor_var = "size_genome.bp", 
#                data_frame = grand_table, 
#                tree_full = tree_full, 
#                remove_outliers = T,
#                plot_model = T,
#                transform = "non")

grand_table$log10_repeats_total.perc <- log10(grand_table$repeats_total.perc) 

model_and_plot(response_var = "log10_repeats_total.perc", 
               predictor_var = "size_genome.bp", 
               data_frame = grand_table, 
               tree_full = tree_full, 
               remove_outliers = T,
               plot_model = T,
               transform = "Logarithmic")

# model_and_plot(response_var = "TE_total.perc", 
#                predictor_var = "size_genome.bp", 
#                data_frame = grand_table, 
#                tree_full = tree_full, 
#                remove_outliers = T,
#                plot_model = T,
#                transform = "non")

grand_table$log10_TE_total.perc <- log10(grand_table$TE_total.perc) 

model_and_plot(response_var = "log10_TE_total.perc", 
               predictor_var = "size_genome.bp", 
               data_frame = grand_table, 
               tree_full = tree_full, 
               remove_outliers = T,
               plot_model = T,
               transform = "Logarithmic")

################################################################################



################################################################################
### GC of the genome vs of repeats and centromeric repeats
################################################################################
model_and_plot(response_var = "GC_repeats_total.perc", 
               predictor_var = "GC_genome.perc", 
               data_frame = grand_table, 
               tree_full = tree_full, 
               remove_outliers = TRUE,
               plot_model = T,
               transform = "non",
               add_xisy_line = TRUE)

model_and_plot(response_var = "cen_rep_1_GC.perc", 
               predictor_var = "GC_genome.perc", 
               data_frame = grand_table, 
               tree_full = tree_full, 
               remove_outliers = TRUE,
               plot_model = T,
               transform = "non",
               add_xisy_line = TRUE)


################################################################################



################################################################################
### cen_rep_1_mean_HOR_score and various parameters
################################################################################

model_and_plot(response_var = "cen_rep_1_mean_HOR_score", 
               predictor_var = "gaps_in_centromeric_arrays_size_against_arrays.perc", 
               data_frame = grand_table, 
               tree_full = tree_full, 
               remove_outliers = T,
               plot_model = T,
               transform = "non")

model_and_plot(response_var = "cen_rep_1_mean_HOR_score", 
               predictor_var = "gaps_in_centromeric_arrays_size_against_arrays.perc", 
               data_frame = grand_table, 
               tree_full = tree_full, 
               remove_outliers = T,
               plot_model = T,
               transform = "Logarithmic")

model_and_plot(response_var = "cen_rep_1_mean_HOR_score", 
               predictor_var = "gaps_in_centromeric_arrays_size_against_arrays.perc", 
               data_frame = grand_table, 
               tree_full = tree_full, 
               remove_outliers = T,
               plot_model = T,
               transform = "Square root")

model_and_plot(response_var = "cen_rep_1_mean_HOR_score", 
               predictor_var = "gaps_in_centromeric_arrays_TE_coverage.perc", 
               data_frame = grand_table, 
               tree_full = tree_full, 
               remove_outliers = T,
               plot_model = T,
               transform = "non")
################################################################################



################################################################################
### similarity within chromosome and between
################################################################################


model_and_plot(response_var = "cen_rep_1_similarity_within_chromosomes_mean_pairwise.perc", 
               predictor_var = "cen_rep_1_similarity_between_chromosomes_mean_pairwise.perc", 
               data_frame = grand_table, 
               tree_full = tree_full, 
               remove_outliers = T,
               plot_model = T,
               transform = "non",
               add_xisy_line = T)

################################################################################





################################################################################
### EVERYTHING
################################################################################


library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(reshape2)


numeric_columns <- c(9:106,111:124, 143:156)

matrix_coef <- matrix(data = NA, nrow = length(numeric_columns), ncol = length(numeric_columns))
matrix_pval <- matrix(data = NA, nrow = length(numeric_columns), ncol = length(numeric_columns))

for(i in 1 : length(numeric_columns[-length(numeric_columns)]))
{
  print(paste0(i, " / ", length(numeric_columns)))
  for(j in (i + 1) : length(numeric_columns)) {
    results <- model_and_plot(response_var = colnames(grand_table)[numeric_columns[i]], 
                              predictor_var = colnames(grand_table)[numeric_columns[j]], 
                              data_frame = grand_table, 
                              tree_full = tree_full, 
                              remove_outliers = F,
                              plot_model = F,
                              transform = "non",
                              plot_dir = "/home/pwlodzimierz/ToL/upload_files/plots_from_grand_table/all_plots")
    matrix_coef[i,j] <- results[[1]][1]
    matrix_pval[i,j] <- results[[2]][1]
  }
}

rownames(matrix_coef) <- colnames(grand_table)[numeric_columns]
colnames(matrix_coef) <- colnames(grand_table)[numeric_columns]
rownames(matrix_pval) <- colnames(grand_table)[numeric_columns]
colnames(matrix_pval) <- colnames(grand_table)[numeric_columns]

write.csv(matrix_coef, file = "/home/pwlodzimierz/ToL/upload_files/plots_from_grand_table/matrix_coef.csv")
write.csv(matrix_pval, file = "/home/pwlodzimierz/ToL/upload_files/plots_from_grand_table/matrix_pval.csv")




matrix_coef <- as.matrix(read.csv(file = "/home/pwlodzimierz/ToL/upload_files/plots_from_grand_table/matrix_coef.csv", row.names = 1))
matrix_pval <- as.matrix(read.csv(file = "/home/pwlodzimierz/ToL/upload_files/plots_from_grand_table/matrix_pval.csv", row.names = 1))

for(i in 1 : nrow(matrix_coef)) {
  matrix_coef[i,i] = 1
  matrix_pval[i,i] = 0
  for(j in 1 : nrow(matrix_coef)) {
    matrix_coef[j,i] = matrix_coef[i,j]
    matrix_pval[j,i] = matrix_pval[i,j]
  }
}


# Define color gradient for significant p-values
sig_colors <- colorRampPalette(c("#CCFFCC", "#006600"))(100)
grey_color <- "grey"
neg_colors <- colorRampPalette(c("#FFCCCC", "#660000"))(100)

# Create color matrix based on p-values and coefficients
color_matrix <- matrix(grey_color, nrow=nrow(matrix_pval), ncol=ncol(matrix_pval))
for(i in 1 : nrow(matrix_pval)) {
  for(j in 1 : ncol(matrix_pval)) {
    if(matrix_pval[i,j] <= 0.05) {
      p_scaled <- min(max((0.05 - matrix_pval[i,j]) / (0.05 - 0.001), 0), 1)
      if(matrix_coef[i,j] > 0) {
        color_matrix[i,j] <- sig_colors[round(p_scaled * 99) + 1]
      } else {
        color_matrix[i,j] <- neg_colors[round(p_scaled * 99) + 1]
      }
    }
  }
}


# Hierarchical clustering based on transformed p-values
transformed_pval <- -log10(matrix_pval + 1e-10)
transformed_pval <- 10 - transformed_pval
dist_matrix <- as.dist(transformed_pval)
hclust_result <- hclust(dist_matrix, method="complete")
order <- hclust_result$order
rownames(matrix_coef)[order]

# Reorder matrices
matrix_coef <- matrix_coef[order, order]
color_matrix <- color_matrix[order, order]

# Prepare data for ggplot
data <- melt(matrix_coef, varnames=c("Row", "Col"), value.name="Coef")
data$Pval <- melt(matrix_pval)$value
data$Color <- melt(color_matrix)$value
data$Label <- formatC(data$Coef, digits=0, format="e")

# Plot heatmap
pdf(file = "/home/pwlodzimierz/ToL/upload_files/plots_from_grand_table/matrix_coef_coloured_by_pval.pdf", width = 40, height = 40)
ggplot(data, aes(x=Col, y=Row, fill=Color, label=Label)) +
  geom_tile(color="white") +
  geom_text(size=2.4, angle=45) +
  scale_fill_identity() +
  theme_minimal() +
  labs(title="Heatmap of Coefficients with P-value Coloring") +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_blank())
dev.off()

png(filename = "/home/pwlodzimierz/ToL/upload_files/plots_from_grand_table/matrix_coef_coloured_by_pval.png", width = 8000, height = 8000)
ggplot(data, aes(x=Col, y=Row, fill=Color, label=Label)) +
  geom_tile(color="white") +
  geom_text(size=5, angle=45) +
  scale_fill_identity() +
  theme_minimal() +
  labs(title="Heatmap of Coefficients with P-value Coloring") +
  theme(axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_blank())
dev.off()

################################################################################















################################################################################
### 
################################################################################

model_and_plot(response_var = "centromeric_repeats_fraction_on_dominant_strand.perc", 
               predictor_var = "gaps_in_centromeric_arrays_size_against_arrays.perc", 
               data_frame = grand_table, 
               tree_full = tree_full, 
               remove_outliers = T,
               plot_model = T,
               transform = "non",
               add_xisy_line = F)





################################################################################











###########################################################
pdf(file = "boxplot exons coverage in centromeres, pericentromeres and arms.pdf")
boxplot(grand_table$exons_centromeres.perc, grand_table$exons_pericentromeres.perc, grand_table$exons_arms.perc, 
        names = c("exons cen %", "exons pericen %", "exons arms %"))
dev.off()


###########################################################
pdf(file = "boxplot TE coverage in centromeric gaps, centromeres, pericentromeres and arms.pdf")
boxplot(grand_table$gaps_in_centromeric_arrays_TE_coverage.perc[grand_table$gaps_in_centromeric_arrays_TE_coverage.perc > 0], grand_table$TE_total_centromere.perc, grand_table$TE_total_pericentromere.perc, grand_table$TE_total_arms.perc, 
        names = c("TE gaps %", "TE cen %", "TE pericen %", "TE arms %"))
dev.off()



###########################################################
pdf(file = "boxplots centromere pericentromere arms TE families.pdf", height = 12)
par(mfrow = c(3,1))
par(mar = c(8,4,4,2) + 0.1)
boxplot(grand_table$TE_C1_LTR_centromere.perc[is.finite(grand_table$TE_C1_LTR_centromere.perc)], 
        grand_table$TE_COPIA_centromere.perc[is.finite(grand_table$TE_COPIA_centromere.perc)], 
        grand_table$TE_GYPSY_centromere.perc[is.finite(grand_table$TE_GYPSY_centromere.perc)], 
        grand_table$TE_C1_nonLTR_centromere.perc[is.finite(grand_table$TE_C1_nonLTR_centromere.perc)], 
        grand_table$TE_C2_TIR_centromere.perc[is.finite(grand_table$TE_C2_TIR_centromere.perc)], 
        grand_table$TE_C2_TIR_centromere.perc[is.finite(grand_table$TE_C2_TIR_centromere.perc)],
        grand_table$TE_HELITRON_centromere.perc[is.finite(grand_table$TE_HELITRON_centromere.perc)], 
        names = c("Class 1 LTR", "COPIA", "GYPSY", "Class 1 non-LTR", "Class 2 TIR", "Class 2 non-TIR", "HELITRON"), 
        ylab = "% of centromeres covered by the TE family",
        main = "centromeres transposon % coverage",
        las = 2,
        ylim = c(0,75))

par(mar = c(8,4,4,2) + 0.1)
boxplot(grand_table$TE_C1_LTR_pericentromere.perc[is.finite(grand_table$TE_C1_LTR_pericentromere.perc)], 
        grand_table$TE_COPIA_pericentromere.perc[is.finite(grand_table$TE_COPIA_pericentromere.perc)], 
        grand_table$TE_GYPSY_pericentromere.perc[is.finite(grand_table$TE_GYPSY_pericentromere.perc)], 
        grand_table$TE_C1_nonLTR_pericentromere.perc[is.finite(grand_table$TE_C1_nonLTR_pericentromere.perc)], 
        grand_table$TE_C2_TIR_pericentromere.perc[is.finite(grand_table$TE_C2_TIR_pericentromere.perc)], 
        grand_table$TE_C2_TIR_pericentromere.perc[is.finite(grand_table$TE_C2_TIR_pericentromere.perc)],
        grand_table$TE_HELITRON_pericentromere.perc[is.finite(grand_table$TE_HELITRON_pericentromere.perc)], 
        names = c("Class 1 LTR", "COPIA", "GYPSY", "Class 1 non-LTR", "Class 2 TIR", "Class 2 non-TIR", "HELITRON"), 
        ylab = "% of pericentromeres covered by the TE family",
        main = "pericentromeres transposon % coverage",
        las = 2,
        ylim = c(0,75))

par(mar = c(8,4,4,2) + 0.1)
boxplot(grand_table$TE_C1_LTR_arms.perc[is.finite(grand_table$TE_C1_LTR_arms.perc)], 
        grand_table$TE_COPIA_arms.perc[is.finite(grand_table$TE_COPIA_arms.perc)], 
        grand_table$TE_GYPSY_arms.perc[is.finite(grand_table$TE_GYPSY_arms.perc)], 
        grand_table$TE_C1_nonLTR_arms.perc[is.finite(grand_table$TE_C1_nonLTR_arms.perc)], 
        grand_table$TE_C2_TIR_arms.perc[is.finite(grand_table$TE_C2_TIR_arms.perc)], 
        grand_table$TE_C2_nonTIR.perc[is.finite(grand_table$TE_C2_nonTIR.perc)],
        grand_table$TE_HELITRON_arms.perc[is.finite(grand_table$TE_HELITRON_arms.perc)], 
        names = c("Class 1 LTR", "COPIA", "GYPSY", "Class 1 non-LTR", "Class 2 TIR", "Class 2 non-TIR", "HELITRON"), 
        ylab = "% of arms covered by the TE family",
        main = "arms transposon % coverage",
        las = 2,
        ylim = c(0,75))
dev.off()




###########################################################
cor_val <- cor.test(grand_table$repeats_total.perc[grand_table$repeats_total.perc > 0], 
                    grand_table$size_chromosomes.bp[grand_table$repeats_total.perc > 0])

pdf(file = "genome size vs repeats fraction.pdf")
plot(x = grand_table$size_genome.bp[grand_table$repeats_total.perc > 0], 
     y = grand_table$repeats_total.perc[grand_table$repeats_total.perc > 0], 
     pch = 16, 
     col = grand_table$col[grand_table$repeats_total.perc > 0],
     main = "genome size vs repeats fraction",
     xlab = "genome size",
     ylab = "repeats fraction")
# text(2500000000, 20,paste0("p-value: ", round(cor_val$p.value, 2)))
# text(2500000000, 22,paste0("    cor: ", round(cor_val$estimate, 2)))
dev.off()




###########################################################
cor_val <- cor.test(grand_table$exons_chr.perc[grand_table$exons_chr.perc > 0], 
                    grand_table$size_chromosomes.bp[grand_table$exons_chr.perc > 0])

pdf(file = "genome size vs exons fraction.pdf")
plot(x = grand_table$size_chromosomes.bp[grand_table$exons_chr.perc > 0], 
     y = grand_table$exons_chr.perc[grand_table$exons_chr.perc > 0], 
     pch = 16, 
     col = grand_table$col[grand_table$exons_chr.perc > 0],
     main = "genome size vs exons fraction",
     xlab = "genome size",
     ylab = "exons fraction")
# text(2500000000, 20,paste0("p-value: ", formatC(cor_val$p.value, format = "e", digits = 2)))
# text(2500000000, 24,paste0("    cor: ", round(cor_val$estimate, 2)))
dev.off()




###########################################################
cor_val <- cor.test(grand_table$TE_total.perc[grand_table$TE_total.perc > 0], 
                    grand_table$size_chromosomes.bp[grand_table$TE_total.perc > 0])

pdf(file = "genome size vs transposon fraction.pdf")
plot(x = grand_table$size_chromosomes.bp[grand_table$TE_total.perc > 0], 
     y = grand_table$TE_total.perc[grand_table$TE_total.perc > 0], 
     pch = 16, 
     col = grand_table$col[grand_table$TE_total.perc > 0],
     main = "genome size vs transposon fraction",
     xlab = "genome size",
     ylab = "transposon fraction")
# text(2300000000, 45,paste0("p-value: ", formatC(cor_val$p.value, format = "e", digits = 2)))
# text(2300000000, 50,paste0("    cor: ", round(cor_val$estimate, 2)))
dev.off()





###########################################################
pdf(file = "repeat size SD vs repeat similarity.pdf")
plot(x = grand_table$cen_rep_1_mean_pairwise_similarity.perc, 
     y = grand_table$cen_rep_1_mean_size_SD.perc, 
     pch = 16, 
     col = grand_table$col,
     main = "repeat size SD vs repeat similarity",
     xlab = "mean_pairwise_similarity",
     ylab = "mean_size_SD")
dev.off()




###########################################################
pdf(file = "Genome sizes ordered.pdf")
plot(x = 1 : nrow(grand_table), 
     y = grand_table$size_genome.bp[order(grand_table$size_genome.bp, decreasing = F)], 
     pch = 16, 
     col = grand_table$col[order(grand_table$size_genome.bp, decreasing = F)],
     main = "Genome sizes vs chromosome number",
     xlab = "Genome",
     ylab = "Genome size bp")
legend(160,0,c("Satellite", "Transposon", "Holocentric"), c("#ff00ff", "#ffff00", "#00ffff"), xjust = 1, yjust = 0)
dev.off()




numeric_columns <- c(9:108,111:156)














