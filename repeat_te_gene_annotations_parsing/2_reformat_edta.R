#!/usr/bin/env Rscript
.libPaths(c(.libPaths(), "/home/pwlodzimierz/TRASH_dev/R_libs"))
suppressMessages(library(seqinr))
suppressMessages(library(msa))

setwd("/home/pwlodzimierz/ToL/git_ToL")
source("./aux_fun.R")

replace_existing_analysis = FALSE

data_directories <- list.dirs(path = "/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs", recursive = FALSE, full.names = TRUE)
data_directories <- data_directories[!grepl(pattern = "templated_", data_directories)]
data_directories <- data_directories[grepl(pattern = ".fa", data_directories)]
assembly_files <- list.files(path = "/home/pwlodzimierz/ToL/Assemblies/fastas_2021_Michael", recursive = FALSE, full.names = TRUE)
assembly_files <- assembly_files[!grepl(".fai", assembly_files)]
edta_files <- list.files(path = "/home/pwlodzimierz/ToL/EDTA_files", recursive = FALSE, full.names = TRUE)



taskid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
i = as.numeric(taskid)# 1 to 15
print(i)

# for(i in seq_along(data_directories)) 
{
  print(paste0("A: ", i, " / ", length(data_directories)))
  ### Load data ================================================================
  setwd(data_directories[i])
  
  assembly_name = strsplit(data_directories[i], split = "v2_out_for_HORs/")[[1]][2]
  
  if(assembly_name == "rosCan_S27_v1.fasta") {
    edta_og <- read.table(file = "./rosCan_S27_v1.fasta.F2B.noS1H2.mod.EDTA.TEanno.split.gff3.reassigned", header = FALSE, sep = "\t")
    edta_og <- edta_og[, c(1,4,5,8,6,7,2,3,8,9,10)]
    names(edta_og) = paste0("V", 1:11)
  } else {
    
    edta_file = grep(assembly_name, edta_files)
    
    if(!replace_existing_analysis) {
      if(file.exists(paste0(assembly_name, "_edta_modified.csv"))) quit(save = "no", status = 1)
    }
    if(length(edta_file) == 1) {
      edta_og = read.table(edta_files[edta_file], header = FALSE, sep = "\t")
    } else {
      print(" No edta found, cannot proceed")
      next
    }
    
    
  }
  
  edta_og_write <- edta_og[, c("V1", "V7", "V8", "V2", "V3", "V4", "V6", "V4", "V10", "V11")]
  write.table(x = edta_og_write, file = paste0(assembly_name, "_edta_original.gff"), 
              sep = "\t", row.names = F, col.names = F, quote = FALSE)
  
  ### Reformat EDTA gff table ==================================================
  v10 = lapply(edta_og$V10, function(X) strsplit(X, split = ";")[[1]])
  new_cols = lapply(v10, function(X) unlist(lapply(X, function(x) strsplit(x, split = "=")[[1]][1])))
  new_cols = unique(unlist(new_cols))
  edta_full = edta_og[, 1:9]
  for (j in seq_along(new_cols)) {
    cat(j, "/", length(new_cols), "\n")
    
    new_data = lapply(edta_og$V10, function(X) strsplit(X, split = paste0(new_cols[j], "=", collapse = ""))[[1]][2])
    new_data = unlist(lapply(new_data, function(X) strsplit(X, split = ";")[[1]][1]))
    
    edta_full = cbind(edta_full, new_data)
    names(edta_full)[ncol(edta_full)] = new_cols[j]
  }
  
  
  edta_full$oldV8 = edta_full$V8
  
  edta_repeat_region = edta_full[edta_full$V8 == "repeat_region", ]
  edta_non_rep_reg = edta_full[edta_full$V8 != "repeat_region", ]
  
  names = unique(edta_repeat_region$Name)
  names = sort(names)
  
  names = data.frame(names)
  names$short = unlist(lapply(names$names, function(X) {
    pos = grep("[^A-Za-z]", strsplit(X, split = "")[[1]])
    if(length(pos) == 0 ) return(X)
    return(substr(X, 1, min(pos) - 1))
    } ))
  
  names_short = unique(names$short)
  
  # for(j in seq_len(nrow(edta_repeat_region))) {
  #   print(j)
  #   edta_repeat_region$V8[j] = names$short[names$names == edta_repeat_region$Name[j]]
  # }
  edta_repeat_region$oldV8 = edta_repeat_region$V8
  edta_repeat_region$V8 = unlist(lapply(seq_len(nrow(edta_repeat_region)), function(X) names$short[names$names == edta_repeat_region$Name[X]]))
  
  edta_modified = rbind(edta_non_rep_reg, edta_repeat_region)
  edta_modified <- edta_modified[,c("V1", "V7", "V8", "V2", "V3", "V4", "V6", "V4", 
                                   "ID", "Name", "Classification", "Sequence_ontology", "Identity", "Method", "TSD", "TIR", "V4", "V4", "oldV8")]
  names(edta_modified) <- c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", 
                            "ID", "Name", "Classification", "Sequence_ontology", "Identity", "Method", "TSD", "TIR", "motif", "tsd", "oldV3")
  
  export_gff(annotations.data.frame = edta_modified, output = ".", 
             file.name = paste0(assembly_name, "_edta_modified"), seqid = 1, source = 2, type = 3, 
             start = 4, end = 5, strand = 7, score = ".")
  
  write.csv(edta_modified, paste0(assembly_name, "_edta_modified.csv"))
  
  
}












