
library(seqinr)

chr_sizes <- read.csv(file = "/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full.Ian.csv")

species_to_analyse <- c("lpCarDepa1.1", "lpSchLacu1", "lpLuzSylv1.1", "iiLimLuna2.1", "iiLimMarm1.1", "iiLimRhom1.1")


repeats_to_analyse <- list("81_2", "183_4", c("124_1", "174_2"), "353_2", '166_3', "161_5")


min_distances <- c(70000, 60000, 50000, 111111, 70000, 70000, 25000)

chr_sizes <- read.csv(file = "/home/pwlodzimierz/ToL/Metadata/chr.no.and.sizes.full.Ian.csv")

# TODO: don't do scaffolds! For now, I removed them manually

#set up scripts

for(i in 1 : length(species_to_analyse)) {
  cat(i, "/", length(species_to_analyse), species_to_analyse[i],  "\n")
  
  
  for(irep in 1 : length(repeats_to_analyse[[i]])) {
    
    setwd(paste0("/home/pwlodzimierz/ToL/Repeats_HOR_TRASH/v2_out_for_HORs/", species_to_analyse[i], ".fa"))
    
    arrays <- read.csv(file = "./holocentric_arrays_data.csv")
    
    chr_sizes_genome <- chr_sizes[chr_sizes$assembly.name == paste0(species_to_analyse[i], ".fa"), ]
    chr_sizes_genome <- chr_sizes_genome[chr_sizes_genome$is.chr == 1, ]
    
    arrays <- arrays[arrays$chromosome %in% chr_sizes_genome$chromosome.name, ]
    
    arrays_chr1 <- arrays[arrays$chromosome == chr_sizes_genome$chromosome.name[1], ]
    
    pdf(file = paste0("/home/pwlodzimierz/ToL/upload_files/42_holocentric_plots_figure/chromosome_landscape_", species_to_analyse[i], ".pdf"), width = 12, height = 3)
    plot(NA,NA, xlim = c(0, chr_sizes_genome$size[1]), ylim = c(0,100), xlab = "Coordinates, Mbp", ylab = "", main = paste0(species_to_analyse[i], " ", repeats_to_analyse[[i]][irep]))
    
    for(j in 1 : nrow(arrays_chr1)) {
      rect(xleft = arrays_chr1$start[j], ybottom = 0, xright = arrays_chr1$end[j], ytop = 100, col = "black")
    }
    dev.off()
    
    
    
    
    
    
    
  }
}