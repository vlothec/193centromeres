

library(cowplot); library(ggnewscale); library(ggpubr); library(grid); library(gridExtra); library(RColorBrewer); library(tidyverse);

file.path <- "data/"

# data from cent_prop_combined.sh
list.files.v <- list.files(file.path, pattern = "_centromere_nt_prop_combined")

list.files.v <- list.files.v[!grepl("TRUE_NA",list.files.v)]
species.list <- list.files.v %>%
  {gsub("_.*","",.)} %>%
  {unique(.)}

cent.prop.raw <- data.frame()
for(species in 1:length(species.list)){
  species.i <- species.list[species]
  list.files.species <- list.files.v[grepl(species.i,list.files.v)]
  
  print(species.i)
  
  cent.prop.raw.i <- read.delim(paste0(file.path,"/",list.files.species), header=FALSE, comment.char = "")
  colnames(cent.prop.raw.i) <-c("nt.count","filter_edge","min.gap","side","feature","TE.cls","step")
  
  cent.prop.raw.i <- cent.prop.raw.i[which(cent.prop.raw.i$feature == "cent_nt"),c(2,1)] %>%
    distinct()
  cent.prop.raw.i$species <- species.i
  cent.prop.raw <- rbind(cent.prop.raw,cent.prop.raw.i)
}

species_order_list <- read.table("df_order.txt", quote="\"", comment.char="")
cent.prop.raw$species <- cent.prop.raw$species %>%
  {gsub("\\..*","",.)}

identical(sort(species_order_list$V1),sort(unique(cent.prop.raw$species)))
setdiff(sort(unique(cent.prop.raw$species)),sort(species_order_list$V1)) %>% length()


cent.prop.raw$species <- cent.prop.raw$species %>%
  {factor(.,levels = species_order_list$V1)}

cent.prop.raw %>% head
cent.prop.raw$x.start <- mapply(function(x) x/2,cent.prop.raw$nt.count) %>% {.*(-1)}
cent.prop.raw$x.end <- mapply(function(x) x/2,cent.prop.raw$nt.count)

# next line can be uncommented in case we want to sort yaxis by taxonomic order
{
  # upload phylo data
  data.attributes <- read.delim("dtol_plant_attributes.txt", header = T, sep = "\t")
  data.attributes <- data.attributes %>%
    mutate(species = gsub("\\..*","",fasta))
  
  data.attributes %>% str
  
  cent.prop.raw <- left_join(cent.prop.raw,data.attributes[,c(9,2,5)],by = "species")
  
  cent.prop.raw$group.0 <- mapply(function(x) case_when(x == "Dicots" ~ "plant",
                                                          x == "Monocots" ~ "plant",
                                                          x == "Aves" ~ "chordate",
                                                          x == "Mammalia" ~ "chordate",
                                                          x == "Actinopterygii" ~ "chordate",
                                                          x == "Hymenoptera" ~ "invertebrate",
                                                          x == "Annelida" ~ "invertebrate",
                                                          x == "Coleoptera" ~ "invertebrate",
                                                          x == "Diptera" ~ "invertebrate",
                                                          x == "Hemiptera" ~ "invertebrate",
                                                          x == "Cnidaria" ~ "invertebrate",
                                                          x == "Lepidoptera" ~ "invertebrate",
                                                          x == "Neuroptera" ~ "invertebrate",
                                                          x == "Tunicata" ~ "chordate",
                                                          x == "Reptilia" ~ "chordate",
                                                          x == "Bryozoa" ~ "invertebrate",
                                                          x == "Mollusca" ~ "invertebrate",
                                                          x == "Plecoptera" ~ "invertebrate",
                                                          x == "Echinodermata" ~ "invertebrate",
                                                          .default = NA_character_),cent.prop.raw$Group)
}


cent.prop.raw$group.0 <- cent.prop.raw$group.0 %>%
  {factor(.,levels = c("invertebrate","chordate","plant"))}

cent.prop.raw$col.tag <- "col"
cent.prop.raw$max.nt.count <- max(cent.prop.raw$nt.count)

cent.prop.raw %>% str

cent.prop.raw$species <- cent.prop.raw$species %>%
  {factor(.,levels = species_order_list$V1)}
cent.prop.raw <- cent.prop.raw[with(cent.prop.raw,order(species)),]

# text.size <- 8
text.size <- 6.5
bar.width <- 0.425

# no color code
r1.p2 <- ggplot(cent.prop.raw[which(cent.prop.raw$filter_edge == "TRUE_FALSE"),]) +
  geom_vline(xintercept = 150, linetype = "dotted", color = "transparent") +
  geom_col(aes(y = species, x = max.nt.count/1000000), fill = "grey95", width = bar.width, color = "transparent") +
  geom_col(aes(y = species, x = nt.count/1000000, fill = col.tag), width = bar.width, color = "transparent") +
  scale_fill_manual(values = c("col" = "black")) + 
  scale_y_discrete(expand = expansion(add = 0.85)) +
  scale_x_continuous(breaks = c(0,125,250), labels = c(0,125,250)) +
  facet_grid(group.0 ~ ., scales = "free_y", space = "free_y") +
  coord_cartesian(expand = T) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.key.size = unit(.35,"cm"),
        legend.spacing = unit(2,"cm"),
        legend.text = element_text(size = text.size),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "grey95"),
        plot.margin = unit(c(t = 5,r = 5,b = 5,l = 5),"mm"),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.title = element_blank(),
        axis.text.y = element_text(),
        strip.text.y = element_blank(),
        text = element_text(size = text.size))

pdf(paste0("cent_space.phylotaxa.pdf"),onefile = T, width = 1.75,height = 7.15)
print(r1.p2)
dev.off()

scf_name <- mapply(function (Genus,Species) paste0(gsub("(^[A-Za-z]{1}).*","\\1",Genus),". ", Species), data.attributes$Genus, data.attributes$Species)
data.attributes$scf_name <- scf_name

cent.prop.raw <- left_join(cent.prop.raw,data.attributes[,c(9,10)],by = "species")
cent.prop.raw$species <- cent.prop.raw$species %>%
  {factor(.,levels = species_order_list$V1)}
cent.prop.raw <- cent.prop.raw[with(cent.prop.raw,order(species)),]


# text.size <- 8
text.size <- 6.5
bar.width <- 0.425
# no color code
r1.p2 <- ggplot(cent.prop.raw[which(cent.prop.raw$filter_edge == "TRUE_FALSE"),]) +
  geom_vline(xintercept = 150, linetype = "dotted", color = "transparent") +
  geom_col(aes(y = fct_inorder(scf_name), x = nt.count/1000000, fill = col.tag), width = bar.width, color = "transparent") +
  scale_fill_manual(values = c("col" = "black")) + 
  scale_y_discrete(expand = expansion(add = 0.85)) +
  scale_x_continuous(breaks = c(0,125,250), labels = c(0,125,250)) +
  facet_grid(group.0 ~ ., scales = "free_y", space = "free_y") +
  coord_cartesian(expand = T) +
  theme_minimal() +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.key.size = unit(.35,"cm"),
        legend.spacing = unit(2,"cm"),
        legend.text = element_text(size = text.size),
        panel.background = element_rect(fill = "grey97", color = NA),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "grey95"),
        plot.margin = unit(c(t = 5,r = 5,b = 5,l = 5),"mm"),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        axis.title = element_blank(),
        axis.text.y = element_text(face = "italic"),
        strip.text.y = element_blank(),
        text = element_text(size = text.size))

pdf(paste0("cent_space.phylotaxa.scf_name.pdf"),onefile = T, width = 1.75,height = 7.15)
print(r1.p2)
dev.off()
