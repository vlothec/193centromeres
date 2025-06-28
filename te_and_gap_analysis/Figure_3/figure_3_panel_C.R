
library(cowplot); library(ggnewscale); library(ggpubr); library(grid); library(gridExtra); library(RColorBrewer); library(tidyverse);

file.path <- "data/"
# data from structural_homology_ratio.R
list.file.v <- list.files(path = file.path,
                          pattern = "_TRUE_FALSE_structural_homology_ratio$")

species.list <- gsub("_.*","",list.file.v)


df.complet <- data.frame()
for(species in 1:length(species.list)){
  species.i <- species.list[species]
  data.i <- read.delim(paste0(file.path,list.file.v[grepl(species.i,list.file.v)]))
  data.i$species <- species.i
  df.complet <- rbind(df.complet,data.i)
}

df.complet %>% colnames()

df.complet <- df.complet[with(df.complet,order(prop.nt.sum,decreasing = T)),]
df.complet$species <- df.complet$species %>%
  {gsub("\\..*","",.)}

species_order_list <- read.table("df_order.txt", quote="\"", comment.char="")

df.complet$species <- df.complet$species %>%
  {factor(.,levels = species_order_list$V1)}
df.complet <- df.complet[with(df.complet,order(species)),]


# next line can be uncommented in case we want to sort yaxis by taxonomic order
{
  # upload phylo data
  data.attributes <- read.delim("dtol_plant_attributes.txt", header = T, sep = "\t")
  data.attributes <- data.attributes %>%
    mutate(species = gsub("\\..*","",fasta))
  
  data.attributes %>% str
  
  df.complet <- left_join(df.complet,data.attributes[,c(9,2,5)],by = "species")
  
  df.complet$group.0 <- mapply(function(x) case_when(x == "Dicots" ~ "plant",
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
                                                     .default = NA_character_),df.complet$Group)
  
  df.complet$species <- df.complet$species %>%
    {factor(.,levels = species_order_list$V1)}
  df.complet$group.0 <- df.complet$group.0 %>%
    {factor(.,levels = c("invertebrate","chordate","plant"))}
  
}

taxa.mean <- df.complet %>%
  group_by(group.0,ft.method,TRUE_FALSE.red) %>%
  summarise(mean.cent = mean(prop.nt.sum)) %>%
  filter(ft.method == "structural")

all.mean <- df.complet %>%
  group_by(ft.method,TRUE_FALSE.red) %>%
  summarise(mean.cent = mean(prop.nt.sum)) %>%
  filter(ft.method == "structural")


# classic display
# text.size <- 8
text.size <- 6.5
bar.width <- 0.425

r1.p3 <- ggplot(df.complet) +
  geom_col(aes(x=prop.nt.sum, y=fct_inorder(species), fill = ft.method), color = "transparent", width = bar.width) +
  scale_fill_manual(values = c("structural" = "black", "homology" = "grey95")) +
  scale_x_continuous(breaks = c(0,0.5,1), labels = c(0,0.5,1)) +
  scale_y_discrete(expand = expansion(add = 0.85)) +
  facet_grid(group.0 ~ TRUE_FALSE.red, scales = "free_y", space = "free_y") +
  coord_cartesian(expand = T, xlim = c(0,1)) +
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
        axis.text.y = element_text(face = "plain"),
        strip.text = element_blank(),
        text = element_text(size = text.size))

pdf(paste0("structural_homology_ratio.phylotaxa.pdf"),onefile = T, width = 3,height = 7.15)
print(r1.p3)
dev.off()

scf_name <- mapply(function (Genus,Species) paste0(gsub("(^[A-Za-z]{1}).*","\\1",Genus),". ", Species), data.attributes$Genus, data.attributes$Species)
data.attributes$scf_name <- scf_name

df.complet <- left_join(df.complet,data.attributes[,c(9,10)],by = "species")
df.complet$species <- df.complet$species %>%
  {factor(.,levels = species_order_list$V1)}
df.complet <- df.complet[with(df.complet,order(species)),]


# text.size <- 8
text.size <- 6.5
bar.width <- 0.425

r1.p3 <- ggplot(df.complet) +
  geom_col(aes(x=prop.nt.sum, y=fct_inorder(scf_name), fill = ft.method), color = "transparent", width = bar.width, size = .25) +
  scale_fill_manual(values = c("structural" = "seagreen", "homology" = "#DDDDDD")) +
  scale_x_continuous(breaks = c(0,0.5,1), labels = c(0,0.5,1)) +
  scale_y_discrete(expand = expansion(add = 0.85)) +
  facet_grid(group.0 ~ TRUE_FALSE.red, scales = "free_y", space = "free_y") +
  coord_cartesian(expand = T, xlim = c(0,1)) +
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
        strip.text = element_blank(),
        text = element_text(size = text.size))

pdf(paste0("structural_homology_ratio.phylotaxa.te_color.scf_name.pdf"),onefile = T, width = 3,height = 7.15)
print(r1.p3)
dev.off()
