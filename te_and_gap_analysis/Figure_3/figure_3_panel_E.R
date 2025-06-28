
library(cowplot); library(ggnewscale); library(ggpubr); library(grid); library(gridExtra); library(RColorBrewer); library(tidyverse);

file.path <- "data/"
# data from te_features_tally.R
list.file.v <- list.files(path = file.path,
                          pattern = ".fa_tally_filt_TRUE_FALSE$")

species.list <- gsub("_.*","",list.file.v)

df.complet <- data.frame()
for(species in 1:length(species.list)){
  species.i <- species.list[species]
  data.i <- read.delim(paste0(file.path,list.file.v[grepl(species.i,list.file.v)]))
  data.i$species <- species.i
  data.i <- data.i[which(data.i$feature != "Sequence_Ontology"),]
  df.complet <- rbind(df.complet,data.i)
}

df.complet$cls.0[which(df.complet$feature %in% "non_LTR_retrotransposon")] <- "nonLTR_retrotransposon"
df.complet$cls.0[which(df.complet$feature == "DIRS_YR_retrotransposon")] <- "Class_I_Retrotransposon"
df.complet$cls.0[which(df.complet$feature == "repeat_region")] <- "TE_unclass"
df.complet$feature[which(df.complet$feature == "repeat_region")] <- "TE_unclass"


df.complet$feature.raw <- df.complet$feature
df.complet$feature <- mapply(function(x) case_when(length(grep("_TIR_",x)) > 0 ~ "TIR_transposon",
                                                   length(grep("olinton",x)) > 0 ~ "Polinton",
                                                   length(grep("Tyrosine_Recombinase_Elements",x)) > 0 ~ "Recombinase_Elements",
                                                   length(grep("DIRS_YR_retrotransposon",x)) > 0 ~ "Recombinase_Elements",
                                                   .default = x),df.complet$feature.raw)

feature.order <- unique(df.complet$feature)[c(4,14,1,12,13,11,5,17,9,16,18,3,15,10,6,7,2,8,19)]

df.complet <- df.complet[which(df.complet$cls.0 != "rRNA_gene"),]
df.complet <- df.complet[which(df.complet$TRUE_FALSE.red == "IN"),]

df.complet %>% str()

df.complet <- df.complet %>%
  group_by(species) %>%
  mutate(cent.TE.space = sum(nt.sum)) %>%
  ungroup() %>%
  mutate(prop.nt.sum.cent = nt.sum/cent.TE.space)

# here the data set get reduced as per those features that include multiple levels as TIR transposons
df.complet$prop.nt.sum.cent.raw <- df.complet$prop.nt.sum.cent
df.complet <- df.complet %>%
  group_by(species,cls.0,feature) %>%
  summarise(prop.nt.sum.cent = sum(prop.nt.sum.cent.raw),
            total.nt.sum.cent = sum(nt.sum)) %>%
  ungroup()

df.complet$feature %>% unique() %>% length()
df.complet$feature[which(df.complet$prop.nt.sum.cent > 0.01)] %>% unique() %>% length()
df.complet$feature[which(df.complet$prop.nt.sum.cent > 0.05)] %>% unique() %>% length()

species_order_list <- read.table("df_order.txt", quote="\"", comment.char="")

df.complet$species <- df.complet$species %>%
  {gsub("\\..*","",.)}

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
}


df.complet$species <- df.complet$species %>%
  {factor(.,levels = species_order_list$V1)}
df.complet <- df.complet[with(df.complet,order(species,cls.0,feature)),]

cls.order <- df.complet %>%
  group_by(cls.0) %>%
  summarise(total.nt.sum = sum(total.nt.sum.cent)) %>%
  arrange(desc(total.nt.sum)) %>%
  select(cls.0) %>%
  as.vector()


df.complet.all <- df.complet[which(is.na(df.complet$species) == F),]


# text.size <- 8
text.size <- 7.5
bar.width <- 0.425


# vertical layout
df.complet <- df.complet.all[,]
df.complet <- left_join(df.complet,distinct(cls.0.prop.order[,c(1)]), by = "species")

max.total.nt.sum.cent <- df.complet %>%
  group_by(species) %>%
  summarise(total.nt.sum.cent = sum(total.nt.sum.cent))
max.total.nt.sum.cent <- max(max.total.nt.sum.cent$total.nt.sum.cent)

df.complet$group.0 <- df.complet$group.0 %>%
  {factor(.,levels = c("invertebrate","chordate","plant"))}

myCol.dark <-c("LTR_retrotransposon" = "#3E88AE", "Class_I_Retrotransposon" = "orange", "nonLTR_retrotransposon" = "#DD5571", "Class_II_DNA_Transposon" = "#F8D36B", "repeat_region" = "red", "rRNA_gene" = "cyan", "TE_unclass" = "grey70")
df.complet$cls.0 <- df.complet$cls.0 %>%
  {factor(.,levels = names(myCol.dark))}


# tag non relevant features
df.complet$feature.rel <- mapply(function(x) ifelse(x >= 0.01, "PASS", "NO_PASS"), df.complet$prop.nt.sum.cent)
df.complet$feature.rel %>%
  table()
df.complet$feature[which(df.complet$feature.rel == "PASS")] %>%
  table()

df.complet$feature <- df.complet$feature %>%
  {factor(.,levels = feature.order)}
df.complet <- df.complet[with(df.complet,order(species,feature)),]


labels.v <- c("TIR_DDE","LINE","Retrovirus","LTR_Gypsy","LTR_unk","TE_unclass","Penelope","LTR_Copia","Helitron","other nonLTR","BEL/Pao","YR","Polinton","SINE","Crypton")
names(labels.v) <- df.complet$feature[which(df.complet$feature.rel == "PASS")] %>%
  unique() %>%
  as.character()
df.complet$feature.0 <- df.complet$feature %>%
  as.character()
df.complet$feature <- NA_character_

for(ft in 1:length(labels.v)){
  ft.i <- names(labels.v)[ft]
  label.i <- labels.v[ft]
  df.complet$feature[which(df.complet$feature.0 == ft.i)] <- label.i
}
feature.order <- unique(df.complet$feature[which(df.complet$feature.rel == "PASS")])[c(9,13,1,15,12,7,2,10,14,3,11,8,4,5,6)]

# bubbles
df.complet.subset <- df.complet[which(df.complet$feature.rel == "PASS"),]
df.complet.subset$feature <- df.complet.subset$feature %>%
  {factor(.,levels = feature.order)}
df.complet.subset <- df.complet.subset[with(df.complet.subset,order(species,cls.0,feature)),]


p.1 <- ggplot(df.complet.subset) +
  geom_point(aes(x = feature, y = fct_inorder(species), fill = cls.0, size = prop.nt.sum.cent), shape = 21, color = "black",stroke = .35) +
  facet_grid(group.0 ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = myCol.dark, labels = c("LTR_retrotransposon" = "Class I LTR", "nonLTR_retrotransposon" = "Class I non-LTR TPRT", "Class_I_Retrotransposon" = "Class I non-LTR other", "Class_II_DNA_Transposon" = "Class II DNA element", "TE_unclass" = "TE unclass")) +
  scale_size_area(max_size = 10, n.breaks = 6) +
  scale_x_discrete(expand = expansion(mult = 0.05), breaks = feature.order, position = "top") + 
  guides(fill = guide_legend(override.aes = list(
    color = NA, # Make square transparent
    byrow = TRUE
  ))) +
  coord_cartesian(expand = T,clip = "off") +
  theme_minimal() +
  theme(legend.position = "none",
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.key.size = unit(.35,"cm"),
        legend.spacing = unit(2,"cm"),
        legend.text = element_text(size = text.size),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "grey95"),
        plot.margin = unit(c(t = 5,r = 5,b = 5,l = 5),"mm"),
        axis.text.x.top = element_text(angle = 45,hjust = 0),
        axis.title = element_blank(),
        axis.text.y = element_text(),
        strip.text = element_blank(),
        text = element_text(size = text.size))


pdf(paste0("p1_TRUE_FALSE_cent_space_prop_cutoff.all_in.palette.combined.bubbles.pdf"),onefile = T, width = 6,height = 7.15)
p.1
dev.off()


scf_name <- mapply(function (Genus,Species) paste0(gsub("(^[A-Za-z]{1}).*","\\1",Genus),". ", Species), data.attributes$Genus, data.attributes$Species)
data.attributes$scf_name <- scf_name

df.complet.subset <- left_join(df.complet.subset,data.attributes[,c(9,10)],by = "species")
df.complet.subset$species <- df.complet.subset$species %>%
  {factor(.,levels = species_order_list$V1)}
df.complet.subset <- df.complet.subset[with(df.complet.subset,order(species)),]


p.1 <- ggplot(df.complet.subset) +
  geom_point(aes(x = feature, y = fct_inorder(scf_name), fill = cls.0, size = prop.nt.sum.cent), shape = 21, color = "black",stroke = .25) +
  facet_grid(group.0 ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = myCol.dark, labels = c("LTR_retrotransposon" = "Class I LTR", "nonLTR_retrotransposon" = "Class I non-LTR TPRT", "Class_I_Retrotransposon" = "Class I non-LTR other", "Class_II_DNA_Transposon" = "Class II DNA element", "TE_unclass" = "TE unclass")) +
  scale_size_area(max_size = 9, n.breaks = 6) +
  scale_x_discrete(expand = expansion(mult = 0.05), breaks = feature.order, position = "top") + 
  guides(fill = guide_legend(override.aes = list(
    color = NA, # Make square transparent
    byrow = TRUE
  ))) +
  coord_cartesian(expand = T,clip = "off") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.title = element_blank(),
        legend.key.size = unit(.4,"cm"),
        legend.spacing = unit(2,"cm"),
        legend.text = element_text(size = text.size),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = "grey95"),
        panel.background = element_rect(fill = "grey97", color = NA),
        panel.spacing.y = unit(c(t = 4.15,b = 4.15),"mm"),
        plot.margin = unit(c(t = 5,r = 5,b = 5,l = 5),"mm"),
        axis.text.x.top = element_text(angle = 45,hjust = 0),
        axis.title = element_blank(),
        axis.text.y = element_text(face = "italic"),
        strip.text = element_blank(),
        text = element_text(size = text.size))


pdf(paste0("p1_TRUE_FALSE_cent_space_prop_cutoff.all_in.palette.combined.bubbles.scf_name.pdf"),onefile = T, width = 6.65,height = 7.15)
p.1
dev.off()
