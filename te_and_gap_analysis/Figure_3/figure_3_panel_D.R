
library(cowplot); library(ggnewscale); library(ggpubr); library(grid); library(gridExtra); library(RColorBrewer); library(tidyverse);

file.path <- "data/"

# data from te_pident_stats_intact.R
list.file.v <- list.files(path = file.path,
                          pattern = "_tally_filt_TRUE_FALSE_pident$")

species.list <- gsub("_.*","",list.file.v)


df.complet <- data.frame()
for(species in 1:length(species.list)){
  species.i <- species.list[species]
  data.i <- read.delim(paste0(file.path,list.file.v[grepl(species.i,list.file.v)]))
  data.i$species <- species.i
  df.complet <- rbind(df.complet,data.i)
}

df.complet %>% colnames()
df.complet %>% head()
df.complet$species <- df.complet$species %>%
  {gsub("\\..*","",.)}

species_order_list <- read.table("df_order.txt", quote="\"", comment.char="")

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

df.complet.all <- df.complet[which(df.complet$feature.0 != "TIR" & is.na(df.complet$species) == F),]

phylotaxa.count <- df.complet.all %>%
  select(species,group.0) %>%
  distinct() %>%
  summarise(n.species = n(), .by = "group.0")
phylotaxa.count <- phylotaxa.count %>%
  arrange(desc(n.species))
phylotaxa.count.max <- phylotaxa.count$n.species[which.max(phylotaxa.count$n.species)]

cls.0.prop.order <- df.complet.all %>%
  group_by(species,group.0) %>%
  select(species,group.0) %>%
  distinct() %>% data.frame %>% complete(species)

cls.0.prop.order <- cls.0.prop.order[with(cls.0.prop.order,order(cls.0.prop.order$species)),]
cls.0.prop.order <- cls.0.prop.order %>%
  group_by(group.0) %>%
  mutate(species.num = seq(1,n(),1))

# extract obs.value
df.complet.all.obs.values <- df.complet.all %>%
  select(species,TRUE_FALSE.red,ft.method,group.0,feature.0,obs.values) %>%
  mutate(obs.values.list = lapply(strsplit(obs.values, ";", TRUE), as.numeric))

df.complet.all.obs.values <- do.call('rbind', do.call('Map', c(data.frame, df.complet.all.obs.values)))

# vertical layout
# text.size <- 8
text.size <- 6.5
bar.width <- 0.425

df.complet <- df.complet.all[,]
df.complet <- left_join(df.complet,cls.0.prop.order[,c(1,3)], by = "species")

df.complet.obs.values <- df.complet.all.obs.values[,]
df.complet.obs.values <- left_join(df.complet.obs.values,cls.0.prop.order[,c(1,3)], by = "species")


r1.p1 <- ggplot(df.complet) +
  geom_vline(xintercept = 0.90, linetype = "dotted", color = "transparent") +
  geom_jitter(data = df.complet.obs.values[which(df.complet.obs.values$TRUE_FALSE.red == "OUT"),], aes(y = fct_inorder(species), x = obs.values.list, color = TRUE_FALSE.red), size = .1,width = 0,height = 0.2) + 
  geom_jitter(data = df.complet.obs.values[which(df.complet.obs.values$TRUE_FALSE.red != "OUT"),], aes(y = fct_inorder(species), x = obs.values.list, color = TRUE_FALSE.red), size = .1,width = 0,height = 0.2) + 
  scale_color_manual(values = c("OUT" = "#DDDDDD","IN" = "#3E88AE")) +
  scale_x_continuous(breaks = seq(0.85,1,0.05),labels = seq(0.85,1,0.05),expand = expansion(0,0)) +
  scale_y_discrete(expand = expansion(add = 0.85)) +
  facet_grid(group.0 ~ .,scales = "free_y",space = "free_y") +
  coord_cartesian(expand = T, xlim = c(0.835,1.015)) +
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

pdf(paste0("te_pident_intact_plot_col_LTR.pdf"),onefile = T, width = 2.75,height = 7.15)
print(r1.p1)
dev.off()


scf_name <- mapply(function (Genus,Species) paste0(gsub("(^[A-Za-z]{1}).*","\\1",Genus),". ", Species), data.attributes$Genus, data.attributes$Species)
data.attributes$scf_name <- scf_name

df.complet <- left_join(df.complet,data.attributes[,c(9,10)],by = "species")
df.complet$species <- df.complet$species %>%
  {factor(.,levels = species_order_list$V1)}
df.complet <- df.complet[with(df.complet,order(species)),]

df.complet.obs.values <- left_join(df.complet.obs.values,data.attributes[,c(9,10)],by = "species")
df.complet.obs.values$species <- df.complet.obs.values$species %>%
  {factor(.,levels = species_order_list$V1)}
df.complet.obs.values <- df.complet.obs.values[with(df.complet.obs.values,order(species)),]

r1.p1 <- ggplot(df.complet) +
  geom_vline(xintercept = 0.90, linetype = "dotted", color = "transparent") +
  geom_jitter(data = df.complet.obs.values[which(df.complet.obs.values$TRUE_FALSE.red == "OUT"),], aes(y = fct_inorder(scf_name), x = obs.values.list, color = TRUE_FALSE.red), size = .1,width = 0,height = 0.2) + 
  geom_jitter(data = df.complet.obs.values[which(df.complet.obs.values$TRUE_FALSE.red != "OUT"),], aes(y = fct_inorder(scf_name), x = obs.values.list, color = TRUE_FALSE.red), size = .1,width = 0,height = 0.2) + 
  scale_color_manual(values = c("OUT" = "#DDDDDD","IN" = "#3E88AE")) +
  scale_x_continuous(breaks = seq(0.85,1,0.05),labels = seq(0.85,1,0.05),expand = expansion(0,0)) +
  scale_y_discrete(expand = expansion(add = 0.85)) +
  facet_grid(group.0 ~ .,scales = "free_y",space = "free_y") +
  coord_cartesian(expand = T, xlim = c(0.835,1.015)) +
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
        axis.text.y = element_text(face = "italic"),
        strip.text = element_blank(),
        text = element_text(size = text.size))

pdf(paste0("te_pident_intact_plot_col_LTR.scf_name.pdf"),onefile = T, width = 2.75,height = 7.15)
print(r1.p1)
dev.off()
