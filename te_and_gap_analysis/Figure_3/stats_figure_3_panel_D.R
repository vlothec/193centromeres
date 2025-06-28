
file.path <- "data/"
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




# extract obs.value
df.complet.all.obs.values <- df.complet.all %>%
  select(species,TRUE_FALSE.red,ft.method,group.0,feature.0,obs.values) %>%
  mutate(obs.values.list = lapply(strsplit(obs.values, ";", TRUE), as.numeric))

df.complet.all.obs.values <- do.call('rbind', do.call('Map', c(data.frame, df.complet.all.obs.values)))


taxa.groups <- df.complet[,c(14,17)] %>%
  distinct()

df.complet.red <- df.complet.all.obs.values %>%
  group_by(species,TRUE_FALSE.red) %>%
  summarise(mean.obs.values = mean(obs.values.list,na.rm = T))

df.complet.red <- left_join(df.complet.red,taxa.groups,by = "species")


# following lines are intended to assess that variance differs
{df.complet.red %>%
    group_by(TRUE_FALSE.red) %>%
    summarise(avg.mean.obs.values = mean(mean.obs.values),
              median.mean.obs.values = median(mean.obs.values),
              st.mean.obs.values = sd(mean.obs.values))
  
  boxplot(mean.obs.values ~ TRUE_FALSE.red,
          data = df.complet.red)
  
  var.test(mean.obs.values ~ TRUE_FALSE.red,
           data = df.complet.red)
}

sink("stats_figure_3_panel_D_all.txt", split = TRUE)

print("all")

t.i <- t.test(df.complet.red$mean.obs.values ~ df.complet.red$TRUE_FALSE.red,
       paired = FALSE)
print(t.i)
rm(t.i)

for(tax in 1:3){
  tax.i <- unique(df.complet$group.0)[tax]
  df.complet.red.subset <- df.complet.red[which(df.complet.red$group.0 == tax.i),] 
  print(tax.i)
  t.i <- t.test(df.complet.red.subset$mean.obs.values ~ df.complet.red.subset$TRUE_FALSE.red,
                paired = FALSE)
  print(t.i)
  rm(t.i)
}

sink()
