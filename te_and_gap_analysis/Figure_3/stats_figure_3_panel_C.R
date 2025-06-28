
file.path <- "data/"
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

taxa.groups <- df.complet[,c(9,12)] %>%
  distinct()

df.complet.red <- complete(df.complet[which(df.complet$ft.method == "structural"),c(9,2,6)],species,TRUE_FALSE.red,fill = list(value = 0))
df.complet.red$prop.nt.sum[which(is.na(df.complet.red$prop.nt.sum) == TRUE)] <- 0
df.complet.red <- left_join(df.complet.red,taxa.groups,by = "species")

# following lines are intended to assess that variance differs
{df.complet.red %>%
  group_by(TRUE_FALSE.red) %>%
  summarise(avg.prop.nt.sum = mean(prop.nt.sum),
            median.prop.nt.sum = median(prop.nt.sum),
            st.prop.nt.sum = sd(prop.nt.sum))

boxplot(prop.nt.sum ~ TRUE_FALSE.red,
        data = df.complet.red)

var.test(prop.nt.sum ~ TRUE_FALSE.red,
                 data = df.complet.red)
}

sink("stats_figure_3_panel_C_all.txt", split = TRUE)

print("all")

t.i <- t.test(df.complet.red$prop.nt.sum ~ df.complet.red$TRUE_FALSE.red,
       paired = FALSE)
print(t.i)

rm(t.i)
for(tax in 1:3){
  tax.i <- unique(df.complet$group.0)[tax]
  df.complet.red.subset <- df.complet.red[which(df.complet.red$group.0 == tax.i),] 
  print(tax.i)
  t.i <- t.test(df.complet.red.subset$prop.nt.sum ~ df.complet.red.subset$TRUE_FALSE.red,
         paired = FALSE)
  print(t.i)
  rm(t.i)
}

sink()
