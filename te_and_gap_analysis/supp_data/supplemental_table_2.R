
file.path <- "data/"

# data from cent_prop_combined.sh
list.files.v <- list.files(file.path, pattern = "_centromere_nt_prop_combined")

list.files.v <- list.files.v[!grepl("TRUE_NA",list.files.v)]
species.list <- list.files.v %>%
  {gsub("\\..*","",.)} %>%
  {unique(.)}


TE_general_cls <- read.table("TE_general_cls.txt", quote="\"", comment.char="")
colnames(TE_general_cls) <- c("cls.0","EDTA.feature")

gap <- 2
gap.i <- c("TRUE_TRUE","TRUE_FALSE","FALSE_TRUE","FALSE_FALSE")[gap]

non.filter.data <- data.frame()
for(species in 1:length(species.list)){
  species.i <- species.list[species]
  list.files.species <- list.files.v[grepl(species.i,list.files.v)]
  
  print(species.i)
  
  cent.prop.raw <- read.delim(paste0(file.path,"/",list.files.species), header=FALSE, comment.char = "")
  colnames(cent.prop.raw) <-c("nt.count","filter_edge","min.gap","side","feature","TE.cls","step")
  
  
  if(length(grep("TRUE_NA",cent.prop.raw$filter_edge)) == 0){
    
    # first condition: 0 && upper
    min.gap.length.v <- c(0,50,100,150,200,250)
    min.gap <- 1
    min.gap.length <- min.gap.length.v[min.gap]
    
    side <- 1
    side.i <- c("upper","lower")[side]
    
    cent.prop.i <- cent.prop.raw[which(cent.prop.raw$filter_edge == gap.i & cent.prop.raw$side == side.i & cent.prop.raw$min.gap == min.gap.length),]
    cent.prop.i$feature[which(cent.prop.i$step == 4 & cent.prop.i$feature != "gap_nt")] <- "EDTA_nt"
    
    # panel 3 combined
    data.i.sum <- cent.prop.i[c(which(cent.prop.i$step == 4 & cent.prop.i$feature == "EDTA_nt"),
                                which(cent.prop.i$step == 2 & cent.prop.i$feature == "TE_nt"),
                                which(cent.prop.i$step == 3 & cent.prop.i$feature == "TE_nt"),
                                which(cent.prop.i$step == 1)
                                ),]
    data.i.sum$TE.cls[which(data.i.sum$step == 2 & data.i.sum$feature == "TE_nt")] <- "all_TEs"
    data.i.sum$TE.cls[which(data.i.sum$step == 1 & data.i.sum$feature == "cent_nt")] <- "centromere"
    data.i.sum$TE.cls[which(data.i.sum$step == 1 & data.i.sum$feature == "gap_nt")] <- "gaps"
    
    table(data.i.sum$TE.cls)
    # infer TRASH contribution
    data.i.sum.TRASH <- data.frame(nt.count = data.i.sum$nt.count[which(data.i.sum$step == 1 & data.i.sum$feature == "cent_nt")] - data.i.sum$nt.count[which(data.i.sum$step == 1 & data.i.sum$feature == "gap_nt")],
                                   filter_edge = data.i.sum$filter_edge[which(data.i.sum$step == 1 & data.i.sum$feature == "cent_nt")],
                                   min.gap = data.i.sum$min.gap[which(data.i.sum$step == 1 & data.i.sum$feature == "cent_nt")],
                                   side = data.i.sum$side[which(data.i.sum$step == 1 & data.i.sum$feature == "cent_nt")],
                                   feature = "TR_nt",
                                   TE.cls = "TRASH",
                                   step = 1)
    
    data.i.sum <- rbind(data.i.sum,data.i.sum.TRASH)
    
    # infer unaccounted nt
    data.i.sum.unk <- data.frame(nt.count = data.i.sum$nt.count[which(data.i.sum$step == 1 & data.i.sum$feature == "gap_nt")] - data.i.sum$nt.count[which(data.i.sum$step == 2 & data.i.sum$feature == "TE_nt")],
                                 filter_edge = data.i.sum$filter_edge[which(data.i.sum$step == 2 & data.i.sum$feature == "TE_nt")],
                                 min.gap = data.i.sum$min.gap[which(data.i.sum$step == 2 & data.i.sum$feature == "TE_nt")],
                                 side = data.i.sum$side[which(data.i.sum$step == 2 & data.i.sum$feature == "TE_nt")],
                                 feature = "unk_nt",
                                 TE.cls = "unknown",
                                 step = 2)
    
    data.i.sum <- rbind(data.i.sum,data.i.sum.unk)
    
    data.i.sum <- data.i.sum %>%
      mutate(prop = round((nt.count/max(nt.count))*100,digits = 5),
             axis.label = ifelse(TE.cls != "centromere", paste0(TE.cls, ": ", str_pad(format(round((nt.count/max(nt.count))*100,digits = 2),nsmall = 2),width = 5,side = "left"),"%"), TE.cls),
             nt.count = nt.count)
    
    data.i.sum <- data.i.sum %>%
      mutate(species = species.i)
    
    non.filter.data <- rbind(non.filter.data,data.i.sum) 
  }
  
  rm(list = setdiff(ls(),c("species.list","list.files.v","file.path","TE_general_cls","gap","gap.i","non.filter.data","upper.filter.data","lower.filter.data")))
}


# combined version
non.filter.data <- non.filter.data[with(non.filter.data,order(prop)),]
non.filter.data$species <- non.filter.data$species %>%
  {factor(.,levels = unique(non.filter.data$species[which(non.filter.data$step == 1 & non.filter.data$TE.cls == "TRASH")]))}
non.filter.data <- non.filter.data[with(non.filter.data,order(species)),]

non.filter.data$species <- non.filter.data$species %>%
  {gsub("\\..*","",.)}

# next line can be uncommented in case we want to sort yaxis by taxonomic order
{
  data.attributes <- read.delim("dtol_plant_attributes.txt", header = T, sep = "\t")
  data.attributes <- data.attributes %>%
    mutate(species = gsub("\\..*","",fasta))
  
  data.attributes %>% str
  
  non.filter.data <- left_join(non.filter.data,data.attributes[,c(9,2,5)],by = "species")
  
  non.filter.data$group.0 <- mapply(function(x) case_when(x == "Dicots" ~ "plant",
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
                                                          .default = NA_character_),non.filter.data$Group)
}


non.filter.data$group.0 <- non.filter.data$group.0 %>%
  {factor(.,levels = c("invertebrate","chordate","plant"))}


# tally
distinct(non.filter.data[,c(10,13)]) %>% {table(.[,2])}

data.i <-non.filter.data
data.i <- data.i[,c(10:13,2,3,4,7,5,6,1,8)]
data.i <- data.i[with(data.i,order(species,step)),]

data.i <- data.i %>% 
  group_by(species) %>%
  mutate(prop.gap = (nt.count/nt.count[feature == "gap_nt"])*100,
         prop.te = (nt.count/nt.count[TE.cls == "all_TEs"])*100)

prop.gap <- mapply(function(x) ifelse(x > 100, NA, x),data.i$prop.gap)
data.i$prop.gap <- prop.gap
prop.te <- mapply(function(x) ifelse(x > 100, NA, x),data.i$prop.te)
data.i$prop.te <- prop.te

data.i$prop.te[which(data.i$step == 1)] <- NA
data.i$prop.gap[which(data.i$step == 1)] <- NA
data.i$prop.te[which(data.i$step == 2 & data.i$feature == "unk_nt")] <- NA


data.i <- data.i[with(data.i,order(species,step,feature)),]
write.table(data.i,"species_table_s2.txt",sep = "\t",quote = F,col.names = T,row.names = F)

# group.0
cent.space <- data.i[which(data.i$feature == "cent_nt"),] %>%
  group_by(group.0) %>%
  summarise(cent.space = sum(nt.count))

gap.space <- data.i[which(data.i$feature == "gap_nt"),] %>%
  group_by(group.0) %>%
  summarise(gap.space = sum(nt.count))

te.space <- data.i[which(data.i$TE.cls == "all_TEs"),] %>%
  group_by(group.0) %>%
  summarise(te.space = sum(nt.count))

# step 1
d1.sum <- data.i[which(data.i$step == 1 & data.i$feature != "cent_nt"),] %>%
  group_by(group.0,step,feature,TE.cls) %>%
  summarise(sum.nt.count = sum(nt.count),
            avg.prop = mean(prop),
            median.prop = median(prop)) 

d1.sum <- left_join(d1.sum,cent.space, by = "group.0")
d1.sum <- left_join(d1.sum,gap.space, by = "group.0")
d1.sum <- left_join(d1.sum,te.space, by = "group.0")

d1 <- d1.sum %>%
  mutate(sum.prop.cent = (sum.nt.count/cent.space)*100,
         sum.prop.gap = NA,
         sum.prop.te = NA) %>% 
  as.data.frame() %>%
  {.[,c(1,2,3,4,5,6,7,11,12,13)]}

# step 2
d1.sum <- data.i[which(data.i$step == 2),]  %>%
  group_by(group.0,step,feature,TE.cls) %>%
  summarise(sum.nt.count = sum(nt.count),
            avg.prop = mean(prop),
            median.prop = median(prop)) 

d1.sum <- left_join(d1.sum,cent.space, by = "group.0")
d1.sum <- left_join(d1.sum,gap.space, by = "group.0")
d1.sum <- left_join(d1.sum,te.space, by = "group.0")

d2 <- d1.sum %>%
  mutate(sum.prop.cent = (sum.nt.count/cent.space)*100,
         sum.prop.gap = (sum.nt.count/gap.space)*100,
         sum.prop.te = NA) %>% 
  as.data.frame() %>%
  {.[,c(1,2,3,4,5,6,7,11,12,13)]}

# step 3
d1.sum <- data.i[which(data.i$step == 3),] %>%
  group_by(group.0,step,feature,TE.cls) %>%
  summarise(sum.nt.count = sum(nt.count),
            avg.prop = mean(prop),
            median.prop = median(prop))

d1.sum <- left_join(d1.sum,cent.space, by = "group.0")
d1.sum <- left_join(d1.sum,gap.space, by = "group.0")
d1.sum <- left_join(d1.sum,te.space, by = "group.0")

d3 <- d1.sum %>%
  mutate(sum.prop.cent = (sum.nt.count/cent.space)*100,
         sum.prop.gap = (sum.nt.count/gap.space)*100,
         sum.prop.te = (sum.nt.count/te.space)*100) %>% 
  as.data.frame() %>%
  {.[,c(1,2,3,4,5,6,7,11,12,13)]}


# step 4
d1.sum <- data.i[which(data.i$step == 4),] %>%
  group_by(group.0,step,feature,TE.cls) %>%
  summarise(sum.nt.count = sum(nt.count),
            avg.prop = mean(prop),
            median.prop = median(prop))

d1.sum <- left_join(d1.sum,cent.space, by = "group.0")
d1.sum <- left_join(d1.sum,gap.space, by = "group.0")
d1.sum <- left_join(d1.sum,te.space, by = "group.0")

d4 <- d1.sum %>%
  mutate(sum.prop.cent = (sum.nt.count/cent.space)*100,
         sum.prop.gap = (sum.nt.count/gap.space)*100,
         sum.prop.te = (sum.nt.count/te.space)*100) %>% 
  as.data.frame() %>%
  {.[,c(1,2,3,4,5,6,7,11,12,13)]}

d.complete <- rbind(d1,d2,d3,d4)
write.table(d.complete,"taxa_table_s2.txt",col.names = T,row.names = F,quote = F)


# all in counts
data.i <- data.i %>%
  ungroup()

cent.space <- data.i[which(data.i$feature == "cent_nt"),] %>%
  summarise(cent.space = sum(nt.count))

gap.space <- data.i[which(data.i$feature == "gap_nt"),] %>%
  summarise(gap.space = sum(nt.count))

te.space <- data.i[which(data.i$TE.cls == "all_TEs"),] %>%
  summarise(te.space = sum(nt.count)) %>% as.vector()

d1 <- data.i[which(data.i$step == 1 & data.i$feature != "cent_nt"),] %>%
  group_by(step,feature,TE.cls) %>%
  summarise(sum.nt.count = sum(nt.count),
            avg.prop = mean(prop),
            median.prop = median(prop),
            cent.space.i = cent.space$cent.space) %>%
  mutate(sum.prop.cent = (sum.nt.count/cent.space.i)*100,
         sum.prop.gap = NA,
         sum.prop.te = NA) %>% 
  as.data.frame() %>%
  {.[,c(1,2,3,4,5,6,8,9,10)]}

d2 <- data.i[which(data.i$step == 2),] %>%
  group_by(step,feature,TE.cls) %>%
  summarise(sum.nt.count = sum(nt.count),
            avg.prop = mean(prop),
            median.prop = median(prop),
            cent.space.i = cent.space$cent.space,
            gap.space.i = gap.space$gap.space) %>%
  mutate(sum.prop.cent = (sum.nt.count/cent.space.i)*100,
         sum.prop.gap = (sum.nt.count/gap.space.i)*100,
         sum.prop.te = NA) %>% 
  as.data.frame() %>%
  {.[,c(1,2,3,4,5,6,9,10,11)]}


d3 <- data.i[which(data.i$step == 3),] %>%
  group_by(step,feature,TE.cls) %>%
  summarise(sum.nt.count = sum(nt.count),
            avg.prop = mean(prop),
            median.prop = median(prop),
            cent.space.i = cent.space$cent.space,
            gap.space.i = gap.space$gap.space,
            te.space.i = te.space$te.space) %>%
  mutate(sum.prop.cent = (sum.nt.count/cent.space.i)*100,
         sum.prop.gap = (sum.nt.count/gap.space.i)*100,
         sum.prop.te = (sum.nt.count/te.space.i)*100) %>% 
  as.data.frame() %>%
  {.[,c(1,2,3,4,5,6,10,11,12)]}


d4 <- data.i[which(data.i$step == 4),] %>%
  group_by(step,feature,TE.cls) %>%
  summarise(sum.nt.count = sum(nt.count),
            avg.prop = mean(prop),
            median.prop = median(prop),
            cent.space.i = cent.space$cent.space,
            gap.space.i = gap.space$gap.space,
            te.space.i = te.space$te.space) %>%
  mutate(sum.prop.cent = (sum.nt.count/cent.space.i)*100,
         sum.prop.gap = (sum.nt.count/gap.space.i)*100,
         sum.prop.te = (sum.nt.count/te.space.i)*100) %>% 
  as.data.frame() %>%
  {.[,c(1,2,3,4,5,6,10,11,12)]}

d.complete <- rbind(d1,d2,d3,d4)

write.table(d.complete,"all_table_s2.txt",col.names = T,row.names = F,quote = F)
