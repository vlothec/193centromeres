
options(scipen=999)

# data from drGeuUrba1_cent_prop_step_3.R
cent_feature_prop_by_chr <- read.delim("drGeuUrba1.1_cent_feature_prop_by_chr_clade_and_method_simplified.txt")
colnames(cent_feature_prop_by_chr)

cent_feature_prop_by_chr$chr <- cent_feature_prop_by_chr$chr %>%
  {factor(.,levels = unique(cent_feature_prop_by_chr$chr)[order(as.numeric(gsub("SUPER_","",unique(cent_feature_prop_by_chr$chr))))])}

cent_feature_prop_by_chr <- cent_feature_prop_by_chr[with(cent_feature_prop_by_chr,order(cent.feature.space)),]

cent_feature_prop_by_chr.sum <- cent_feature_prop_by_chr %>%
  group_by(chr,Clade.red) %>%
  summarise(total.cent.feature.space = sum(cent.feature.space))

cent_feature_prop_by_chr.sum$Clade.red <- cent_feature_prop_by_chr.sum$Clade.red %>%
  {factor(.,levels = c("other_TEs","Athila"))}
cent_feature_prop_by_chr.sum <- cent_feature_prop_by_chr.sum[with(cent_feature_prop_by_chr.sum,order(total.cent.feature.space,Clade.red)),]

# get PROT bit
{

# data from drGeuUrba1_cent_prop_step_5.R
cent_gap <- read.delim("drGeuUrba1.1_cent_gap_coord_all_TEs.txt")
cent_gap %>% str()

cent_gap.prot <- cent_gap[which(cent_gap$gap.length >= 625 & cent_gap$gap.length <= 631),] %>%
  group_by(chr) %>%
  summarise(sum.nt = sum(gap.length))

length_cent_chr <- cent_feature_prop_by_chr[,c(1,6)] %>% distinct()
cent_gap.prot <- left_join(length_cent_chr,cent_gap.prot, by = "chr")

cent_gap.prot <- cent_gap.prot %>%
  mutate(prot.prop = sum.nt/length.cent)
}

cent_feature_prop_by_chr.sum <- left_join(cent_feature_prop_by_chr.sum,cent_gap.prot[,c(1,4)], by = "chr")

cent_feature_prop_by_chr.sum <- cent_feature_prop_by_chr.sum %>%
  mutate(total.cent.feature.space.corrected = ifelse(Clade.red == "Athila",total.cent.feature.space+prot.prop,total.cent.feature.space))

cent_feature_prop_by_chr.sum$total.cent.feature.space.corrected[which(cent_feature_prop_by_chr.sum$chr == "SUPER_7" & cent_feature_prop_by_chr.sum$Clade.red == "Athila")] <- cent_feature_prop_by_chr.sum$total.cent.feature.space.corrected[which(cent_feature_prop_by_chr.sum$chr == "SUPER_7" & cent_feature_prop_by_chr.sum$Clade.red == "Athila")]-0.005

pdf("figure_s13_panel_A.pdf", width = 7.15, height = 5.5, onefile = T)

ggplot(cent_feature_prop_by_chr.sum) +
  geom_col(data = data.frame("chr" = unique(cent_feature_prop_by_chr.sum$chr)[order(as.numeric(gsub("SUPER_","",unique(cent_feature_prop_by_chr.sum$chr))))], "total.cent.feature.space.corrected" = 1),aes(y = chr, x = total.cent.feature.space.corrected), fill = "grey92", color = "black", width = 0.5) +
  geom_col(aes(y = chr, x = total.cent.feature.space.corrected, fill = Clade.red), color = "black", width = 0.5) +
  geom_col(data = cent_feature_prop_by_chr[which(cent_feature_prop_by_chr$Clade.red == "Athila" & cent_feature_prop_by_chr$method.ann == "structural"),],aes(y = chr, x = cent.feature.space), fill = NA, color = "black", width = 0.5) +
  scale_x_continuous(expand = expansion(0.01)) +
  scale_fill_manual(values = c("Athila" = "orange", "other_TEs" = "grey50")) +
  coord_cartesian(expand = T) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.line.x = element_line(),
        axis.ticks.x = element_line(),
        legend.key.size = unit(.85,"line"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.justification = 1,
        legend.title = element_blank(),
        panel.grid = element_blank())

dev.off()
