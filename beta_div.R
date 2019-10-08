rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/Forensic_Pig_Combined/")

#packages
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(gridExtra)
library(phyloseq)
library(plyr)
library(PMCMR)
library(tidyverse)
library(vegan)

#make a table in excel for rare data
tax_group <- read.csv("tax_group_fp_combo.csv")
tax_group <- tax_group %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group <- tax_group %>%  ungroup()

regroup_physeq_object <-function(table) {
  tax <- table %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
  tax <- as.matrix(tax)
  otu <- table %>% select(contains("P"))
  otu <- otu[,-1]
  OTU=otu_table(otu, taxa_are_rows=TRUE)
  TAX=tax_table(tax)
  taxa_names(TAX)=row.names(OTU)
  physeq_all=phyloseq(OTU,TAX)
  return(physeq_all)
}

physeq <- regroup_physeq_object(tax_group)

metadata=(read.csv("AquaticPigMetadataCombined.csv",header=TRUE))
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
sampdat$Decomp_Stage <- factor(sampdat$Decomp_Stage, levels = c('Submerged_Fresh', 'Early_Floating',
                                                          'Floating_Decay', 'Advanced_Floating_Decay', 
                                                          'Sunken_Remains'))
physeq=merge_phyloseq(physeq, sampdat)

physeq_benbow <- subset_samples(physeq, Study == "Benbow")
physeq_wallace <- subset_samples(physeq, Study == "Wallace")


  phy_SF <- subset_samples(physeq_benbow, Decomp_Stage == 'Submerged_Fresh')
  phy_EF <- subset_samples(physeq_benbow, Decomp_Stage == 'Early_Floating')
  phy_FD <- subset_samples(physeq_benbow, Decomp_Stage == 'Floating_Decay')
  phy_AFD <- subset_samples(physeq_benbow, Decomp_Stage == 'Advanced_Floating_Decay')

  otu <- data.frame(otu_table(phy_SF))
  otu= otu[rowSums(otu)!=0,] 
  beta <- vegdist(otu, method = "jaccard")
  mean(beta)
  sd(beta)
  
  otu <- data.frame(otu_table(phy_EF))
  otu= otu[rowSums(otu)!=0,] 
  beta <- vegdist(otu, method = "jaccard")
  mean(beta)
  sd(beta)
  
  otu <- data.frame(otu_table(phy_FD))
  otu= otu[rowSums(otu)!=0,] 
  beta <- vegdist(otu, method = "jaccard")
  mean(beta)
  sd(beta)
  
  otu <- data.frame(otu_table(phy_AFD))
  otu= otu[rowSums(otu)!=0,] 
  beta <- vegdist(otu, method = "jaccard")
  mean(beta)
  sd(beta)
  
  phy_SF <- subset_samples(physeq_wallace, Decomp_Stage == 'Submerged_Fresh')
  phy_EF <- subset_samples(physeq_wallace, Decomp_Stage == 'Early_Floating')
  phy_FD <- subset_samples(physeq_wallace, Decomp_Stage == 'Floating_Decay')
  phy_AFD <- subset_samples(physeq_wallace, Decomp_Stage == 'Advanced_Floating_Decay')
  phy_SR <- subset_samples(physeq_wallace, Decomp_Stage == 'Sunken_Remains')
  
  otu <- data.frame(otu_table(phy_SF))
  otu= otu[rowSums(otu)!=0,] 
  beta <- vegdist(otu, method = "jaccard")
  mean(beta)
  sd(beta)
  
  otu <- data.frame(otu_table(phy_EF))
  otu= otu[rowSums(otu)!=0,] 
  beta <- vegdist(otu, method = "jaccard")
  mean(beta)
  sd(beta)
  
  otu <- data.frame(otu_table(phy_FD))
  otu= otu[rowSums(otu)!=0,] 
  beta <- vegdist(otu, method = "jaccard")
  mean(beta)
  sd(beta)
  
  otu <- data.frame(otu_table(phy_AFD))
  otu= otu[rowSums(otu)!=0,] 
  beta <- vegdist(otu, method = "jaccard")
  mean(beta)
  sd(beta)
  
  otu <- data.frame(otu_table(phy_SR))
  otu= otu[rowSums(otu)!=0,] 
  beta <- vegdist(otu, method = "jaccard")
  mean(beta)
  sd(beta)
  
theme_set(theme_bw(base_size = 30))
ord = ordinate(physeq_benbow, method="PCoA", distance="(A+B-2*J)/(A+B-J)")
ordplot=plot_ordination(physeq_benbow, ord, color="Decomp_Stage")
ordplot
tiff("beta_jaccard_FP_benbow.TIF", width = 800, height = 800)
ordplot_b_tax <- ordplot + geom_point(size = 4) + 
  stat_ellipse(alpha = 0.4, geom = "polygon", aes(fill = Decomp_Stage)) +
  scale_fill_manual(values = c("#FF0000","#006699", "#FF6600", "#FFCC00")) +
  scale_color_manual(values = c("#FF0000","#006699", "#FF6600", "#FFCC00")) + 
  theme(legend.position="none")
dev.off()

ordw = ordinate(physeq_wallace, method="PCoA", distance="(A+B-2*J)/(A+B-J)")
ordplotw=plot_ordination(physeq, ordw, color="Decomp_Stage")
ordplotw 
tiff("beta_jaccard_FP_wallace.TIF", width = 800, height = 800)
ordplot_w_tax <- ordplotw + geom_point(size = 4) + 
  stat_ellipse(alpha = 0.4, geom = "polygon", aes(fill = Decomp_Stage)) +
  scale_fill_manual(values = c("#FF0000", "#006699", "#FF6600", "#FFCC00", "#00A08A"), 
                    name="Decomposition Stage",
                    labels=c('Submerged Fresh', 'Early Floating',
                             'Floating Decay', 'Advanced Floating Decay', 
                             'Sunken Remains')) +
  scale_color_manual(values = c("#FF0000", "#006699", "#FF6600", "#FFCC00", "#00A08A"))
dev.off()
ordplot_w_tax

ord = ordinate(physeq_benbow, method="PCoA", distance="(A+B-2*J)/(A+B-J)")
ordplot=plot_ordination(physeq_benbow, ord, color="Day")
ordplot
tiff("beta_jaccard_FP_benbow_day.TIF", width = 1000, height = 800)
ordplot + geom_point(size = 4) + 
  stat_ellipse(alpha = 0.4, geom = "polygon", aes(fill = Day)) +
  scale_fill_manual(values = c("#66023C", "#CB2314", "#273046", "#354823", "#E1BD6D", "#A9A9A9")) +
  scale_color_manual(values = c("#66023C", "#CB2314", "#273046", "#354823", "#E1BD6D", "#A9A9A9"))
dev.off()

ordw = ordinate(physeq_wallace, method="PCoA", distance="(A+B-2*J)/(A+B-J)")
ordplotw=plot_ordination(physeq, ordw, color="Day")
ordplotw
tiff("beta_jaccard_FP_wallace_day.TIF", width = 1000, height = 800)
ordplotw + geom_point(size = 4) + 
  stat_ellipse(alpha = 0.4, geom = "polygon", aes(fill = Day)) +
  scale_fill_manual(values = c("#66023C", "#CB2314", "#273046", "#CC5500",
                               "#354823", "#E1BD6D", "#A9A9A9")) +
  scale_color_manual(values = c("#66023C", "#CB2314", "#273046", "#CC5500",
                               "#354823", "#E1BD6D", "#A9A9A9"))
dev.off()


#stats

GPdist=phyloseq::distance(physeq, "(A+B-2*J)/(A+B-J)")
sampledf <- data.frame(sample_data(physeq))
adonis(GPdist ~ Study, data = sampledf)
GPdist=phyloseq::distance(physeq, "(A+B-2*J)/(A+B-J)")
sampledf <- data.frame(sample_data(physeq))
beta <- betadisper(GPdist, sampledf$Study)
permutest(beta)

adonis(GPdist ~ Study*Decomp_Stage, data = sampledf)


beta_diversity_calc_decomp <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "(A+B-2*J)/(A+B-J)")
  sampledf <- data.frame(sample_data(physeq))
  print(return(adonis(GPdist ~ Decomp_Stage, data = sampledf)))
}

beta_dispersion_calc_decomp <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "(A+B-2*J)/(A+B-J)")
  sampledf <- data.frame(sample_data(physeq))
  beta <- betadisper(GPdist, sampledf$Decomp_Stage)
  print(return(permutest(beta)))
}

beta_diversity_calc_day <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "(A+B-2*J)/(A+B-J)")
  sampledf <- data.frame(sample_data(physeq))
  print(return(adonis(GPdist ~ Day, data = sampledf)))
}

beta_dispersion_calc_day <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "(A+B-2*J)/(A+B-J)")
  sampledf <- data.frame(sample_data(physeq))
  beta <- betadisper(GPdist, sampledf$Day)
  print(return(permutest(beta)))
}

beta_diversity_calc_decomp(physeq_benbow)
beta_dispersion_calc_decomp(physeq_benbow)
beta_diversity_calc_day(physeq_benbow)
beta_dispersion_calc_day(physeq_benbow)
beta_diversity_calc_decomp(physeq_wallace)
beta_dispersion_calc_decomp(physeq_wallace)
beta_diversity_calc_day(physeq_wallace)
beta_dispersion_calc_day(physeq_wallace)

decomp_pairs_benbow <- function(physeq) {
  phy_SF <- subset_samples(physeq, Decomp_Stage == 'Submerged_Fresh')
  phy_EF <- subset_samples(physeq, Decomp_Stage == 'Early_Floating')
  phy_FD <- subset_samples(physeq, Decomp_Stage == 'Floating_Decay')
  phy_AFD <- subset_samples(physeq, Decomp_Stage == 'Advanced_Floating_Decay')
  return(list(phy_EF.FD <- merge_phyloseq(phy_EF, phy_FD), 
              phy_FD.AFD <- merge_phyloseq(phy_AFD, phy_FD),
              phy_EF.AFD <- merge_phyloseq(phy_EF, phy_AFD),
              phy_SF.EF <- merge_phyloseq(phy_SF, phy_EF),
              phy_SF.FD <- merge_phyloseq(phy_SF, phy_FD),
              phy_SF.AFD <- merge_phyloseq(phy_SF, phy_AFD)))
}

decomp_list_benbow <- decomp_pairs_benbow(physeq_benbow)

decomp_pairs_wallace <- function(physeq) {
  phy_SF <- subset_samples(physeq, Decomp_Stage == 'Submerged_Fresh')
  phy_EF <- subset_samples(physeq, Decomp_Stage == 'Early_Floating')
  phy_FD <- subset_samples(physeq, Decomp_Stage == 'Floating_Decay')
  phy_AFD <- subset_samples(physeq, Decomp_Stage == 'Advanced_Floating_Decay')
  phy_SR <- subset_samples(physeq, Decomp_Stage == 'Sunken_Remains')
  return(list(phy_EF.FD <- merge_phyloseq(phy_EF, phy_FD), 
              phy_FD.AFD <- merge_phyloseq(phy_AFD, phy_FD),
              phy_EF.AFD <- merge_phyloseq(phy_EF, phy_AFD),
              phy_SF.EF <- merge_phyloseq(phy_SF, phy_EF),
              phy_SF.FD <- merge_phyloseq(phy_SF, phy_FD),
              phy_SF.AFD <- merge_phyloseq(phy_SF, phy_AFD),
              phy_SF.SR <- merge_phyloseq(phy_SF, phy_SR),
              phy_EF.SR <- merge_phyloseq(phy_EF, phy_SR),
              phy_AFD.SR <- merge_phyloseq(phy_AFD, phy_SR),
              phy_FD.SR <- merge_phyloseq(phy_FD, phy_SR)))
}

decomp_list_wallace <- decomp_pairs_wallace(physeq_wallace)

for(i in 1:length(decomp_list_benbow)) {
  print(beta_diversity_calc_decomp(decomp_list_benbow[[i]]))
  print(beta_dispersion_calc_decomp(decomp_list_benbow[[i]]))
}

for(i in 1:length(decomp_list_wallace)) {
  print(beta_diversity_calc_decomp(decomp_list_wallace[[i]]))
  print(beta_dispersion_calc_decomp(decomp_list_wallace[[i]]))
}

day_pairs_benbow <- function(physeq) {
  physeq_1 <- subset_samples(physeq, Day == '1')
  physeq_4 <- subset_samples(physeq, Day == '4')
  physeq_8 <- subset_samples(physeq, Day == '8')
  physeq_12 <- subset_samples(physeq, Day == '12')
  physeq_16 <- subset_samples(physeq, Day == '16')
  physeq_21 <- subset_samples(physeq, Day == '21')
  return(list(physeq_1.4 <- merge_phyloseq(physeq_1, physeq_4),
              physeq_1.8 <- merge_phyloseq(physeq_1, physeq_8),
              physeq_1.12 <- merge_phyloseq(physeq_1, physeq_12),
              physeq_1.16 <- merge_phyloseq(physeq_1, physeq_16),
              physeq_1.21 <- merge_phyloseq(physeq_1, physeq_21),
              physeq_4.8 <- merge_phyloseq(physeq_8, physeq_4),
              physeq_4.12 <- merge_phyloseq(physeq_12, physeq_4),
              physeq_4.16 <- merge_phyloseq(physeq_16, physeq_4),
              physeq_4.21 <- merge_phyloseq(physeq_21, physeq_4),
              physeq_8.12 <- merge_phyloseq(physeq_8, physeq_12),
              physeq_8.16 <- merge_phyloseq(physeq_8, physeq_16),
              physeq_8.21 <- merge_phyloseq(physeq_8, physeq_21),
              physeq_12.16 <- merge_phyloseq(physeq_16, physeq_12),
              physeq_12.21 <- merge_phyloseq(physeq_21, physeq_12),
              physeq_16.21 <- merge_phyloseq(physeq_16, physeq_21)))
}

day_list_benbow <- day_pairs_benbow(physeq_benbow)

day_pairs_wallace <- function(physeq) {
  physeq_1 <- subset_samples(physeq, Day == '1')
  physeq_4 <- subset_samples(physeq, Day == '4')
  physeq_7 <- subset_samples(physeq, Day == '7')
  physeq_10 <- subset_samples(physeq, Day == '10')
  physeq_13 <- subset_samples(physeq, Day == '13')
  physeq_16 <- subset_samples(physeq, Day == '16')
  physeq_19 <- subset_samples(physeq, Day == '19')
  return(list(physeq_1.4 <- merge_phyloseq(physeq_1, physeq_4),
              physeq_1.7 <- merge_phyloseq(physeq_1, physeq_7),
              physeq_1.10 <- merge_phyloseq(physeq_1, physeq_10),
              physeq_1.13 <- merge_phyloseq(physeq_1, physeq_13),
              physeq_1.16 <- merge_phyloseq(physeq_1, physeq_16),
              physeq_1.19 <- merge_phyloseq(physeq_1, physeq_19),
              physeq_4.7 <- merge_phyloseq(physeq_7, physeq_4),
              physeq_4.10 <- merge_phyloseq(physeq_4, physeq_10),
              physeq_4.13 <- merge_phyloseq(physeq_13, physeq_4),
              physeq_4.16 <- merge_phyloseq(physeq_16, physeq_4),
              physeq_4.19 <- merge_phyloseq(physeq_19, physeq_4),
              physeq_7.10 <- merge_phyloseq(physeq_7, physeq_10),
              physeq_7.13 <- merge_phyloseq(physeq_7, physeq_13),
              physeq_7.16 <- merge_phyloseq(physeq_7, physeq_16),
              physeq_7.19 <- merge_phyloseq(physeq_7, physeq_19),
              physeq_10.13 <- merge_phyloseq(physeq_10, physeq_13),
              physeq_10.16 <- merge_phyloseq(physeq_10, physeq_16),
              physeq_10.19 <- merge_phyloseq(physeq_10, physeq_19),
              physeq_13.16 <- merge_phyloseq(physeq_16, physeq_13),
              physeq_13.19 <- merge_phyloseq(physeq_19, physeq_13),
              physeq_16.19 <- merge_phyloseq(physeq_16, physeq_19)))
}

day_list_wallace <- day_pairs_wallace(physeq_wallace)

for(i in 1:length(day_list_benbow)) {
  print(beta_diversity_calc_day(day_list_benbow[[i]]))
  print(beta_dispersion_calc_day(day_list_benbow[[i]]))
}

for(i in 1:length(day_list_wallace)) {
  print(beta_diversity_calc_day(day_list_wallace[[i]]))
  print(beta_dispersion_calc_day(day_list_wallace[[i]]))
}

#picrust

benbow_pi <- import_biom("table-with-taxonomy_benbow.biom")
benbow_pi <- merge_phyloseq(benbow_pi, sampdat)
wallace_pi <- import_biom("table-with-taxonomy_wallace.biom")
wallace_pi <- merge_phyloseq(wallace_pi, sampdat)

phy_SF <- subset_samples(benbow_pi, Decomp_Stage == 'Submerged_Fresh')
phy_EF <- subset_samples(benbow_pi, Decomp_Stage == 'Early_Floating')
phy_FD <- subset_samples(benbow_pi, Decomp_Stage == 'Floating_Decay')
phy_AFD <- subset_samples(benbow_pi, Decomp_Stage == 'Advanced_Floating_Decay')

otu <- data.frame(otu_table(phy_SF))
otu= otu[rowSums(otu)!=0,] 
beta <- vegdist(otu, method = "jaccard")
mean(beta)
sd(beta)

otu <- data.frame(otu_table(phy_EF))
otu= otu[rowSums(otu)!=0,] 
beta <- vegdist(otu, method = "jaccard")
mean(beta)
sd(beta)

otu <- data.frame(otu_table(phy_FD))
otu= otu[rowSums(otu)!=0,] 
beta <- vegdist(otu, method = "jaccard")
mean(beta)
sd(beta)

otu <- data.frame(otu_table(phy_AFD))
otu= otu[rowSums(otu)!=0,] 
beta <- vegdist(otu, method = "jaccard")
mean(beta)
sd(beta)

phy_SF <- subset_samples(wallace_pi, Decomp_Stage == 'Submerged_Fresh')
phy_EF <- subset_samples(wallace_pi, Decomp_Stage == 'Early_Floating')
phy_FD <- subset_samples(wallace_pi, Decomp_Stage == 'Floating_Decay')
phy_AFD <- subset_samples(wallace_pi, Decomp_Stage == 'Advanced_Floating_Decay')
phy_SR <- subset_samples(wallace_pi, Decomp_Stage == 'Sunken_Remains')

otu <- data.frame(otu_table(phy_SF))
otu= otu[rowSums(otu)!=0,] 
beta <- vegdist(otu, method = "jaccard")
mean(beta)
sd(beta)

otu <- data.frame(otu_table(phy_EF))
otu= otu[rowSums(otu)!=0,] 
beta <- vegdist(otu, method = "jaccard")
mean(beta)
sd(beta)

otu <- data.frame(otu_table(phy_FD))
otu= otu[rowSums(otu)!=0,] 
beta <- vegdist(otu, method = "jaccard")
mean(beta)
sd(beta)

otu <- data.frame(otu_table(phy_AFD))
otu= otu[rowSums(otu)!=0,] 
beta <- vegdist(otu, method = "jaccard")
mean(beta)
sd(beta)

otu <- data.frame(otu_table(phy_SR))
otu= otu[rowSums(otu)!=0,] 
beta <- vegdist(otu, method = "jaccard")
mean(beta)
sd(beta)

theme_set(theme_bw(base_size = 18))
ord = ordinate(benbow_pi, method="PCoA", distance="(A+B-2*J)/(A+B-J)")
ordplot_pi_b=plot_ordination(benbow_pi, ord, color="Decomp_Stage")
ordplot_pi_b  + theme(legend.position="none")
tiff("beta_jaccard_FP_benbow_picr.TIF", width = 800, height = 800)
ordplot_b_func<- ordplot_pi_b + geom_point(size = 4) + 
  stat_ellipse(alpha = 0.4, geom = "polygon", aes(fill = Decomp_Stage)) +
  scale_fill_manual(values = c("#FF0000","#006699", "#FF6600", "#FFCC00")) +
  scale_color_manual(values = c("#FF0000","#006699", "#FF6600", "#FFCC00")) + 
  theme(legend.position="none")
dev.off()

ordw = ordinate(wallace_pi, method="PCoA", distance="(A+B-2*J)/(A+B-J)")
ordplot_pi_W=plot_ordination(wallace_pi, ordw, color="Decomp_Stage")
ordplot_pi_W + theme(legend.position="none")
tiff("beta_jaccard_FP_wallace_picrus.TIF", width = 800, height = 800)
ordplot_w_func <- ordplot_pi_W + geom_point(size = 4) + 
  stat_ellipse(alpha = 0.4, geom = "polygon", aes(fill = Decomp_Stage)) +
  scale_fill_manual(values = c("#FF0000", "#006699", "#FF6600", "#FFCC00", "#00A08A")) +
  scale_color_manual(values = c("#FF0000", "#006699", "#FF6600", "#FFCC00", "#00A08A")) + 
  theme(legend.position="none") 
dev.off()

beta_diversity_calc_decomp(benbow_pi)
beta_dispersion_calc_decomp(benbow_pi)
beta_diversity_calc_decomp(wallace_pi)
beta_dispersion_calc_decomp(wallace_pi)


decomp_list_benbow <- decomp_pairs_benbow(benbow_pi)

for(i in 1:length(decomp_list_benbow)) {
  print(beta_diversity_calc_decomp(decomp_list_benbow[[i]]))
  print(beta_dispersion_calc_decomp(decomp_list_benbow[[i]]))
}

decomp_list_wallace <- decomp_pairs_wallace(wallace_pi)

for(i in 1:length(decomp_list_wallace)) {
  print(beta_diversity_calc_decomp(decomp_list_wallace[[i]]))
  print(beta_dispersion_calc_decomp(decomp_list_wallace[[i]]))
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(ordplot_w_tax)

tiff("beta_total.TIF", width = 1200, height = 1200)
grid.arrange(ordplot_b_tax, ordplot_w_tax + theme(legend.position="none"), ordplot_b_func, 
                         ordplot_w_func, nrow=2)
dev.off()




