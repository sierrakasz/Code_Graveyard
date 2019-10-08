rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/Forensic_Pig_Combined/")

#packages
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(multcompView)
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

physeq=merge_phyloseq(physeq, sampdat)

erich <- estimate_richness(physeq, measures = c("Observed"))
erich <- add_rownames(erich, "SampleID")
erich <- erich %>%
  gather(Index, Observation, c("Observed"), na.rm = TRUE)
write.table(erich, "alpha_div_forensicpig.tsv")
rich <- read.csv("metadata_alpha_div.csv")
rich <- rich %>% select(id, Index, Observation, Study,
                        Decomp_Stage)

rich$Index <- factor(rich$Index, levels = c("Observed Richness", "Shannon", "Faith PD"))
rich$Decomp_Stage <- factor(rich$Decomp_Stage, levels = c('Submerged_Fresh', 'Early_Floating',
                                                                  'Floating_Decay', 'Advanced_Floating_Decay', 
                                                                  'Sunken_Remains'))

rich_benbow <- subset(rich, Study == "Benbow")
rich_wallace <- subset(rich, Study == "Wallace")

theme_set(theme_bw(base_size = 30))
tiff("alphadiv_FP_benbow.TIF", width = 1200, height = 1200)
p <- ggplot(rich_benbow, aes(x=Decomp_Stage, y=Observation, fill=Decomp_Stage)) +
  geom_boxplot() + xlab("Decomposition Stage") +ylab("Alpha Diversity Measure")
p + facet_wrap(~Index, scales="free") + theme(axis.text.x = element_blank(),
                                              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                              panel.background = element_blank(),
                                              legend.position="none") +
  scale_fill_manual(values = c("#FF0000", "#006699", "#FF6600", "#FFCC00"),
                    name="Decomposition Stage",
                    labels=c('Submerged Fresh', 'Early Floating',
                             'Floating Decay', 'Advanced Floating Decay', 
                             'Sunken Remains'))
dev.off()

tiff("alphadiv_FP_wallace.TIF", width = 1200, height = 1200)
p <- ggplot(rich_wallace, aes(x=Decomp_Stage, y=Observation, fill=Decomp_Stage)) +
  geom_boxplot() + xlab("Decomposition Stage") +ylab("Alpha Diversity Measure")
p + facet_wrap(~Index, scales="free") + theme(axis.text.x = element_blank(),
                                              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                              panel.background = element_blank(),
                                              legend.position="none") +
  scale_fill_manual(values = c("#FF0000", "#006699", "#FF6600", "#FFCC00", "#00A08A"), 
                    name="Decomposition Stage",
                    labels=c('Submerged Fresh', 'Early Floating',
                             'Floating Decay', 'Advanced Floating Decay', 
                             'Sunken Remains'))

dev.off()

#legend dummy plot
tiff("legend_dummy.TIF", width = 1600, height = 1200)
dummy <- ggplot(rich_wallace, aes(x=Decomp_Stage, y=Observation, color=Decomp_Stage)) +
  geom_point(size = 8, shape = 15) +
  scale_color_manual(values = c("#FF0000", "#006699", "#FF6600", "#FFCC00", "#00A08A"), 
                    name="Decomposition Stage",
                    labels=c('Submerged Fresh', 'Early Floating',
                             'Floating Decay', 'Advanced Floating Decay', 
                             'Sunken Remains')) +
  theme(legend.position = 'bottom')
dummy
dev.off()

tiff("alphadiv_FP_benbow_day.TIF", width = 1000, height = 800)
pl <- ggplot(rich_benbow, aes(x=Day, y=Observation, fill=Day)) +
  geom_boxplot()
pl + facet_wrap(~Index, scales="free") + theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = c("#66023C", "#CB2314", "#273046", "#354823", "#E1BD6D", "#A9A9A9"))
dev.off()

tiff("alphadiv_FP_wallace_day.TIF", width = 1000, height = 800)
pl <- ggplot(rich_wallace, aes(x=Day, y=Observation, fill=Day)) +
  geom_boxplot()
pl + facet_wrap(~Index, scales="free") + theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values = c("#66023C", "#CB2314", "#273046", "#CC5500",
                               "#354823", "#E1BD6D", "#A9A9A9"))
dev.off()

#stats
prepare_samples_kw <- function(rich) {
  return(list(rich_obs <- rich %>% filter(Index == 'Observed'),
  rich_sha <- rich %>% filter(Index == 'Shannon'),
  rich_inv <- rich %>% filter(Index == 'Faith_PD')))
}

kw_values_benbow <- prepare_samples_kw(rich_benbow)
kw_values_wallace <- prepare_samples_kw(rich_wallace)


bm <- compare_means(Observation~Decomp_Stage, data=rich_benbow, group.by = "Index", method = "kruskal.test",p.adjust.method="bonferroni")
write.csv(bm, "benbow_alpha_kw_stats.csv")
Means=compare_means(Observation~Decomp_Stage, data = rich_benbow, 
                    group.by = "Index", p.adjust.method = "bonferroni")

NComparisons<-length(unique(rich$Decomp_Stage))*length(unique(rich$Index))
SigList<-length(unique(rich$Index))
SigLetters2<-vector(length=NComparisons)
#vec<-unlist(lst)
sink("benbow_alpha_wil_stats.csv")
for (i in levels(Means$Index)){
  Tax<-i
  TaxAbundance<-subset(Means,Index==i )
  print(TaxAbundance)
  Hyphenated<-as.character(paste0(TaxAbundance$group1,"-",TaxAbundance$group2))
  difference<-TaxAbundance$p.adj
  names(difference)<-Hyphenated
  Letters<-multcompLetters(difference)
  print(Letters)
  
}
sink()

bm <- compare_means(Observation~Decomp_Stage, data=rich_wallace, group.by = "Index", method = "kruskal.test",p.adjust.method="bonferroni")
write.csv(bm, "wallace_alpha_kw_stats.csv")
Means=compare_means(Observation~Decomp_Stage, data = rich_wallace, 
                    group.by = "Index", p.adjust.method = "bonferroni")

NComparisons<-length(unique(rich$Decomp_Stage))*length(unique(rich$Index))
SigList<-length(unique(rich$Index))
SigLetters2<-vector(length=NComparisons)
#vec<-unlist(lst)
sink("wallace_alpha_wil_stats.csv")
for (i in levels(Means$Index)){
  Tax<-i
  TaxAbundance<-subset(Means,Index==i )
  print(TaxAbundance)
  Hyphenated<-as.character(paste0(TaxAbundance$group1,"-",TaxAbundance$group2))
  difference<-TaxAbundance$p.adj
  names(difference)<-Hyphenated
  Letters<-multcompLetters(difference)
  print(Letters)
  
}
sink()

for(i in 1:length(kw_values_benbow)) {
  print(kruskal.test(Observation ~ Decomp_Stage, data = kw_values_benbow[[i]]))
  out <- posthoc.kruskal.nemenyi.test(x=kw_values_benbow[[i]]$Observation, g=kw_values_benbow[[i]]$Decomp_Stage, dist='Tukey', p.adjust.method = 'bonf' )
  print(out$p.value)
}

for(i in 1:length(kw_values_benbow)) {
  print(kruskal.test(Observation ~ Day, data = kw_values_benbow[[i]]))
  out <- posthoc.kruskal.nemenyi.test(x=kw_values_benbow[[i]]$Observation, g=kw_values_benbow[[i]]$Day, dist='Tukey', p.adjust.method = 'bonf' )
  print(out$p.value)
}

for(i in 1:length(kw_values_wallace)) {
  print(kruskal.test(Observation ~ Decomp_Stage, data = kw_values_wallace[[i]]))
  out <- posthoc.kruskal.nemenyi.test(x=kw_values_wallace[[i]]$Observation, g=kw_values_wallace[[i]]$Decomp_Stage, dist='Tukey', p.adjust.method = 'bonf' )
  print(out$p.value)
}

for(i in 1:length(kw_values_wallace)) {
  print(kruskal.test(Observation ~ Day, data = kw_values_wallace[[i]]))
  out <- posthoc.kruskal.nemenyi.test(x=kw_values_wallace[[i]]$Observation, g=kw_values_wallace[[i]]$Day, dist='Tukey', p.adjust.method = 'bonf' )
  print(out$p.value)
}


