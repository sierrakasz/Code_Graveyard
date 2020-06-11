rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/ThesisProject_AKP/")

#packages
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(phyloseq)
library(plyr)
library(PMCMR)
library(tidyverse)
library(vegan)

#upload the data
tax_group <- read.csv("tax_group_AKP.csv")
tax_group <- tax_group %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group <- tax_group %>%  ungroup()

regroup_physeq_object <-function(table) {
  tax <- table %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
  tax <- as.matrix(tax)
  otu <- table %>% select(contains("WCME"))
  OTU=otu_table(otu, taxa_are_rows=TRUE)
  TAX=tax_table(tax)
  taxa_names(TAX)=row.names(OTU)
  physeq_all=phyloseq(OTU,TAX)
  return(physeq_all)
}

physeq_all <- regroup_physeq_object(tax_group)

#make a table in excel for rare data
metadata=(read.csv("HPMMMetadata_pt2_cod.csv",header=TRUE))
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID

physeq=merge_phyloseq(physeq_all, sampdat)

physeq_rec <- subset_samples(physeq, Sample_Area == 'Rectum')
physeq_eye <- subset_samples(physeq, Sample_Area == 'Eyes')
physeq_ear <- subset_samples(physeq, Sample_Area == 'Ears')
physeq_mou <- subset_samples(physeq, Sample_Area == 'Nose')
physeq_nos <- subset_samples(physeq, Sample_Area == 'Mouth')


normalize_wout_rarefying <- function(physeq) {
  otu <- data.frame(otu_table(physeq))
  tax <- data.frame(tax_table(physeq))
  sample_nums <- rowSums(otu == 0)
  otu$OTUID <- row.names(otu)
  tax$OTUID <- row.names(tax)
  len <- length(colnames(otu)) -1
  core_perc <- sample_nums / len
  new_otu <- cbind(otu, core_perc)
  new_core_taxa <- new_otu$OTUID[new_otu$core_perc <= 0.99]
  new_core_otu <- filter(new_otu, new_otu$OTUID %in% new_core_taxa)
  new_core_otu <- subset(new_core_otu, select=-c(OTUID, core_perc))
  new_core_tax <- filter(tax, tax$OTUID %in% new_core_taxa)
  dis_core_tax <- filter(tax, !tax$OTUID %in% new_core_tax$OTUID)
  newphyseq <- prune_taxa(new_core_tax$OTUID, physeq)
  return(newphyseq)
}

physeq_rec_norm <- normalize_wout_rarefying(physeq_rec)
physeq_eye_norm <- normalize_wout_rarefying(physeq_eye)
physeq_ear_norm <- normalize_wout_rarefying(physeq_ear)
physeq_mou_norm <- normalize_wout_rarefying(physeq_mou)
physeq_nos_norm <- normalize_wout_rarefying(physeq_nos)

GPfr_rec_24 <-subset_samples(physeq_rec_norm, FinePMI == 'Less24')
GPfr_eye_24 <-subset_samples(physeq_eye_norm, FinePMI == 'Less24')
GPfr_ear_24 <-subset_samples(physeq_ear_norm, FinePMI == 'Less24')
GPfr_nos_24 <-subset_samples(physeq_nos_norm, FinePMI == 'Less24')
GPfr_mou_24 <-subset_samples(physeq_mou_norm, FinePMI == 'Less24')

lets_find_the_core <- function(physeq) {
  otu <- data.frame(otu_table(physeq))
  tax <- data.frame(tax_table(physeq))
  sample_nums <- rowSums(otu == 0)
  otu$OTUID <- row.names(otu)
  tax$OTUID <- row.names(tax)
  len <- length(colnames(otu)) -1
  core_perc <- sample_nums / len
  new_otu <- cbind(otu, core_perc)
  new_core_taxa <- new_otu$OTUID[new_otu$core_perc <= 0.10]
  new_core_otu <- filter(new_otu, new_otu$OTUID %in% new_core_taxa)
  new_core_otu <- subset(new_core_otu, select=-c(OTUID, core_perc))
  new_core_tax <- filter(tax, tax$OTUID %in% new_core_taxa)
  return(new_core_taxa)
}

core_list_rec <- lets_find_the_core(physeq_norm_rec)
core_list_eye <- lets_find_the_core(physeq_norm_eye)
core_list_ear <- lets_find_the_core(physeq_norm_ear)
core_list_nos <- lets_find_the_core(physeq_norm_nos)
core_list_mou <- lets_find_the_core(physeq_norm_mou)

proportion_of_core <-function(physeq, core_list) {
  otu <- data.frame(otu_table(physeq))
  tax <- data.frame(tax_table(physeq))
  otu$OTUID <- row.names(otu)
  tax$OTUID <- row.names(tax)
  new_core_otu <- filter(otu, otu$OTUID %in% core_list)
  new_core_otu[,'OTUID'] <- list(NULL)
  core_percent <- c()
  for(i in 1:length(new_core_otu)) {
    core_percent <- c(core_percent, sum(new_core_otu[,i])/sum(otu[,i]))
  }
  core_percent_df <- data.frame(colnames(new_core_otu), core_percent)
  colnames(core_percent_df) <- c('SampleID', 'Core_percent')
  return(core_percent_df)
}

core_percent_rec <- proportion_of_core(physeq_norm_rec, core_list_rec)
core_percent_eye <- proportion_of_core(physeq_norm_eye, core_list_eye)
core_percent_ear <- proportion_of_core(physeq_norm_ear, core_list_ear)
core_percent_mou <- proportion_of_core(physeq_norm_mou, core_list_mou)
core_percent_nos <- proportion_of_core(physeq_norm_nos, core_list_nos)


core_df <- rbind(core_percent_rec, core_percent_eye, core_percent_ear, core_percent_mou, core_percent_nos)

core_df = merge(core_df, metadata, by= 'SampleID')
core_df <- core_df %>% select(SampleID, Core_percent, MoD, CoD_Simple, Violent,
                          Sample_Area, FinePMI)

kruskal.test(Core_percent ~ Sample_Area, data = core_df)
out <- posthoc.kruskal.nemenyi.test(x=core_df$Core_percent, g=core_df$Sample_Area, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Core_percent ~ MoD, data = core_df)
out <- posthoc.kruskal.nemenyi.test(x=core_df$Core_percent, g=core_df$MoD, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Core_percent ~ CoD_Simple, data = core_df)
out <- posthoc.kruskal.nemenyi.test(x=core_df$Core_percent, g=core_df$CoD_Simple, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Core_percent ~ Violent, data = core_df)



core_df_rec = merge(core_percent_rec, metadata, by= 'SampleID')
core_df_rec <- core_df_rec %>% select(SampleID, Core_percent, MoD, CoD_Simple, Violent,
                              Sample_Area, FinePMI)


kruskal.test(Core_percent ~ MoD, data = core_df_rec)
out <- posthoc.kruskal.nemenyi.test(x=core_df_rec$Core_percent, g=core_df_rec$MoD, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Core_percent ~ CoD_Simple, data = core_df_rec)
out <- posthoc.kruskal.nemenyi.test(x=core_df_rec$Core_percent, g=core_df_rec$CoD_Simple, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Core_percent ~ Violent, data = core_df_rec)

core_df_eye = merge(core_percent_eye, metadata, by= 'SampleID')
core_df_eye <- core_df_eye %>% select(SampleID, Core_percent, MoD, CoD_Simple, Violent,
                                      Sample_Area, FinePMI)


kruskal.test(Core_percent ~ MoD, data = core_df_eye)
out <- posthoc.kruskal.nemenyi.test(x=core_df_eye$Core_percent, g=core_df_eye$MoD, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Core_percent ~ CoD_Simple, data = core_df_eye)
out <- posthoc.kruskal.nemenyi.test(x=core_df_eye$Core_percent, g=core_df_eye$CoD_Simple, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Core_percent ~ Violent, data = core_df_eye)

core_df_ear = merge(core_percent_ear, metadata, by= 'SampleID')
core_df_ear <- core_df_ear %>% select(SampleID, Core_percent, MoD, CoD_Simple, Violent,
                                      Sample_Area, FinePMI)


kruskal.test(Core_percent ~ MoD, data = core_df_ear)
out <- posthoc.kruskal.nemenyi.test(x=core_df_ear$Core_percent, g=core_df_ear$MoD, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Core_percent ~ CoD_Simple, data = core_df_ear)
out <- posthoc.kruskal.nemenyi.test(x=core_df_ear$Core_percent, g=core_df_ear$CoD_Simple, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Core_percent ~ Violent, data = core_df_ear)


core_df_nos = merge(core_percent_nos, metadata, by= 'SampleID')
core_df_nos <- core_df_nos %>% select(SampleID, Core_percent, MoD, CoD_Simple, Violent,
                                      Sample_Area, FinePMI)


kruskal.test(Core_percent ~ MoD, data = core_df_nos)
out <- posthoc.kruskal.nemenyi.test(x=core_df_nos$Core_percent, g=core_df_nos$MoD, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Core_percent ~ CoD_Simple, data = core_df_nos)
out <- posthoc.kruskal.nemenyi.test(x=core_df_nos$Core_percent, g=core_df_nos$CoD_Simple, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Core_percent ~ Violent, data = core_df_nos)

core_df_mou = merge(core_percent_mou, metadata, by= 'SampleID')
core_df_mou <- core_df_mou %>% select(SampleID, Core_percent, MoD, CoD_Simple, Violent,
                                      Sample_Area, FinePMI)


kruskal.test(Core_percent ~ MoD, data = core_df_mou)
out <- posthoc.kruskal.nemenyi.test(x=core_df_mou$Core_percent, g=core_df_mou$MoD, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Core_percent ~ CoD_Simple, data = core_df_mou)
out <- posthoc.kruskal.nemenyi.test(x=core_df_mou$Core_percent, g=core_df_mou$CoD_Simple, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Core_percent ~ Violent, data = core_df_mou)


core_df <- rbind(core_percent_eye, core_percent_ear, core_percent_nos)

core_df = merge(core_df, metadata, by= 'SampleID')
core_df <- core_df %>% select(SampleID, Core_percent, MoD, CoD_Simple, Violent,
                              Sample_Area, FinePMI)

kruskal.test(Core_percent ~ Sample_Area, data = core_df)
out <- posthoc.kruskal.nemenyi.test(x=core_df$Core_percent, g=core_df$Sample_Area, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Core_percent ~ MoD, data = core_df)
out <- posthoc.kruskal.nemenyi.test(x=core_df$Core_percent, g=core_df$MoD, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Core_percent ~ CoD_Simple, data = core_df)
out <- posthoc.kruskal.nemenyi.test(x=core_df$Core_percent, g=core_df$CoD_Simple, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Core_percent ~ Violent, data = core_df)


View(core_df %>% group_by(MoD) %>% summarise_all(funs(mean, sd)))
core_df %>% group_by(MoD) %>% summarise(no_rows = length(MoD))


theme_set(theme_classic(base_size = 20))
tiff("MoD_coreprop.TIF", width = 2000, height = 2000, res=300)
ggplot(core_df, aes(x=MoD, y=Core_percent, color= MoD)) + geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2)) + xlab('Manner of Death') +
  ylab('Core taxa sequences/ total sequences') + theme(legend.position = "none")
dev.off()






core_df <- merge(core_df, metadata, by = 'SampleID')
write.csv(core_df, 'core_percentage_norare.csv')

physeq_ear_norm_Core <- prune_taxa(core_list_ear, physeq_ear_norm)
physeq_eye_norm_Core <- prune_taxa(core_list_eye, physeq_eye_norm)
physeq_nos_norm_Core <- prune_taxa(core_list_nos, physeq_nos_norm)
physeq_mou_norm_Core <- prune_taxa(core_list_mou, physeq_mou_norm)
physeq_rec_norm_Core <- prune_taxa(core_list_rec, physeq_rec_norm)

noncore_physeq <- function(physeq, core_list) {
  otu <- data.frame(otu_table(physeq))
  tax <- data.frame(tax_table(physeq))
  otu$OTUID <- row.names(otu)
  tax$OTUID <- row.names(tax)
  dis_core_tax <- filter(tax, !tax$OTUID %in% core_list)
  dis_core_otu <- dis_core_tax$OTUID
  newphyseq <- prune_taxa(dis_core_otu, physeq)
  return(newphyseq)
}

physeq_ear_norm_non <- noncore_physeq(physeq_ear_norm, core_list_ear)
physeq_eye_norm_non <- noncore_physeq(physeq_eye_norm, core_list_eye)
physeq_rec_norm_non <- noncore_physeq(physeq_rec_norm, core_list_rec)
physeq_mou_norm_non <- noncore_physeq(physeq_mou_norm, core_list_mou)
physeq_nos_norm_non <- noncore_physeq(physeq_nos_norm, core_list_nos)

beta_plots_list <- list(physeq_ear_norm_Core, physeq_ear_norm_non,
                        physeq_eye_norm_Core, physeq_eye_norm_non,
                        physeq_rec_norm_Core, physeq_rec_norm_non,
                        physeq_mou_norm_Core, physeq_mou_norm_non,
                        physeq_nos_norm_Core, physeq_nos_norm_non)

make_pcoa_plot <- function(physeq) {
  ord = ordinate(physeq, method="PCoA", distance="jaccard")
  ordplot=plot_ordination(physeq, ord, color = 'FinePMI')
  return(ordplot)
}


beta_diversity_calc<- function(physeq) {
  GPdist=phyloseq::distance(physeq, "(A+B-2*J)/(A+B-J)")
  sampledf <- data.frame(sample_data(physeq))
  print(return(adonis(GPdist ~ BroadPMI, data = sampledf)))
}

beta_dispersion_calc <- function(physeq) {
  GPdist=phyloseq::distance(physeq, "(A+B-2*J)/(A+B-J)")
  sampledf <- data.frame(sample_data(physeq))
  beta <- betadisper(GPdist, sampledf$BroadPMI)
  print(return(permutest(beta)))
}