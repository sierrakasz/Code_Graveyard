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

piphy <- import_biom("table-with-taxonomy.biom")

metadata=(read.csv("HPMMMetadata_pt2_cod.csv",header=TRUE))
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID

piphy=merge_phyloseq(piphy, sampdat)

piphy_rec <- subset_samples(piphy, Sample_Area == 'Rectum')
piphy_eye <- subset_samples(piphy, Sample_Area == 'Eyes')
piphy_ear <- subset_samples(piphy, Sample_Area == 'Ears')
piphy_mou <- subset_samples(piphy, Sample_Area == 'Nose')
piphy_nos <- subset_samples(piphy, Sample_Area == 'Mouth')

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

piphy_rec_norm <- normalize_wout_rarefying(piphy_rec)
piphy_eye_norm <- normalize_wout_rarefying(piphy_eye)
piphy_ear_norm <- normalize_wout_rarefying(piphy_ear)
piphy_mou_norm <- normalize_wout_rarefying(piphy_mou)
piphy_nos_norm <- normalize_wout_rarefying(piphy_nos)

GPfr_rec_24 <-subset_samples(piphy_rec_norm, FinePMI == 'Less24')
GPfr_eye_24 <-subset_samples(piphy_eye_norm, FinePMI == 'Less24')
GPfr_ear_24 <-subset_samples(piphy_ear_norm, FinePMI == 'Less24')
GPfr_nos_24 <-subset_samples(piphy_nos_norm, FinePMI == 'Less24')
GPfr_mou_24 <-subset_samples(piphy_mou_norm, FinePMI == 'Less24')

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

core_list_rec <- lets_find_the_core(GPfr_rec_24)
core_list_eye <- lets_find_the_core(GPfr_eye_24)
core_list_ear <- lets_find_the_core(GPfr_ear_24)
core_list_nos <- lets_find_the_core(GPfr_nos_24)
core_list_mou <- lets_find_the_core(GPfr_mou_24)

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

core_percent_rec <- proportion_of_core(piphy_rec_norm, core_list_rec)
core_percent_eye <- proportion_of_core(piphy_eye_norm, core_list_eye)
core_percent_ear <- proportion_of_core(piphy_ear_norm, core_list_ear)
core_percent_mou <- proportion_of_core(piphy_mou_norm, core_list_mou)
core_percent_nos <- proportion_of_core(piphy_nos_norm, core_list_nos)


core_df <- rbind(core_percent_rec, core_percent_eye, core_percent_ear, core_percent_mou, core_percent_nos)

core_df <- merge(core_df, metadata, by = 'SampleID')
write.csv(core_df, 'core_percentage_func.csv')

piphy_ear_norm_Core <- prune_taxa(core_list_ear, piphy_ear_norm)
piphy_eye_norm_Core <- prune_taxa(core_list_eye, piphy_eye_norm)
piphy_nos_norm_Core <- prune_taxa(core_list_nos, piphy_nos_norm)
piphy_mou_norm_Core <- prune_taxa(core_list_mou, piphy_mou_norm)
piphy_rec_norm_Core <- prune_taxa(core_list_rec, piphy_rec_norm)

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

piphy_ear_norm_non <- noncore_physeq(piphy_ear_norm, core_list_ear)
piphy_eye_norm_non <- noncore_physeq(piphy_eye_norm, core_list_eye)
piphy_rec_norm_non <- noncore_physeq(piphy_rec_norm, core_list_rec)
piphy_mou_norm_non <- noncore_physeq(piphy_mou_norm, core_list_mou)
piphy_nos_norm_non <- noncore_physeq(piphy_nos_norm, core_list_nos)

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

