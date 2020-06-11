rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/ThesisProject_AKP/")

#packages
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(lme4)
library(phyloseq)
library(piecewiseSEM)
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

core_percent_rec <- proportion_of_core(physeq_rec_norm, core_list_rec)
core_percent_eye <- proportion_of_core(physeq_eye_norm, core_list_eye)
core_percent_ear <- proportion_of_core(physeq_ear_norm, core_list_ear)
core_percent_mou <- proportion_of_core(physeq_mou_norm, core_list_mou)
core_percent_nos <- proportion_of_core(physeq_nos_norm, core_list_nos)


core_df <- rbind(core_percent_rec, core_percent_eye, core_percent_ear, core_percent_mou, core_percent_nos)

core_df <- merge(core_df, metadata, by = 'SampleID')
write.csv(core_df, 'core_percentage_norare.csv')

core_data <- read.csv('core_percentage_norare.csv', header= T)

core_data$FinePMI <- factor(core_data$FinePMI, levels = c("Less24", "25-48", "49-72", "Great73"))

#models
m0 <- lmer(Core_percent ~ 1 + (1|Pack_ID), data = core_data, REML = F)
m1 <- lmer(Core_percent ~ FinePMI + (1|Pack_ID), data = core_data, REML = F)
m1.1 <- lmer(Core_percent ~ Sample_Area + (1|Pack_ID), data = core_data, REML = F)
m2 <- lmer(Core_percent ~ FinePMI + Sample_Area + (1|Pack_ID), data = core_data, REML = F)
m3 <- lmer(Core_percent ~ FinePMI + Sample_Area + CoD_Simple + (1|Pack_ID), data = core_data, REML = F)
m4 <- lmer(Core_percent ~ FinePMI + Sample_Area + Violent + (1|Pack_ID), data = core_data, REML = F)
m5 <- lmer(Core_percent ~ FinePMI + Sample_Area + Violent + CoD_Simple +(1|Pack_ID), data = core_data, REML = F)
m6 <- lmer(Core_percent ~ FinePMI + Sample_Area + Violent * CoD_Simple + (1|Pack_ID), data = core_data, REML = F)

anova(m0, m1, m1.1, m2, m3, m4, m5, m6, test = 'Chisq')

summary(m2)
rsquared(m2)

#functions
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

core_data <- read.csv('core_percentage_func.csv', header= T)

core_data$FinePMI <- factor(core_data$FinePMI, levels = c("Less24", "25-48", "49-72", "Great73"))

ggplot(data = core_data, aes(x=FinePMI, y = Core_percent, color = Sample_Area)) +
  geom_boxplot()

#add random effects of family/genus 

m0 <- lmer(Core_percent ~ 1 + (1|Pack_ID), data = core_data, REML = F)
m1 <- lmer(Core_percent ~ FinePMI + (1|Pack_ID), data = core_data, REML = F)
m1.1 <- lmer(Core_percent ~ Sample_Area + (1|Pack_ID), data = core_data, REML = F)
m2 <- lmer(Core_percent ~ FinePMI + Sample_Area + (1|Pack_ID), data = core_data, REML = F)
m3 <- lmer(Core_percent ~ FinePMI + Sample_Area + CoD_Simple + (1|Pack_ID), data = core_data, REML = F)
m4 <- lmer(Core_percent ~ FinePMI + Sample_Area + Violent + (1|Pack_ID), data = core_data, REML = F)
m5 <- lmer(Core_percent ~ FinePMI + Sample_Area + Violent + CoD_Simple +(1|Pack_ID), data = core_data, REML = F)
m6 <- lmer(Core_percent ~ FinePMI + Sample_Area + Violent * CoD_Simple + (1|Pack_ID), data = core_data, REML = F)

anova(m0, m1, m1.1, m2, m3, m4, m5, m6, test = 'Chisq')
anova(m0, m2, test = 'Chisq')

summary(m2)
rsquared(m2)

