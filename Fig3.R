rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/ThesisProject_AKP/")

#packages
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(microbiome)
library(phyloseq)
library(dplyr)
library(PMCMR)
library(tidyverse)
library(vegan)

#upload the data
tax_group <- read.csv("tax_group_AKP_edit.csv")
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

physeq_ear_norm_Core <- prune_taxa(core_list_ear_tax, physeq_ear_norm)
physeq_eye_norm_Core <- prune_taxa(core_list_eye_tax, physeq_eye_norm)
physeq_nos_norm_Core <- prune_taxa(core_list_nos_tax, physeq_nos_norm)
physeq_mou_norm_Core <- prune_taxa(core_list_mou_tax, physeq_mou_norm)
physeq_rec_norm_Core <- prune_taxa(core_list_rec_tax, physeq_rec_norm)

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

physeq_ear_norm_non <- noncore_physeq(physeq_ear_norm, core_list_ear_tax)
physeq_eye_norm_non <- noncore_physeq(physeq_eye_norm, core_list_eye_tax)
physeq_rec_norm_non <- noncore_physeq(physeq_rec_norm, core_list_rec_tax)
physeq_mou_norm_non <- noncore_physeq(physeq_mou_norm, core_list_mou_tax)
physeq_nos_norm_non <- noncore_physeq(physeq_nos_norm, core_list_nos_tax)


ord = ordinate(physeq_ear_norm_Core, method="PCoA", distance="jaccard")
ordplot_a=plot_ordination(physeq_ear_norm_Core, ord, color = 'FinePMI')
ordplot_a

ord_non <- ordinate(physeq_ear_norm_non, method="PCoA", distance="jaccard")
ordplot_b =plot_ordination(physeq_ear_norm_non, ord_non, color = 'FinePMI')
ordplot_b

theme_set(theme_classic(base_size = 18))
tiff("Supp_ear.TIF", width = 3000, height = 1500, res=300)
ggarrange(ordplot_a,ordplot_b, 
          labels = c("A", "B"), ncol = 2)
dev.off()

ord = ordinate(physeq_nos_norm_Core, method="PCoA", distance="jaccard")
ordplot_a=plot_ordination(physeq_nos_norm_Core, ord, color = 'FinePMI')
ordplot_a

ord_non <- ordinate(physeq_nos_norm_non, method="PCoA", distance="jaccard")
ordplot_b =plot_ordination(physeq_nos_norm_non, ord_non, color = 'FinePMI')
ordplot_b

theme_set(theme_classic(base_size = 18))
tiff("Supp_nos.TIF", width = 3000, height = 1500, res=300)
ggarrange(ordplot_a,ordplot_b, 
          labels = c("A", "B"), ncol = 2)
dev.off()

samps_mou <- as.character(sample_names(physeq_mou_norm_Core))
samps_mou_not <- c('WCME_513')
samps_mou_want <- samps_mou[!samps_mou %in% samps_mou_not]
physeq_mou_norm_Core <- prune_samples(samps_mou_want, physeq_mou_norm_Core)

ord = ordinate(physeq_mou_norm_Core, method="PCoA", distance="jaccard")
ordplot_a=plot_ordination(physeq_mou_norm_Core, ord, color = 'FinePMI')
ordplot_a

ord_non <- ordinate(physeq_mou_norm_non, method="PCoA", distance="jaccard")
ordplot_b =plot_ordination(physeq_mou_norm_non, ord_non, color = 'FinePMI')
ordplot_b

theme_set(theme_classic(base_size = 18))
tiff("Supp_mou.TIF", width = 3000, height = 1500, res=300)
ggarrange(ordplot_a,ordplot_b, 
          labels = c("A", "B"), ncol = 2)
dev.off()

samps_rec <- as.character(sample_names(physeq_rec_norm_Core))
samps_rec_not <- c('WCME_1084', 'WCME_1086', 'WCME_1106', 'WCME_1164',
                   'WCME_1309',  'WCME_578',  'WCME_601',  'WCME_813',
                   'WCME_866',  'WCME_925')
samps_rec_want <- samps_rec[!samps_rec %in% samps_rec_not]
physeq_rec_norm_Core <- prune_samples(samps_rec_want, physeq_rec_norm_Core)

ord = ordinate(physeq_rec_norm_Core, method="PCoA", distance="jaccard")
ordplot_a=plot_ordination(physeq_rec_norm_Core, ord, color = 'FinePMI')
ordplot_a

ord_non <- ordinate(physeq_rec_norm_non, method="PCoA", distance="jaccard")
ordplot_b =plot_ordination(physeq_rec_norm_non, ord_non, color = 'FinePMI')
ordplot_b

theme_set(theme_classic(base_size = 18))
tiff("Supp_rec.TIF", width = 3000, height = 1500, res=300)
ggarrange(ordplot_a,ordplot_b, 
          labels = c("A", "B"), ncol = 2)
dev.off()

samps_eye <- as.character(sample_names(physeq_eye_norm_Core))
samps_eye_not <- c('WCME_1235', 'WCME_167', 'WCME_512', 'WCME_970')
samps_eye_want <- samps_eye[!samps_eye %in% samps_eye_not]
physeq_eye_norm_Core <- prune_samples(samps_eye_want, physeq_eye_norm_Core)
      
ord = ordinate(physeq_eye_norm_Core, method="PCoA", distance="jaccard")
ordplot_a=plot_ordination(physeq_eye_norm_Core, ord, color = 'FinePMI')
ordplot_a

ord_non <- ordinate(physeq_eye_norm_non, method="PCoA", distance="jaccard")
ordplot_b =plot_ordination(physeq_eye_norm_non, ord_non, color = 'FinePMI')
ordplot_b

theme_set(theme_classic(base_size = 18))
tiff("Supp_eye.TIF", width = 3000, height = 1500, res=300)
ggarrange(ordplot_a,ordplot_b, 
          labels = c("A", "B"), ncol = 2)
dev.off()

physeq_core_taxa <- merge_phyloseq(physeq_ear_norm_Core, physeq_eye_norm_Core, physeq_rec_norm_Core, 
                                   physeq_mou_norm_Core, physeq_nos_norm_Core)

ord = ordinate(physeq_core_taxa, method="PCoA", distance="jaccard")
a=plot_ordination(physeq_core_taxa, ord, color = 'FinePMI') 
a

physeq_non_taxa <- merge_phyloseq(physeq_ear_norm_non, physeq_eye_norm_non, physeq_rec_norm_non, 
                                   physeq_mou_norm_non, physeq_nos_norm_non)

ord_non <- ordinate(physeq_non_taxa, method="PCoA", distance="jaccard")
b =plot_ordination(physeq_non_taxa, ord_non, color = 'FinePMI')
b


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

ord = ordinate(piphy_ear_norm_Core, method="PCoA", distance="jaccard")
ordplot_a=plot_ordination(piphy_ear_norm_Core, ord, color = 'FinePMI')
ordplot_a

ord_non <- ordinate(piphy_ear_norm_non, method="PCoA", distance="jaccard")
ordplot_b =plot_ordination(piphy_ear_norm_non, ord_non, color = 'FinePMI')
ordplot_b

theme_set(theme_classic(base_size = 18))
tiff("Supp_ear_func.TIF", width = 3000, height = 1500, res=300)
ggarrange(ordplot_a,ordplot_b, 
          labels = c("A", "B"), ncol = 2)
dev.off()

ord = ordinate(piphy_nos_norm_Core, method="PCoA", distance="jaccard")
ordplot_a=plot_ordination(piphy_nos_norm_Core, ord, color = 'FinePMI')
ordplot_a

ord_non <- ordinate(piphy_nos_norm_non, method="PCoA", distance="jaccard")
ordplot_b =plot_ordination(piphy_nos_norm_non, ord_non, color = 'FinePMI')
ordplot_b

theme_set(theme_classic(base_size = 18))
tiff("Supp_nos_func.TIF", width = 3000, height = 1500, res=300)
ggarrange(ordplot_a,ordplot_b, 
          labels = c("A", "B"), ncol = 2)
dev.off()

ord = ordinate(piphy_mou_norm_Core, method="PCoA", distance="jaccard")
ordplot_a=plot_ordination(piphy_mou_norm_Core, ord, color = 'FinePMI')
ordplot_a

ord_non <- ordinate(piphy_mou_norm_non, method="PCoA", distance="jaccard")
ordplot_b =plot_ordination(piphy_mou_norm_non, ord_non, color = 'FinePMI')
ordplot_b

theme_set(theme_classic(base_size = 18))
tiff("Supp_mou_func.TIF", width = 3000, height = 1500, res=300)
ggarrange(ordplot_a,ordplot_b, 
          labels = c("A", "B"), ncol = 2)
dev.off()

ord = ordinate(piphy_rec_norm_Core, method="PCoA", distance="jaccard")
ordplot_a=plot_ordination(piphy_rec_norm_Core, ord, color = 'FinePMI')
ordplot_a

ord_non <- ordinate(piphy_rec_norm_non, method="PCoA", distance="jaccard")
ordplot_b =plot_ordination(piphy_rec_norm_non, ord_non, color = 'FinePMI')
ordplot_b

theme_set(theme_classic(base_size = 18))
tiff("Supp_rec_func.TIF", width = 3000, height = 1500, res=300)
ggarrange(ordplot_a,ordplot_b, 
          labels = c("A", "B"), ncol = 2)
dev.off()

ord = ordinate(piphy_eye_norm_Core, method="PCoA", distance="jaccard")
ordplot_a=plot_ordination(piphy_eye_norm_Core, ord, color = 'FinePMI')
ordplot_a

ord_non <- ordinate(piphy_eye_norm_non, method="PCoA", distance="jaccard")
ordplot_b =plot_ordination(piphy_eye_norm_non, ord_non, color = 'FinePMI')
ordplot_b

theme_set(theme_classic(base_size = 18))
tiff("Supp_eye_func.TIF", width = 3000, height = 1500, res=300)
ggarrange(ordplot_a,ordplot_b, 
          labels = c("A", "B"), ncol = 2)
dev.off()

piphy_core_taxa <- merge_phyloseq(piphy_ear_norm_Core, piphy_eye_norm_Core, piphy_rec_norm_Core, 
                                   piphy_mou_norm_Core, piphy_nos_norm_Core)

ord = ordinate(piphy_core_taxa, method="PCoA", distance="jaccard")
c=plot_ordination(piphy_core_taxa, ord, color = 'FinePMI') 
c

piphy_non_taxa <- merge_phyloseq(piphy_ear_norm_non, piphy_eye_norm_non, piphy_rec_norm_non, 
                                  piphy_mou_norm_non, piphy_nos_norm_non)

ord_non <- ordinate(piphy_non_taxa, method="PCoA", distance="jaccard")
d =plot_ordination(piphy_non_taxa, ord_non, color = 'FinePMI')
d


#Final Figure
theme_set(theme_classic(base_size = 18))
tiff("Fig3.TIF", width = 4000, height = 4000, res=300)
ggarrange(a,b,c,d, 
          labels = c("A", "B", "C", "D"),
          nrow = 2, ncol = 2)
dev.off()
