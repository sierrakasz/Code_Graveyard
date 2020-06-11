rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/ThesisProject_AKP/")

#packages
library(arm)
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

GPfr_rec_24_tax <-subset_samples(physeq_rec_norm, FinePMI == 'Less24')
GPfr_eye_24_tax <-subset_samples(physeq_eye_norm, FinePMI == 'Less24')
GPfr_ear_24_tax <-subset_samples(physeq_ear_norm, FinePMI == 'Less24')
GPfr_nos_24_tax <-subset_samples(physeq_nos_norm, FinePMI == 'Less24')
GPfr_mou_24_tax <-subset_samples(physeq_mou_norm, FinePMI == 'Less24')

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

core_list_rec_tax <- lets_find_the_core(GPfr_rec_24_tax)
core_list_eye_tax <- lets_find_the_core(GPfr_eye_24_tax)
core_list_ear_tax <- lets_find_the_core(GPfr_ear_24_tax)
core_list_nos_tax <- lets_find_the_core(GPfr_nos_24_tax)
core_list_mou_tax <- lets_find_the_core(GPfr_mou_24_tax)

#Panel A
taxa_for_fig <- data.frame(tax_table(physeq)) 
taxa_for_fig$OTU <- rownames(taxa_for_fig)
core_data_frame_fig <- data.frame(c(core_list_rec_tax, core_list_eye_tax,
                                    core_list_ear_tax, core_list_nos_tax, core_list_mou_tax), 
                                  c(rep('Rectum', times = length(core_list_rec_tax)), rep('Eyes', times = length(core_list_eye_tax)),
                                    rep('Ears', times = length(core_list_ear_tax)), rep('Nose', times = length(core_list_nos_tax)),
                                    rep('Mouth', times = length(core_list_mou_tax))))
colnames(core_data_frame_fig) <- c('OTU', 'Body_Site')
core_merge_df <- merge(core_data_frame_fig, taxa_for_fig, by = 'OTU')

mean_of_taxa_core <- function(physeq, core_list) {
  otu_for_fig <- data.frame(otu_table(physeq))
  otu_for_fig$OTU <- rownames(otu_for_fig)
  otu_for_fig_cor <- filter(otu_for_fig, otu_for_fig$OTU %in% core_list)
  otu_for_fig_cor$Mean <- rowMeans(otu_for_fig_cor[,-length(otu_for_fig_cor)])
  new_df <- otu_for_fig_cor[,c('OTU', 'Mean')]
  return(new_df)
}

mean_core_ear <- mean_of_taxa_core(physeq_ear_norm, core_list_ear_tax)
mean_core_ear$Body_Site <- 'Ears'
mean_core_eye <- mean_of_taxa_core(physeq_eye_norm, core_list_eye_tax)
mean_core_eye$Body_Site <- 'Eyes'
mean_core_nos <- mean_of_taxa_core(physeq_nos_norm, core_list_nos_tax)
mean_core_nos$Body_Site <- 'Nose'
mean_core_mou <- mean_of_taxa_core(physeq_mou_norm, core_list_mou_tax)
mean_core_mou$Body_Site <- 'Mouth'
mean_core_rec <- mean_of_taxa_core(physeq_rec_norm, core_list_rec_tax)
mean_core_rec$Body_Site <- 'Rectum'

mean_core_merge <- rbind(mean_core_ear, mean_core_eye, mean_core_mou, mean_core_rec, mean_core_nos)
total_merge <- merge(core_merge_df, mean_core_merge, by = c('OTU', 'Body_Site'))


a <- ggplot(data = total_merge, aes(x=Genus, y = Mean, fill = Body_Site)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  facet_grid(Body_Site~., scales = "free", space = "free") + 
  coord_flip() + xlab("Core Genera") + ylab("Mean Number of Sequences Per Sample") +
  scale_fill_manual("Body Site", values = c('Ears' = '#ECA72F', 'Eyes' = '#3DC5B6', 
                                            'Nose' = '#D45924', 'Mouth' = '#F65E5A', 'Rectum' = '#90839F'))
a

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
m2 <- lmer(Core_percent ~ FinePMI + Sample_Area + (1|Pack_ID), data = core_data, REML = F)

#Panel B
core_data$fit <- predict(m2)


b <- ggplot(data = core_data, aes(x=FinePMI, y = Core_percent, fill = Sample_Area)) +
  geom_boxplot() + facet_wrap(~Sample_Area, nrow = 1) + 
  geom_point(aes(y=fit), shape = 2, show.legend = F) + ylab("Proportion of Core Taxa") + 
  xlab("Postmortem Interval") +
  scale_fill_manual("Body Site", values = c('Ears' = '#ECA72F', 'Eyes' = '#3DC5B6', 
                                            'Nose' = '#D45924', 'Mouth' = '#F65E5A', 'Rectum' = '#90839F')) +
  theme(axis.text.x = element_text(angle = 90))
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

#Panel C
func_for_fig <- data.frame(tax_table(piphy)) 
func_for_fig$KEGG <- rownames(func_for_fig)
core_data_frame_fig_func <- data.frame(c(core_list_rec, core_list_eye,
                                    core_list_ear, core_list_nos, core_list_mou), 
                                  c(rep('Rectum', times = length(core_list_rec)), rep('Eyes', times = length(core_list_eye)),
                                    rep('Ears', times = length(core_list_ear)), rep('Nose', times = length(core_list_nos)),
                                    rep('Mouth', times = length(core_list_mou))))
colnames(core_data_frame_fig_func) <- c('KEGG', 'Body_Site')
core_merge_df_func <- merge(core_data_frame_fig_func, func_for_fig, by = 'KEGG')

mean_of_func_core <- function(physeq, core_list) {
  otu_for_fig <- data.frame(otu_table(physeq))
  otu_for_fig$KEGG <- rownames(otu_for_fig)
  otu_for_fig_cor <- filter(otu_for_fig, otu_for_fig$KEGG %in% core_list)
  otu_for_fig_cor$Mean <- rowMeans(otu_for_fig_cor[,-length(otu_for_fig_cor)])
  new_df <- otu_for_fig_cor[,c('KEGG', 'Mean')]
  return(new_df)
}

mean_core_func_ear <- mean_of_func_core(piphy_ear_norm, core_list_ear)
mean_core_func_ear$Body_Site <- 'Ears'
mean_core_func_eye <- mean_of_func_core(piphy_eye_norm, core_list_eye)
mean_core_func_eye$Body_Site <- 'Eyes'
mean_core_func_nos <- mean_of_func_core(piphy_nos_norm, core_list_nos)
mean_core_func_nos$Body_Site <- 'Nose'
mean_core_func_mou <- mean_of_func_core(piphy_mou_norm, core_list_mou)
mean_core_func_mou$Body_Site <- 'Mouth'
mean_core_func_rec <- mean_of_func_core(piphy_rec_norm, core_list_rec)
mean_core_func_rec$Body_Site <- 'Rectum'

mean_core_func_merge <- rbind(mean_core_func_ear, mean_core_func_eye, mean_core_func_mou, mean_core_func_rec, mean_core_func_nos)
total_merge_func <- merge(core_merge_df_func, mean_core_func_merge, by = c('KEGG', 'Body_Site'))

hist(total_merge_func$Mean)
total_merge_func_out <- filter(total_merge_func, Mean < 30000)

c <- ggplot(data = total_merge_func_out, aes(x=Rank1, y = Mean, fill = Body_Site)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  facet_grid(Body_Site~., scales = "free", space = "free") + 
  coord_flip() + xlab("Core Functions") + ylab("Mean Number of KEGG Representatives Per Sample") +
  scale_fill_manual("Body Site", values = c('Ears' = '#ECA72F', 'Eyes' = '#3DC5B6', 
                                            'Nose' = '#D45924', 'Mouth' = '#F65E5A', 'Rectum' = '#90839F')) +
  theme(axis.text.y=element_blank())
c

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

#add random effects of family/genus 
m2 <- lmer(Core_percent ~ FinePMI + Sample_Area + (1|Pack_ID), data = core_data, REML = F)

#Panel D
core_data$fit <- predict(m2)


d <- ggplot(data = core_data, aes(x=FinePMI, y = Core_percent, fill = Sample_Area)) +
  geom_boxplot() + facet_wrap(~Sample_Area,  nrow = 1) + 
  geom_point(aes(y=fit), shape = 2, show.legend = F) + ylab("Proportion of Core Functions") + 
  xlab("Postmortem Interval") +
  scale_fill_manual("Body Site", values = c('Ears' = '#ECA72F', 'Eyes' = '#3DC5B6', 
                                            'Nose' = '#D45924', 'Mouth' = '#F65E5A', 'Rectum' = '#90839F')) +
  theme(axis.text.x = element_text(angle = 90))
d

#Final Figure
theme_set(theme_classic(base_size = 18))
tiff("Fig2.TIF", width = 4500, height = 4500, res=300)
ggarrange(a,b,c,d, 
          labels = c("A", "B", "C", "D"),
          nrow = 2, ncol = 2)
dev.off()

