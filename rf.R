rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/Forensic_Pig/Combined")

#packages
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(phyloseq)
library(plyr)
library(PMCMR)
library(randomForest)
library(rsample)
library(tidyverse)

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
sampdat$Day <- factor(sampdat$Day, levels = c('1', '4', '7', '8', '10', '12', 
                                              '13', '16', '19', '21'))
sampdat$Decomp_Stage <- factor(sampdat$Decomp_Stage, levels = c('Submerged_Fresh', 'Early_Floating',
                                                                'Floating_Decay', 'Advanced_Floating_Decay', 
                                                                'Sunken_Remains'))
physeq=merge_phyloseq(physeq, sampdat)

physeq_benbow <- subset_samples(physeq, Study == "Benbow")
physeq_wallace <- subset_samples(physeq, Study == "Wallace")

#functions
random_foresting_Decomp <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$Decomp_Stage)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data, ntree=1000)
  return(Forest)
}

random_foresting_Decomp(physeq_benbow)
random_foresting_Decomp(physeq_wallace)

random_foresting_Day <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$Day)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data, ntree=1000)
  return(Forest)
}

random_foresting_Day(physeq_benbow)
random_foresting_Day(physeq_wallace)

forest_predictors <- function(forest) {
  imp <- importance(forest)
  imp <- data.frame(predictors = rownames(imp), imp)
  imp.sort <- arrange(imp, desc(MeanDecreaseGini))
  imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
  imp.20 <- imp.sort[1:20, ]
  tax <- data.frame(tax_table(physeq))
  head(tax)
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
  tax$predictors <- rownames(tax)
  imp.20 <- merge(imp.20, tax)
  write.table(imp.20, sep=",", "random_forest_predictors.csv", append = TRUE)
  return(imp.20)
}

for1 <- random_foresting_Decomp(physeq_benbow)
for2 <- random_foresting_Decomp(physeq_wallace)
for3 <- random_foresting_Day(physeq_benbow)
for4 <- random_foresting_Day(physeq_wallace)

forest_predictors(for1)
forest_predictors(for2)
forest_predictors(for3)
forest_predictors(for4)
