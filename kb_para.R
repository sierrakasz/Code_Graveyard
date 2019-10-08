#packages
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(HMP)
library(phyloseq)
library(plyr)
library(PMCMR)
library(randomForest)
library(tidyverse)


#make a table in excel for rare data
tax_group_nr <- read.csv("tax_group_norare.csv")
tax_group_nr <- tax_group_nr %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group_nr <- tax_group_nr %>%  ungroup()

tax_group_7 <- read.csv("tax_group_7000.csv")
tax_group_7 <- tax_group_7 %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group_7 <- tax_group_7 %>%  ungroup()

tax_group_1 <- read.csv("tax_group_1000.csv")
tax_group_1 <- tax_group_1 %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group_1 <- tax_group_1 %>%  ungroup()

tax_group_rare <- read.csv("tax_group_rare.csv")
tax_group_rare <- tax_group_rare %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group_rare <- tax_group_rare %>%  ungroup()

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

physeq_nr <- regroup_physeq_object(tax_group_nr)
physeq_7 <- regroup_physeq_object(tax_group_7)
physeq_1 <- regroup_physeq_object(tax_group_1)
physeq_rare <- regroup_physeq_object(tax_group_rare)

#import metadata and combine
metadata=(read.csv("Metadata/HPMMMeta_r_merge_rare.csv",header=TRUE))
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID

metadata_rare=(read.csv("Metadata/HPMMMeta_r_merge_rare_levels.csv",header=TRUE))
sampdat_rare=sample_data(metadata_rare)
sampdat_rare$Rarefy <- factor(sampdat_rare$Rarefy, levels = c("No Rarefaction", "7000", "1000"))
sample_names(sampdat_rare)=metadata_rare$SampleID

merger = merge_phyloseq(physeq_rare, sampdat_rare)
sample_data(merger)$Rarefy <- factor(sample_data(merger)$Rarefy, levels = c("No Rarefaction", "7000", "1000"))

physeq_nr=merge_phyloseq(physeq_nr, sampdat)
physeq_7=merge_phyloseq(physeq_7, sampdat)
physeq_1=merge_phyloseq(physeq_1, sampdat)

physeq_nr <- subset_samples(merger, Rarefy == 'No Rarefaction')
physeq_7 <- subset_samples(merger, Rarefy == '7000')
physeq_1 <- subset_samples(merger, Rarefy == '1000')

physeq_60 <- subset_samples(merger, Subsample == '60')
physeq_120 <- subset_samples(merger, Subsample == '120')
physeq_188 <- subset_samples(merger, Subsample == '188')

data_nr <- as.data.frame(t(otu_table(physeq_nr)))
data_7 <- as.data.frame(t(otu_table(physeq_7)))
data_1 <- as.data.frame(t(otu_table(physeq_1)))

group.data <- list(data_nr, data_7, data_1)

## Not run:
kl <- Kullback.Leibler(group.data, plot = FALSE)
kl

data_60 <- as.data.frame(t(otu_table(physeq_60)))
data_120 <- as.data.frame(t(otu_table(physeq_120)))
data_188 <- as.data.frame(t(otu_table(physeq_188)))

group.data <- list(data_60, data_120, data_188)

## Not run:
kl <- Kullback.Leibler(group.data, plot = FALSE)
kl



# random forest -----------------------------------------------------------

random_foresting_PMI <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$FinePMI)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data, ntree=1000)
  return(Forest)
}
