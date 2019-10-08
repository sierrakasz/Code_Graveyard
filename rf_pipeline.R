
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(HMP)
library(phyloseq)
library(plyr)
library(PMCMR)
library(randomForest)
library(rsample)
library(tidyverse)
library(phyloseq)
library(tidyverse)
library(vegan)

###tax group is a file that combines the OTU file and taxonomy file. 
##combining the files made it easier for collapsing differences among the pipeline outputs
#tax.group is available on github

tax_group <- read.csv("tax_group.csv")
tax_group <- tax_group %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group <- tax_group %>%  ungroup()

#this function specifically collapses the OTU table and taxonomy file 
#when the taxonomy is the same name among pipelines
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

#import metadata and combine into one phyloseq object
metadata=(read.csv("Metadata/HPMMMeta_r_merge.csv",header=TRUE))
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
physeq=merge_phyloseq(physeq_all, sampdat)
physeq


physeq_q <- subset_samples(physeq, Pipeline == 'QIIME2')
physeq_m <- subset_samples(physeq, Pipeline == 'mothur')
physeq_mg <- subset_samples(physeq, Pipeline == 'MG-RAST')

physeq_q_p <- tax_glom(physeq_q, taxrank = 'Phylum')
physeq_m_p <- tax_glom(physeq_m, taxrank = 'Phylum')
physeq_mg_p <- tax_glom(physeq_mg, taxrank = 'Phylum')
physeq_q_f <- tax_glom(physeq_q, taxrank = 'Family')
physeq_m_f <- tax_glom(physeq_m, taxrank = 'Family')
physeq_mg_f <- tax_glom(physeq_mg, taxrank = 'Family')


random_foresting_sample_area <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$Sample_Area)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data)
  return(Forest)
}

forest_q_p <- random_foresting_sample_area(physeq_q_p)
forest_m_p <- random_foresting_sample_area(physeq_m_p)
forest_mg_p <- random_foresting_sample_area(physeq_mg_p)
forest_q_f <- random_foresting_sample_area(physeq_q_f)
forest_m_f <- random_foresting_sample_area(physeq_m_f)
forest_mg_f <- random_foresting_sample_area(physeq_mg_f)

forest_predictors_sample_area <- function(forest) {
  imp <- importance(forest)
  imp <- data.frame(predictors = rownames(imp), imp)
  imp.sort <- arrange(imp, desc(MeanDecreaseGini))
  imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
  imp.20 <- imp.sort[1:20, ]
  tax <- data.frame(tax_table(physeq))
  head(tax)
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family)
  tax$predictors <- rownames(tax)
  imp.20 <- merge(imp.20, tax)
  write.table(imp.20, sep=",", "random_forest_predictors_pipe_p.csv", append = TRUE)
  return(imp.20)
}


pred_q_f <- forest_predictors_sample_area(forest_q_f)
pred_m_f <- forest_predictors_sample_area(forest_m_f)
pred_mg_f <- forest_predictors_sample_area(forest_mg_f)

random_foresting_MoD <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$MoD)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data, ntree=1000)
  return(Forest)
}

forest_q_p <- random_foresting_MoD(physeq_q_p)
forest_m_p <- random_foresting_MoD(physeq_m_p)
forest_mg_p <- random_foresting_MoD(physeq_mg_p)
forest_q_f <- random_foresting_MoD(physeq_q_f)
forest_m_f <- random_foresting_MoD(physeq_m_f)
forest_mg_f <- random_foresting_MoD(physeq_mg_f)

forest_predictors_MoD <- function(forest) {
  imp <- importance(forest)
  imp <- data.frame(predictors = rownames(imp), imp)
  imp.sort <- arrange(imp, desc(MeanDecreaseGini))
  imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
  imp.20 <- imp.sort[1:20, ]
  tax <- data.frame(tax_table(physeq))
  head(tax)
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family)
  tax$predictors <- rownames(tax)
  imp.20 <- merge(imp.20, tax)
  write.table(imp.20, sep=",", "random_forest_predictors_pipe_md.csv", append = TRUE)
  return(imp.20)
}


pred_q_f <- forest_predictors_MoD(forest_q_f)
pred_m_f <- forest_predictors_MoD(forest_m_f)
pred_mg_f <- forest_predictors_MoD(forest_mg_f)


# KL div ------------------------------------------------------------------
data_mg <- as.data.frame(t(otu_table(physeq_mg)))
data_m <- as.data.frame(t(otu_table(physeq_m)))
data_q <- as.data.frame(t(otu_table(physeq_q)))

group.data <- list(data_mg, data_m, data_q)

## Not run:
kl <- Kullback.Leibler(group.data, plot = FALSE)
kl

