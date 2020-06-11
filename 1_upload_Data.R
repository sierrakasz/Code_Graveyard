rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/ThesisProject_AKP/")

#packages
library(phyloseq)
library(plyr)
library(tidyverse)

otu <- read.csv("table.csv")
tax <- read.csv("taxonomy.csv")
tree <- read_tree('tree.nwk')

rownames(otu) <- otu$OTUID
otu <- otu[,-1]

OTU=otu_table(otu, taxa_are_rows=TRUE)
physeq_all=phyloseq(OTU,tree)
physeq <- rarefy_even_depth(physeq_all, sample.size = 5000)


tax <- tax %>% 
  mutate(Kingdom = str_replace(Kingdom, "D_0__", "")) %>% 
  mutate(Phylum = str_replace(Phylum, "D_1__", "")) %>% 
  mutate(Class = str_replace(Class, "D_2__", "")) %>% 
  mutate(Order = str_replace(Order, "D_3__", "")) %>% 
  mutate(Family = str_replace(Family, "D_4__", "")) %>% 
  mutate(Genus = str_replace(Genus, "D_5__", "")) %>% 
  mutate(Species = str_replace(Species, "D_6__", ""))
tax[tax==""]<-"Unassigned"
tax_group <- merge(otu, tax, by = 'OTUID')

write.csv(tax_group, "tax_group_AKP.csv")

#upload the data
tax_group <- read.csv("tax_group_AKP.csv")
tax_group <- tax_group %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group <- tax_group %>%  ungroup()

write.csv(tax_group, 'tax_group_collapse.csv')

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
physeq <- rarefy_even_depth(physeq, sample.size = 5000)
otu <- data.frame(otu_table(physeq))
tax <- data.frame(tax_table(physeq))
otu$OTUID <- row.names(otu)
tax$OTUID <- row.names(tax)
merger <- merge(otu, tax, by = 'OTUID')
write.csv(merger, "tax_group_AKP_edit.csv")
#for editting names
