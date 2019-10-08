#packages
#library(ape)
#library(DESeq2)
#library(dplyr)
#library(dunn.test)
#library(ggforce)
#library(ggplot2)
#library(ggpubr)
#library(gplots)
#library(ggthemes)
#library(knitr)
#library(lme4)
#library(microbiome)
library(phyloseq)
#library(plyr)
#library(plotly)
#library(randomForest)
library(tidyr)
library(tidyverse)
library(tidyselect)
#library(UpSetR)
library(vegan)

rm(list=ls())

##import qiime files
#subthirty
biom_qiime_30 <- import_biom("qiime biom files/qiime.thirty.biom")

#remove D_ in name
#fixing taxa names from QIIME
a <- data.frame(tax_table(biom_qiime_30))
b <- a %>% 
  mutate(Rank1 = str_replace(Rank1, "D_0__", "")) %>% 
  mutate(Rank2 = str_replace(Rank2, "D_1__", "")) %>% 
  mutate(Rank3 = str_replace(Rank3, "D_2__", "")) %>% 
  mutate(Rank4 = str_replace(Rank4, "D_3__", "")) %>% 
  mutate(Rank5 = str_replace(Rank5, "D_4__", "")) %>% 
  mutate(Rank6 = str_replace(Rank6, "D_5__", ""))
c <- b[,-7]
c <- as.matrix(c)

d <- data.frame(otu_table(biom_qiime_30))

OTUq=otu_table(d, taxa_are_rows=TRUE)

TAXq=tax_table(c)

taxa_names(TAXq)=row.names(OTUq)

physeq_qiime_30 = merge_phyloseq(OTUq, TAXq)

sample_names(physeq_qiime_30) <- paste("Q_Old", sample_names(physeq_qiime_30), sep="")

colnames(tax_table(physeq_qiime_30)) <- c("Kingdom", "Phylum", "Class",
                                          "Order","Family", "Genus")

biom_qiime_30 <- import_biom("qiime biom files/table_tax.biom")

#remove D_ in name
#fixing taxa names from QIIME
a <- data.frame(tax_table(biom_qiime_30))
b <- a %>% 
  mutate(Rank1 = str_replace(Rank1, "D_0__", "")) %>% 
  mutate(Rank2 = str_replace(Rank2, "D_1__", "")) %>% 
  mutate(Rank3 = str_replace(Rank3, "D_2__", "")) %>% 
  mutate(Rank4 = str_replace(Rank4, "D_3__", "")) %>% 
  mutate(Rank5 = str_replace(Rank5, "D_4__", "")) %>% 
  mutate(Rank6 = str_replace(Rank6, "D_5__", ""))
c <- b[,-7]
c <- as.matrix(c)

d <- data.frame(otu_table(biom_qiime_30))

OTUq=otu_table(d, taxa_are_rows=TRUE)

TAXq=tax_table(c)

taxa_names(TAXq)=row.names(OTUq)

physeq_qiime_30n = merge_phyloseq(OTUq, TAXq)

sample_names(physeq_qiime_30n) <- paste("Q_New", sample_names(physeq_qiime_30n), sep="")

colnames(tax_table(physeq_qiime_30n)) <- c("Kingdom", "Phylum", "Class",
                                          "Order","Family", "Genus")

rm(a,b,biom_qiime_30,c,d,OTUq,TAXq)

#mothur
physeq_mothur_30 <- import_mothur(mothur_shared_file = "mothur biom files/mothur.thirty.otu.shared", 
                                  mothur_constaxonomy_file = "mothur biom files/mothur.thirty.tax.taxonomy")

sample_names(physeq_mothur_30) <- paste("M_Old", sample_names(physeq_mothur_30), sep="")

colnames(tax_table(physeq_mothur_30)) <- c("Kingdom", "Phylum", "Class", 
                                           "Order", "Family",  "Genus")

physeq_mothur_30n <- import_mothur(mothur_shared_file = "mothur biom files/mothur.thirty.new.shared", 
                                  mothur_constaxonomy_file = "mothur biom files/mothur.thirty.new.taxonomy")

sample_names(physeq_mothur_30n) <- paste("M_New", sample_names(physeq_mothur_30n), sep="")

colnames(tax_table(physeq_mothur_30n)) <- c("Kingdom", "Phylum", "Class", 
                                           "Order", "Family",  "Genus")

#import MG-RAST
#subthirty
otufull=read.csv("MG-Rast biom files\\otu_mgrast_thirty.csv",header=TRUE)
taxa=as.matrix(read.csv("MG-Rast biom files\\taxa_thirty.csv"))

OTU=otu_table(otufull, taxa_are_rows=TRUE)

TAX=tax_table(taxa)

taxa_names(TAX)=row.names(OTU)

sample_names(OTU) <- paste("G_Old", sample_names(OTU), sep="")

physeq_MG_30=phyloseq(OTU,TAX)
colnames(tax_table(physeq_MG_30)) <- c("Kingdom", "Phylum", "Class", 
                                       "Order", "Family",  "Genus")
rm(otufull, taxa, OTU, TAX)

combined=read.csv("MG-Rast biom files\\HPMMS30.csv",header=TRUE)
tax <- combined %>% select(domain, phylum, className, order, family, genus)
otu <- combined %>% select(contains("WCME"))
tax <- as.matrix(tax)

OTU=otu_table(otu, taxa_are_rows=TRUE)
TAX=tax_table(tax)

taxa_names(TAX)=row.names(OTU)

sample_names(OTU) <- paste("G_New", sample_names(OTU), sep="")

physeq_MG_30n=phyloseq(OTU,TAX)
colnames(tax_table(physeq_MG_30n)) <- c("Kingdom", "Phylum", "Class", 
                                       "Order", "Family",  "Genus")
rm(OTU,TAX, otu, tax, combined)

merge_qiime = merge_phyloseq(physeq_qiime_30, physeq_qiime_30n)
merge_mothur = merge_phyloseq(physeq_mothur_30, physeq_mothur_30n)
merge_MG = merge_phyloseq(physeq_MG_30, physeq_MG_30n)
merge = merge_phyloseq(merge_mothur, merge_MG, merge_qiime)

#tidying
x <- data.frame(otu_table(merge))
y <- data.frame(tax_table(merge))
x$names <- rownames(x)
y$names <- rownames(y)

z <- merge(x,y)
head(z)
z <- z[,-1]

tax_group <- z %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))

tax_group <- tax_group %>% filter(Kingdom == "Bacteria")

#fixed taxa by hand
write.csv(tax_group, "tax_group.csv")

#fix it by hand
tax_group <- read.csv("tax_group.csv")
tax_group_new <- tax_group %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
head(tax_group)

tax_group <- tax_group %>% ungroup()

tax <- tax_group %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
otu <- tax_group %>% select(contains("WCME"))
tax<- as.matrix(tax)

OTU=otu_table(otu, taxa_are_rows=TRUE)

TAX=tax_table(tax)

taxa_names(TAX)=row.names(OTU)

physeq_all=phyloseq(OTU,TAX)

rm(otu, tax, tax_group, OTU, TAX)

metadata=(read.csv("Metadata/HPMMMeta_r_oldnew.csv",header=TRUE))

sampdat=sample_data(metadata)
row.names(sampdat)
sample_names(sampdat)=metadata$SampleID

physeq=merge_phyloseq(physeq_all, sampdat)

physeq <- rarefy_even_depth(physeq, replace = FALSE, trimOTUs = TRUE)
physeq

physeq_new <- subset_samples(physeq, Age="New")
physeq_old <- subset_samples(physeq, Age="Old")

ord = ordinate(physeq_new, method="PCoA", distance="jaccard")
ordplot=plot_ordination(physeq_new, ord,"SampleID", color="Pipeline",
                        title ="PCoA of Jaccard Index\nOld Method vs. New Method")
# + stat_ellipse(geom = "polygon", alpha = 1/8, aes(color = Pipeline))
ordplot 

ordplot+ stat_ellipse(type= "norm",geom = "polygon", alpha = 1/4, aes(fill = Sample_Area))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
