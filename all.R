setwd("C:/Users/sierr/Documents/ThesisProject_Pipelines/")

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
#library(vegan)

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

sample_names(physeq_qiime_30) <- paste("Q_30_", sample_names(physeq_qiime_30), sep="")

colnames(tax_table(physeq_qiime_30)) <- c("Kingdom", "Phylum", "Class",
                                          "Order","Family", "Genus")
##import qiime files
#subsixty
biom_qiime_60 <- import_biom("qiime biom files/qiime.sixty.biom")

#remove D_ in name
#fixing taxa names from QIIME
a <- data.frame(tax_table(biom_qiime_60))
b <- a %>% 
  mutate(Rank1 = str_replace(Rank1, "D_0__", "")) %>% 
  mutate(Rank2 = str_replace(Rank2, "D_1__", "")) %>% 
  mutate(Rank3 = str_replace(Rank3, "D_2__", "")) %>% 
  mutate(Rank4 = str_replace(Rank4, "D_3__", "")) %>% 
  mutate(Rank5 = str_replace(Rank5, "D_4__", "")) %>% 
  mutate(Rank6 = str_replace(Rank6, "D_5__", ""))
c <- b[,-7]
c <- as.matrix(c)

d <- data.frame(otu_table(biom_qiime_60))

OTUq=otu_table(d, taxa_are_rows=TRUE)

TAXq=tax_table(c)

taxa_names(TAXq)=row.names(OTUq)

physeq_qiime_60 = merge_phyloseq(OTUq, TAXq)

sample_names(physeq_qiime_60) <- paste("Q_60_", sample_names(physeq_qiime_60), sep="")

colnames(tax_table(physeq_qiime_60)) <- c("Kingdom", "Phylum", "Class",
                                          "Order","Family", "Genus")

##import qiime files
#subonetwenty
biom_qiime_120 <- import_biom("qiime biom files/qiime.onetwenty.biom")

#remove D_ in name
#fixing taxa names from QIIME
a <- data.frame(tax_table(biom_qiime_120))
b <- a %>% 
  mutate(Rank1 = str_replace(Rank1, "D_0__", "")) %>% 
  mutate(Rank2 = str_replace(Rank2, "D_1__", "")) %>% 
  mutate(Rank3 = str_replace(Rank3, "D_2__", "")) %>% 
  mutate(Rank4 = str_replace(Rank4, "D_3__", "")) %>% 
  mutate(Rank5 = str_replace(Rank5, "D_4__", "")) %>% 
  mutate(Rank6 = str_replace(Rank6, "D_5__", ""))
c <- b[,-7]
c <- as.matrix(c)
d <- data.frame(otu_table(biom_qiime_120))

OTUq=otu_table(d, taxa_are_rows=TRUE)

TAXq=tax_table(c)

taxa_names(TAXq)=row.names(OTUq)

physeq_qiime_120 = merge_phyloseq(OTUq, TAXq)

sample_names(physeq_qiime_120) <- paste("Q_120_", sample_names(physeq_qiime_120), sep="")

colnames(tax_table(physeq_qiime_120)) <- c("Kingdom", "Phylum", "Class",
                                           "Order","Family", "Genus")

##import qiime files
#all
biom_qiime_all <- import_biom("qiime biom files/qiime.all.biom")

#remove D_ in name
#fixing taxa names from QIIME
a <- data.frame(tax_table(biom_qiime_all))
b <- a %>% 
  mutate(Rank1 = str_replace(Rank1, "D_0__", "")) %>% 
  mutate(Rank2 = str_replace(Rank2, "D_1__", "")) %>% 
  mutate(Rank3 = str_replace(Rank3, "D_2__", "")) %>% 
  mutate(Rank4 = str_replace(Rank4, "D_3__", "")) %>% 
  mutate(Rank5 = str_replace(Rank5, "D_4__", "")) %>% 
  mutate(Rank6 = str_replace(Rank6, "D_5__", ""))
c <- b[,-7]
c <- as.matrix(c)
d <- data.frame(otu_table(biom_qiime_all))

OTUq=otu_table(d, taxa_are_rows=TRUE)

TAXq=tax_table(c)

taxa_names(TAXq)=row.names(OTUq)

physeq_qiime_all = merge_phyloseq(OTUq, TAXq)

sample_names(physeq_qiime_all) <- paste("Q_188_", sample_names(physeq_qiime_all), sep="")

colnames(tax_table(physeq_qiime_all)) <- c("Kingdom", "Phylum", "Class",
                                           "Order","Family", "Genus")

#remove uneeded variables
rm(a, b, c, d, biom_qiime_30, biom_qiime_60, biom_qiime_120, biom_qiime_all,
   OTUq, TAXq)

#import mothur
#subthirty
physeq_mothur_30 <- import_mothur(mothur_shared_file = "mothur biom files/mothur.thirty.otu.shared", 
              mothur_constaxonomy_file = "mothur biom files/mothur.thirty.tax.taxonomy")

sample_names(physeq_mothur_30) <- paste("M_30_", sample_names(physeq_mothur_30), sep="")

colnames(tax_table(physeq_mothur_30)) <- c("Kingdom", "Phylum", "Class", 
                                           "Order", "Family",  "Genus")

#import mothur
#subsixty
physeq_mothur_60 <- import_mothur(mothur_shared_file = "mothur biom files/mothur.sixty.otu.shared", 
                                  mothur_constaxonomy_file = "mothur biom files/mothur.sixty.tax.taxonomy")

sample_names(physeq_mothur_60) <- paste("M_60_", sample_names(physeq_mothur_60), sep="")

colnames(tax_table(physeq_mothur_60)) <- c("Kingdom", "Phylum", "Class", 
                                           "Order", "Family",  "Genus")

#import mothur
#subonetwenty
physeq_mothur_120 <- import_mothur(mothur_shared_file = "mothur biom files/mothur.onetwenty.otu.shared", 
                                  mothur_constaxonomy_file = "mothur biom files/mothur.onetwenty.tax.taxonomy")

sample_names(physeq_mothur_120) <- paste("M_120_", sample_names(physeq_mothur_120), sep="")

colnames(tax_table(physeq_mothur_120)) <- c("Kingdom", "Phylum", "Class", 
                                            "Order", "Family",  "Genus")

#import mothur
#suball
physeq_mothur_all <- import_mothur(mothur_shared_file = "mothur biom files/mothur.all.otu.shared", 
                                  mothur_constaxonomy_file = "mothur biom files/mothur.all.tax.taxonomy")

sample_names(physeq_mothur_all) <- paste("M_188_", sample_names(physeq_mothur_all), sep="")

colnames(tax_table(physeq_mothur_all)) <- c("Kingdom", "Phylum", "Class", 
                                            "Order", "Family",  "Genus")

#import MG-RAST
#subthirty
otufull=read.csv("MG-Rast biom files\\otu_mgrast_thirty.csv",header=TRUE)
taxa=as.matrix(read.csv("MG-Rast biom files\\taxa_thirty.csv"))

OTU=otu_table(otufull, taxa_are_rows=TRUE)

TAX=tax_table(taxa)

taxa_names(TAX)=row.names(OTU)

sample_names(OTU) <- paste("G_30_", sample_names(OTU), sep="")

physeq_MG_30=phyloseq(OTU,TAX)
colnames(tax_table(physeq_MG_30)) <- c("Kingdom", "Phylum", "Class", 
                                       "Order", "Family",  "Genus")
rm(otufull, taxa, OTU, TAX)

#import MG-RAST
#subsixty
otufull=read.csv("MG-Rast biom files\\otu_mgrast_sixty.csv",header=TRUE)
taxa=as.matrix(read.csv("MG-Rast biom files\\taxa_sixty.csv"))

OTU=otu_table(otufull, taxa_are_rows=TRUE)

TAX=tax_table(taxa)

taxa_names(TAX)=row.names(OTU)

sample_names(OTU) <- paste("G_60_", sample_names(OTU), sep="")

physeq_MG_60=phyloseq(OTU,TAX)
colnames(tax_table(physeq_MG_60)) <- c("Kingdom", "Phylum", "Class", 
                                       "Order", "Family",  "Genus")
rm(otufull, taxa, OTU, TAX)

#import MG-RAST
#subonetwenty
otufull=read.csv("MG-Rast biom files\\otu_mgrast_onetwenty.csv",header=TRUE)
taxa=as.matrix(read.csv("MG-Rast biom files\\taxa_onetwenty.csv"))

OTU=otu_table(otufull, taxa_are_rows=TRUE)

TAX=tax_table(taxa)

taxa_names(TAX)=row.names(OTU)

sample_names(OTU) <- paste("G_120_", sample_names(OTU), sep="")

physeq_MG_120=phyloseq(OTU,TAX)
colnames(tax_table(physeq_MG_120)) <- c("Kingdom", "Phylum", "Class", 
                                        "Order", "Family",  "Genus")
rm(otufull, taxa, OTU, TAX)

#import MG-RAST
#all
otufull=read.csv("MG-Rast biom files\\otu_mgrast_all.csv",header=TRUE)
taxa=as.matrix(read.csv("MG-Rast biom files\\taxa_all.csv"))

OTU=otu_table(otufull, taxa_are_rows=TRUE)

TAX=tax_table(taxa)

taxa_names(TAX)=row.names(OTU)

sample_names(OTU) <- paste("G_188_", sample_names(OTU), sep="")

physeq_MG_all=phyloseq(OTU,TAX)
colnames(tax_table(physeq_MG_all)) <- c("Kingdom", "Phylum", "Class", 
                                        "Order", "Family",  "Genus")
rm(otufull, taxa, OTU, TAX)

#merging
merge_qiime = merge_phyloseq(physeq_qiime_30, physeq_qiime_60, 
                             physeq_qiime_120, physeq_qiime_all)
merge_mothur = merge_phyloseq(physeq_mothur_30, physeq_mothur_60,
                              physeq_mothur_120, physeq_mothur_all)
merge_MG = merge_phyloseq(physeq_MG_30, physeq_MG_60, physeq_MG_120,
                          physeq_MG_all)
merge = merge_phyloseq(merge_mothur, merge_MG, merge_qiime)

rm(physeq_qiime_30, physeq_qiime_60, physeq_qiime_120, physeq_qiime_all,
   physeq_mothur_30, physeq_mothur_60, physeq_mothur_120, physeq_mothur_all,
   physeq_MG_30, physeq_MG_60, physeq_MG_120, physeq_MG_all)

#tidying
x <- data.frame(otu_table(merge))
y <- data.frame(tax_table(merge))
x$names <- rownames(x)
y$names <- rownames(y)

z <- merge(x,y)
z <- z[,-1]

tax_group <- z %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))

#fixed taxa by hand
write.csv(tax_group, "tax_group.csv")

#fix it by hand
tax_group <- read.csv("tax_group.csv")
tax_group <- tax_group %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
head(tax_group)

tax_all <- as.matrix(tax_group[,1:6])

otu_all <- tax_group[,-1:-7]

OTU=otu_table(otu_all, taxa_are_rows=TRUE)

TAX=tax_table(tax_all)

taxa_names(TAX)=row.names(OTU)

physeq_all=phyloseq(OTU,TAX)

rm(otu_all, tax_all, tax_group, OTU, TAX)
