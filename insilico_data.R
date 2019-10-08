rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/ThesisProject_Pipelines/")

#packages
library(phyloseq)
library(plyr)
library(tidyverse)

#import mothur
import_mothur_files <- function(shared, tax) {
  physeq_mothur <- import_mothur(mothur_shared_file = shared,
                                 mothur_constaxonomy_file = tax)
  colnames(tax_table(physeq_mothur)) <- c("Kingdom", "Phylum", "Class", 
                                          "Order", "Family",  "Genus")
  return(physeq_mothur)
}

physeq_mothur <- import_mothur_files("insilico/mothur_mock3.shared", 
                                     "insilico/mothur_mock3.taxonomy")
sample_names(physeq_mothur) <- c('M_HMPMockV1.1.Even1','M_HMPMockV1.1.Even2',
                                 'M_HMPMockV1.2.Staggered1','M_HMPMockV1.2.Staggered2'
)

View(as.data.frame(tax_table(physeq_mothur)))

#import qiime
import_qiime_files <- function(biom) {
  biom_qiime <- import_biom(biom)
  a <- data.frame(tax_table(biom_qiime))
  b <- a %>% 
    mutate(Rank1 = str_replace(Rank1, "D_0__", "")) %>% 
    mutate(Rank2 = str_replace(Rank2, "D_1__", "")) %>% 
    mutate(Rank3 = str_replace(Rank3, "D_2__", "")) %>% 
    mutate(Rank4 = str_replace(Rank4, "D_3__", "")) %>% 
    mutate(Rank5 = str_replace(Rank5, "D_4__", "")) %>% 
    mutate(Rank6 = str_replace(Rank6, "D_5__", ""))
  c <- b[,c("Rank1", "Rank2", "Rank3", "Rank4", "Rank5", "Rank6")]
  c <- as.matrix(c)
  d <- data.frame(otu_table(biom_qiime))
  OTUq=otu_table(d, taxa_are_rows=TRUE)
  TAXq=tax_table(c)
  taxa_names(TAXq)=row.names(OTUq)
  physeq_qiime = merge_phyloseq(OTUq, TAXq)
  colnames(tax_table(physeq_qiime)) <- c("Kingdom", "Phylum", "Class",
                                         "Order","Family", "Genus")
  return(physeq_qiime)
}

physeq_qiime <- import_qiime_files("insilico/table-with-taxonomy.biom")
sample_names(physeq_qiime) <- paste("Q_", sample_names(physeq_qiime), sep="")

#import MG-RAST
import_MGRAST <- function(file) {
  file_name <- read.csv(file)
  otu <- select(file_name, contains('HMP'))
  tax <-select(file_name, domain, phylum, className, order, family, genus)
  tax <- as.matrix(tax)
  OTU=otu_table(otu, taxa_are_rows=TRUE)
  TAX=tax_table(tax)
  taxa_names(TAX)=row.names(OTU)
  physeq_MG=phyloseq(OTU,TAX)
  colnames(tax_table(physeq_MG)) <- c("Kingdom", "Phylum", "Class", 
                                      "Order", "Family",  "Genus")
  return(physeq_MG)
}
physeq_MG <- import_MGRAST("insilico\\Mock3_MG-RAST.csv")
sample_names(physeq_MG) <- paste("G_", sample_names(physeq_MG), sep="")

#merging phyloseq objects together
merge = merge_phyloseq(physeq_qiime,physeq_mothur,physeq_MG)

#tidying data to insure that taxa are consistent 
tidying_data <- function(physeq) {
  x <- data.frame(otu_table(physeq))
  y <- data.frame(tax_table(physeq))
  x$names <- rownames(x)
  y$names <- rownames(y)
  z <- merge(x,y)
  z <- z %>% select(-contains("names"))
  tax_group <- z %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
  tax_group <- tax_group %>%  ungroup()
  return(tax_group)
}
tax_group <- tidying_data(merge)
write.csv(tax_group, "tax_group_insilico.csv")
tax_group <- read.csv("tax_group_insilico.csv")
tax_group <- tax_group %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group <- tax_group %>%  ungroup()

regroup_physeq_object <-function(table) {
  tax <- table %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
  tax <- as.matrix(tax)
  otu <- table %>% select(contains("HMP"))
  OTU=otu_table(otu, taxa_are_rows=TRUE)
  TAX=tax_table(tax)
  taxa_names(TAX)=row.names(OTU)
  physeq_all=phyloseq(OTU,TAX)
  return(physeq_all)
}

physeq_all <- regroup_physeq_object(tax_group)

#import metadata and combine
metadata=read.csv("insilico/mock3_meta_merge.csv",header=TRUE)
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
physeq=merge_phyloseq(physeq_all, sampdat)

#Tax glom
GPrGenus=tax_glom(physeq, "Genus")
GPrGenus <- rarefy_even_depth(GPrGenus)
otu <- as.data.frame(otu_table(GPrGenus))
seq_num <- sum(otu[,1])

#realtive abundance by variable
df <- psmelt(GPrGenus) 
write.csv(df, "rel_abundance_insilico_asis.csv")
#fix names in excel sample - take out prefix
df <- read.csv("rel_abundance_insilico_asis.csv", header=T)
new_df <- df[,c("Sample", "Abundance", "Pipeline", "Kingdom", "Phylum", "Class", 
                "Order", "Family",  "Genus")]

#Linear regression for 
taxa_answers <- read.csv('insilico/taxa_answers.csv', header=T)

#Kingdom, Phylum, Class, Order, Family, Genus,
mo <- new_df %>% filter(Pipeline == 'MG-RAST', Abundance > 0)
mo_true <- taxa_answers %>%  select(Kingdom, Phylum, Class, Order, Family, Genus, Sample, Expected)

mo = subset(mo, select = -c(Pipeline))
mo <- mo %>% group_by(Sample, Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
mo <- mo %>%  ungroup()

mo_true <- mo_true %>% mutate(Expected = Expected*seq_num)
mo_true <- mo_true %>% group_by(Sample, Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
mo_true <- mo_true %>%  ungroup()

mo %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus)
mo_true %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus)

#check names
compare <- merge(mo, mo_true, all= T)
compare[is.na(compare)] <- 0

#percent error
compare <- compare %>% mutate(Error = abs(Abundance - Expected))
percent_error <- sum(compare$Error)/sum(compare$Expected)
percent_error

plot(compare$Expected, compare$Abundance)
                                
#Errors using linear model
fit <- lm(Abundance~Expected, data=compare)
summary(fit)

#MSE
mean(residuals(fit)^2)
#RMSE
sqrt(mean(residuals(fit)^2))
#RSS
sum(residuals(fit)^2)
#RSE
sqrt( sum(residuals(fit)^2) / fit$df.residual ) 


View(vec_ans)
vec_ans <- unique(taxa_answers$Genus)
vec <- unique(mo$Genus)

#pipeline added it
setdiff(vec, vec_ans)

#pipeline missed it
setdiff(vec_ans, vec)

sum(mo$Genus == 'Unclassified')

