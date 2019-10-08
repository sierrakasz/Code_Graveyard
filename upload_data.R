rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/Forensic_Pig_Combined/")

#packages
library(phyloseq)
library(plyr)
library(tidyverse)

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

physeq_qiime <- import_qiime_files("table-with-taxonomy.biom")
physeq_qiime
taxer <- data.frame(otu_table(physeq_qiime))
write.csv(taxer, "count_seqs.csv")

physeq_rare<- rarefy_even_depth(physeq_qiime, 
                                      rngseed = 711, sample.size = 5000)
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

tax_group <- tidying_data(physeq_rare)
write.csv(tax_group, "tax_group_fp_combo.csv")
