rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/ThesisProject_AKP/")

#packages
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(phyloseq)
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

physeq_new <- merge_phyloseq(physeq_ear_norm, physeq_eye_norm, physeq_mou_norm,
                             physeq_rec_norm, physeq_nos_norm)

GPrPhylum=tax_glom(physeq_new, "Phylum")
PhylumLevel = transform_sample_counts(GPrPhylum, function(x) x / sum(x))
PhylumLevel = filter_taxa(PhylumLevel, function(x) mean(x) > 0.0001, TRUE) 

df <- psmelt(PhylumLevel) 
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum", "Sample_Area"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
write.csv(Trtdata, "rel_abundance_phy.csv")

Topfive_trtdata <- filter(Trtdata, Phylum %in% c("Firmicutes", "Bacteroidetes",
                                                 "Actinobacteria", "Proteobacteria",
                                                 "Fusobacteria")) 

a <- ggplot(Topfive_trtdata, aes(x=Sample_Area, y=mean, fill=Phylum)) + 
  geom_bar(stat = 'identity', position = position_stack())
a

GPrFamily=tax_glom(physeq_new, "Family")
FamilyLevel = transform_sample_counts(GPrFamily, function(x) x / sum(x))
FamilyLevel = filter_taxa(FamilyLevel, function(x) mean(x) > 0.0001, TRUE) 

df <- psmelt(FamilyLevel) 
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Family", "Sample_Area"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)

write.csv(Trtdata, "rel_abundance_fam.csv")

Topten_trtdata <- filter(Trtdata, Family %in% c('Staphylococcaceae',
                                                'Carnobacteriaceae',
                                                'Corynebacteriaceae',
                                                'Streptococcaceae',
                                                'Moraxellaceae',
                                                'FamilyXI',
                                                'Prevotellaceae',
                                                'Enterobacteriaceae',
                                                'Pasteurellaceae',
                                                'Neisseriaceae',
                                                'Veillonellaceae',
                                                'Clostridiaceae1',
                                                'Micrococcaceae',
                                                'Fusobacteriaceae',
                                                'Lactobacillaceae',
                                                'Bacteroidaceae',
                                                'Lachnospiraceae',
                                                'Ruminococcaceae',
                                                'Porphyromonadaceae'
)) 

b <- ggplot(Topten_trtdata, aes(x=Sample_Area, y=mean, fill=Family)) + 
  geom_bar(stat = 'identity', position = position_stack())
b
