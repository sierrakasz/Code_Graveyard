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

Topfive_trtdata <- filter(Trtdata, Phylum %in% c("Firmicutes", "Bacteroidetes",
                                                 "Actinobacteria", "Proteobacteria",
                                                 "Fusobacteria")) 

a <- ggplot(Topfive_trtdata, aes(x=Sample_Area, y=mean, fill=Phylum)) + 
  geom_bar(stat = 'identity', position = position_stack()) +
  scale_fill_manual("Phyla", values = c("Actinobacteria" = "#B65143", "Bacteroidetes" = "#E47725", 
                                         "Firmicutes" = "#46CF92", "Fusobacteria" = '#2C7BC2', 
                                         "Proteobacteria" = '#523D96')) + ylab("Relative Abundance") +
  xlab("Body Site") 
a

#Alpha Diversity
erich_ey <- estimate_richness(physeq_eye_norm, measures = c("Observed", "Shannon"))
erich_ea <- estimate_richness(physeq_ear_norm, measures = c("Observed", "Shannon"))
erich_m <- estimate_richness(physeq_mou_norm, measures = c("Observed", "Shannon"))
erich_n <- estimate_richness(physeq_nos_norm, measures = c("Observed", "Shannon"))
erich_r <- estimate_richness(physeq_rec_norm, measures = c("Observed", "Shannon"))

erich <- rbind(erich_ey, erich_ea, erich_m, erich_n, erich_r)

erich <- add_rownames(erich, "SampleID")

#make data tidy 
erich2 <- erich %>%
  gather(Index, Observation, c("Observed", "Shannon"), na.rm = TRUE)

#add metadata
rich = merge(erich2, metadata)
rich <- rich %>% select(SampleID, Index, Observation, Sample_Area)

rich$Index <- factor(rich$Index, levels = c("Observed", "Shannon"))
p <- ggplot(rich, aes(x=Sample_Area, y=Observation, fill=Sample_Area)) +
  geom_boxplot(outlier.size = 3)
b <- p + facet_wrap(~Index, scales= 'free') + xlab("Body Site") + ylab("Counts") +
  scale_fill_manual("Body Site", values = c('Ears' = '#ECA72F', 'Eyes' = '#3DC5B6', 
                    'Nose' = '#D45924', 'Mouth' = '#F65E5A', 'Rectum' = '#90839F')) +
  theme(axis.text.x=element_blank())
b


# funciton

piphy <- import_biom("table-with-taxonomy.biom")

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

piphy_new <- merge_phyloseq(piphy_ear_norm, piphy_eye_norm, piphy_mou_norm,
                            piphy_rec_norm, piphy_nos_norm)

FuncLevel = transform_sample_counts(piphy_new, function(x) x / sum(x))
FuncLevel = filter_taxa(FuncLevel, function(x) mean(x) > 0.0001, TRUE) 

df <- psmelt(FuncLevel) 
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("OTU", 'Rank2', "Sample_Area"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)

Topfive_trtdata <- filter(Trtdata, OTU %in% c('K01990',
                                              'K01992',
                                              'K07024',
                                              'K06147',
                                              'K02015',
                                              'K01277',
                                              'K02004',
                                              'K03088')) 

c <- ggplot(Topfive_trtdata, aes(x=Sample_Area, y=mean, fill=OTU)) + 
  geom_bar(stat = 'identity', position = position_stack()) + xlab("Body Site") + ylab("Relative Abundance") +
  scale_fill_manual("Kegg ID", values = c('#3C8F42', '#8F683C', '#19A49F', '#36376B',
                                          '#D7AFD3', '#D2D56C', '#F27D43', '#91342E'))
c

piphy_new <- merge_phyloseq(piphy_ear_norm, piphy_eye_norm, piphy_mou_norm,
                            piphy_rec_norm, piphy_nos_norm)

FuncLevel = transform_sample_counts(piphy_new, function(x) round(x))

#Alpha Diversity
erich <- estimate_richness(FuncLevel, measures = c("Observed", "Shannon"))

erich <- add_rownames(erich, "SampleID")

#make data tidy 
erich2 <- erich %>%
  gather(Index, Observation, c("Observed", "Shannon"), na.rm = TRUE)

#add metadata
rich = merge(erich2, metadata)
rich <- rich %>% select(SampleID, Index, Observation, Sample_Area)

rich$Index <- factor(rich$Index, levels = c("Observed", "Shannon"))
l <- ggplot(rich, aes(x=Sample_Area, y=Observation, fill=Sample_Area)) +
  geom_boxplot(outlier.size = 3)
d <- l + facet_wrap(~Index, scales= 'free') + xlab("Body Site") + ylab("Counts") +
  scale_fill_manual("Body Site", values = c('Ears' = '#ECA72F', 'Eyes' = '#3DC5B6', 
                                            'Nose' = '#D45924', 'Mouth' = '#F65E5A', 'Rectum' = '#90839F')) + 
  theme(axis.text.x=element_blank())
d


#Final Figure
theme_set(theme_classic(base_size = 18))
tiff("Fig1.TIF", width = 4000, height = 4000, res=300)
ggarrange(a,b,c,d + rremove("x.text"), 
          labels = c("A", "B", "C", "D"),
          nrow = 2, ncol = 2)
dev.off()


