rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/ThesisProject_AKP/")

#packages
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(lme4)
library(microbiome)
library(phyloseq)
library(plyr)
library(PMCMR)
library(RVAideMemoire)
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

#combos
physeq_norm <- normalize_wout_rarefying(physeq)

physeq_norm_rec <- subset_samples(physeq_norm, Sample_Area == 'Rectum')
physeq_norm_ear <- subset_samples(physeq_norm, Sample_Area == 'Ears')
physeq_norm_eye <- subset_samples(physeq_norm, Sample_Area == 'Eyes')
physeq_norm_mou <- subset_samples(physeq_norm, Sample_Area == 'Mouth')
physeq_norm_nos <- subset_samples(physeq_norm, Sample_Area == 'Nose')

diver <- as.data.frame(divergence(physeq_norm, method = 'jaccard'))
diver_rec <- as.data.frame(divergence(physeq_norm_rec, method = 'jaccard'))
diver_ear <- as.data.frame(divergence(physeq_norm_ear, method = 'jaccard'))
diver_eye <- as.data.frame(divergence(physeq_norm_eye, method = 'jaccard'))
diver_mou <- as.data.frame(divergence(physeq_norm_mou, method = 'jaccard'))
diver_nos <- as.data.frame(divergence(physeq_norm_nos, method = 'jaccard'))


# all ---------------------------------------------------------------------


diver$SampleID = rownames(diver)
colnames(diver) <- c('Divergence', 'SampleID')
diver = merge(diver, metadata, by= 'SampleID')
diver <- diver %>% select(SampleID, Divergence, MoD, CoD_Simple, Violent,
                        Sample_Area, FinePMI)

kruskal.test(Divergence ~ Sample_Area, data = diver)
out <- posthoc.kruskal.nemenyi.test(x=diver$Divergence, g=diver$Sample_Area, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Divergence ~ MoD, data = diver)
out <- posthoc.kruskal.nemenyi.test(x=diver$Divergence, g=diver$MoD, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Divergence ~ CoD_Simple, data = diver)
out <- posthoc.kruskal.nemenyi.test(x=diver$Divergence, g=diver$CoD_Simple, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Divergence ~ Violent, data = diver)



# rectum ------------------------------------------------------------------


diver_rec$SampleID = rownames(diver_rec)
colnames(diver_rec) <- c('Divergence', 'SampleID')
diver_rec = merge(diver_rec, metadata, by= 'SampleID')
diver_rec <- diver_rec %>% select(SampleID, Divergence, MoD, CoD_Simple, Violent, Sample_Area)

kruskal.test(Divergence ~ MoD, data = diver_rec)
out <- posthoc.kruskal.nemenyi.test(x=diver_rec$Divergence, g=diver_rec$MoD, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Divergence ~ CoD_Simple, data = diver_rec)
out <- posthoc.kruskal.nemenyi.test(x=diver_rec$Divergence, g=diver_rec$CoD_Simple, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Divergence ~ Violent, data = diver_rec)



# eyes --------------------------------------------------------------------


diver_eye$SampleID = rownames(diver_eye)
colnames(diver_eye) <- c('Divergence', 'SampleID')
diver_eye = merge(diver_eye, metadata, by= 'SampleID')
diver_eye <- diver_eye %>% select(SampleID, Divergence, MoD, CoD_Simple, Violent,  Sample_Area)

kruskal.test(Divergence ~ MoD, data = diver_eye)
out <- posthoc.kruskal.nemenyi.test(x=diver_eye$Divergence, g=diver_eye$MoD, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Divergence ~ CoD_Simple, data = diver_eye)
out <- posthoc.kruskal.nemenyi.test(x=diver_eye$Divergence, g=diver_eye$CoD_Simple, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Divergence ~ Violent, data = diver_eye)


# ears --------------------------------------------------------------------

diver_ear$SampleID = rownames(diver_ear)
colnames(diver_ear) <- c('Divergence', 'SampleID')
diver_ear = merge(diver_ear, metadata, by= 'SampleID')
diver_ear <- diver_ear %>% select(SampleID, Divergence, MoD, CoD_Simple, Violent,  Sample_Area)

kruskal.test(Divergence ~ MoD, data = diver_ear)
out <- posthoc.kruskal.nemenyi.test(x=diver_ear$Divergence, g=diver_ear$MoD, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Divergence ~ CoD_Simple, data = diver_ear)
out <- posthoc.kruskal.nemenyi.test(x=diver_ear$Divergence, g=diver_ear$CoD_Simple, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Divergence ~ Violent, data = diver_ear)

# nose --------------------------------------------------------------------

diver_nos$SampleID = rownames(diver_nos)
colnames(diver_nos) <- c('Divergence', 'SampleID')
diver_nos = merge(diver_nos, metadata, by= 'SampleID')
diver_nos <- diver_nos %>% select(SampleID, Divergence, MoD, CoD_Simple, Violent,  Sample_Area)

kruskal.test(Divergence ~ MoD, data = diver_nos)
out <- posthoc.kruskal.nemenyi.test(x=diver_nos$Divergence, g=diver_nos$MoD, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Divergence ~ CoD_Simple, data = diver_nos)
out <- posthoc.kruskal.nemenyi.test(x=diver_nos$Divergence, g=diver_nos$CoD_Simple, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Divergence ~ Violent, data = diver_nos)



# mouth -------------------------------------------------------------------


diver_mou$SampleID = rownames(diver_mou)
colnames(diver_mou) <- c('Divergence', 'SampleID')
diver_mou = merge(diver_mou, metadata, by= 'SampleID')
diver_mou <- diver_mou %>% select(SampleID, Divergence, MoD, CoD_Simple, Violent,  Sample_Area)

kruskal.test(Divergence ~ MoD, data = diver_mou)
out <- posthoc.kruskal.nemenyi.test(x=diver_mou$Divergence, g=diver_mou$MoD, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Divergence ~ CoD_Simple, data = diver_mou)
out <- posthoc.kruskal.nemenyi.test(x=diver_mou$Divergence, g=diver_mou$CoD_Simple, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Divergence ~ Violent, data = diver_mou)


# combos ------------------------------------------------------------------

diver_combo <-rbind(diver_mou, diver_nos, diver_eye)

kruskal.test(Divergence ~ MoD, data = diver_combo)
out <- posthoc.kruskal.nemenyi.test(x=diver_combo$Divergence, g=diver_combo$MoD, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Divergence ~ CoD_Simple, data = diver_combo)
out <- posthoc.kruskal.nemenyi.test(x=diver_combo$Divergence, g=diver_combo$CoD_Simple, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Divergence ~ Violent, data = diver_combo)


diver_combo %>% group_by(MoD) %>% summarise_all(funs(mean, sd))
diver_combo %>% group_by(MoD) %>% summarise(no_rows = length(MoD))

diver_combo %>% group_by(CoD_Simple) %>% summarise_all(funs(mean, sd))
diver_combo %>% group_by(CoD_Simple) %>% summarise(no_rows = length(CoD_Simple))

diver_combo %>% group_by(Violent) %>% summarise_all(funs(mean, sd))
diver_combo %>% group_by(Violent) %>% summarise(no_rows = length(Violent))

theme_set(theme_classic(base_size = 25))
tiff("MoD_betadisp.TIF", width = 2000, height = 2000, res=300)
ggplot(diver_combo, aes(x=MoD, y=Divergence, color= MoD)) + geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2)) + xlab('Manner of Death') +
  ylab('Beta-dispersion (Jaccard distance)') + theme(legend.position = "none")
dev.off()

tiff("CoD_betadisp.TIF", width = 2000, height = 2000, res=300)
ggplot(diver_combo, aes(x=CoD_Simple, y=Divergence, color= CoD_Simple)) + geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2)) + xlab('Cause of Death') +
  ylab('Beta-dispersion (Jaccard distance)') + theme(legend.position = "none")
dev.off()

tiff("vio_betadisp.TIF", width = 2000, height = 2000, res=300)
ggplot(diver_combo, aes(x=Violent, y=Divergence, color= Violent)) + geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2)) + xlab('Violence') +
  ylab('Beta-dispersion (Jaccard distance)') + theme(legend.position = "none")
dev.off()

# simulated data ----------------------------------------------------------

l0 = lm( Divergence ~ 1, data=diver_combo)
l1 = lm(Divergence ~ Sample_Area, data = diver_combo)
l2 = lm(Divergence ~ MoD, data = diver_combo)
l3 = lm(Divergence ~ CoD_Simple, data = diver_combo)
l4 = lm(Divergence ~ Violent, data = diver_combo)
l5 = lm(Divergence ~ Sample_Area + MoD, data = diver_combo)
l6 = lm(Divergence ~ Sample_Area + CoD_Simple, data = diver_combo)
l7 = lm(Divergence ~ Sample_Area + Violent, data = diver_combo)
l8 = lm(Divergence ~ Sample_Area + CoD_Simple + Violent, data = diver_combo)
l9 = lm(Divergence ~ Sample_Area + MoD + Violent, data = diver_combo)
l10 = lm(Divergence ~ Sample_Area + MoD+ CoD_Simple + Violent, data = diver_combo)

anova(l0, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, test = 'Chisq')

summary(l5)

plot(resid(l6) ~ CoD_Simple, data = diver_combo)

pred.log.obs = rnorm(nrow(diver_combo), predict(l5), sd=sigma(l5))
pred.obs = exp(pred.log.obs)
plot(Divergence ~ MoD, data=diver_combo)
points(pred.log.obs ~ I(MoD+0.1), col=2, pch=20, data=diver_combo)
