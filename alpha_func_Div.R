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

piphy <- import_biom("table-with-taxonomy.biom")

metadata=(read.csv("HPMMMetadata_pt2_cod.csv",header=TRUE))
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID

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
write.csv(rich, 'alpha_div_func.csv')


rich$Index <- factor(rich$Index, levels = c("Observed", "Shannon"))
p <- ggplot(rich, aes(x=Sample_Area, y=Observation, fill=Sample_Area)) +
  geom_boxplot(outlier.size = 3)
p + facet_wrap(~Index, scales= 'free') + xlab("Sample Area") + ylab("Counts") 

#stats
rich_obs <- rich %>% filter(Index == "Observed")
rich_sha <- rich %>% filter(Index == "Shannon")

rich_obs_rec <- rich_obs %>% filter(Sample_Area == "Rectum")
rich_obs_mou <- rich_obs %>% filter(Sample_Area == "Mouth")
rich_obs_eye <- rich_obs %>% filter(Sample_Area == "Eyes")
rich_obs_ear <- rich_obs %>% filter(Sample_Area == "Ears")
rich_obs_nos <- rich_obs %>% filter(Sample_Area == "Nose")

rich_sha_rec <- rich_sha %>% filter(Sample_Area == "Rectum")
rich_sha_mou <- rich_sha %>% filter(Sample_Area == "Mouth")
rich_sha_eye <- rich_sha %>% filter(Sample_Area == "Eyes")
rich_sha_ear <- rich_sha %>% filter(Sample_Area == "Ears")
rich_sha_nos <- rich_sha %>% filter(Sample_Area == "Nose")

rich_list <- list(rich_obs_mou, rich_obs_rec, rich_obs_ear, rich_obs_eye, rich_obs_nos,
                  rich_sha_mou, rich_sha_rec, rich_sha_ear, rich_sha_eye, rich_sha_nos)

for(i in rich_list) {
  print(mean(i$Observation))
  print(sd(i$Observation))
}

#kruskal wallis test
#nemenyi  test

#kruskal.test(alpha_div~pipeline, data = alpha)
#out <- posthoc.kruskal.nemenyi.test(x=alpha_div, g=pipeline, dist="Tukey", p.adjust.method = 'bonf')
#print(otu$statistic)
#https://cran.r-project.org/web/packages/PMCMR/vignettes/PMCMR.pdf



kruskal.test(Observation ~ Sample_Area, data = rich_obs)
out <- posthoc.kruskal.nemenyi.test(x=rich_obs$Observation, g=rich_obs$Sample_Area, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

kruskal.test(Observation ~ Sample_Area, data = rich_sha)
out <- posthoc.kruskal.nemenyi.test(x=rich_sha$Observation, g=rich_sha$Sample_Area, dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)


