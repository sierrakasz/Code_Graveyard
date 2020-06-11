rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/ThesisProject_AKP/")

#packages
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(microbiome)
library(phyloseq)
library(dplyr)
library(PMCMR)
library(tidyverse)
library(vegan)

#upload the data
tax_group <- read.csv("tax_group_AKP_edit.csv")
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

GPfr_rec_24 <-subset_samples(physeq_rec_norm, FinePMI == 'Less24')
GPfr_eye_24 <-subset_samples(physeq_eye_norm, FinePMI == 'Less24')
GPfr_ear_24 <-subset_samples(physeq_ear_norm, FinePMI == 'Less24')
GPfr_nos_24 <-subset_samples(physeq_nos_norm, FinePMI == 'Less24')
GPfr_mou_24 <-subset_samples(physeq_mou_norm, FinePMI == 'Less24')

lets_find_the_core <- function(physeq) {
  otu <- data.frame(otu_table(physeq))
  tax <- data.frame(tax_table(physeq))
  sample_nums <- rowSums(otu == 0)
  otu$OTUID <- row.names(otu)
  tax$OTUID <- row.names(tax)
  len <- length(colnames(otu)) -1
  core_perc <- sample_nums / len
  new_otu <- cbind(otu, core_perc)
  new_core_taxa <- new_otu$OTUID[new_otu$core_perc <= 0.10]
  new_core_otu <- filter(new_otu, new_otu$OTUID %in% new_core_taxa)
  new_core_otu <- subset(new_core_otu, select=-c(OTUID, core_perc))
  new_core_tax <- filter(tax, tax$OTUID %in% new_core_taxa)
  return(new_core_taxa)
}

core_list_rec <- lets_find_the_core(GPfr_rec_24)
core_list_eye <- lets_find_the_core(GPfr_eye_24)
core_list_ear <- lets_find_the_core(GPfr_ear_24)
core_list_nos <- lets_find_the_core(GPfr_nos_24)
core_list_mou <- lets_find_the_core(GPfr_mou_24)

physeq_ear_norm_Core <- prune_taxa(core_list_ear_tax, physeq_ear_norm)
physeq_eye_norm_Core <- prune_taxa(core_list_eye_tax, physeq_eye_norm)
physeq_nos_norm_Core <- prune_taxa(core_list_nos_tax, physeq_nos_norm)
physeq_mou_norm_Core <- prune_taxa(core_list_mou_tax, physeq_mou_norm)
physeq_rec_norm_Core <- prune_taxa(core_list_rec_tax, physeq_rec_norm)

noncore_physeq <- function(physeq, core_list) {
  otu <- data.frame(otu_table(physeq))
  tax <- data.frame(tax_table(physeq))
  otu$OTUID <- row.names(otu)
  tax$OTUID <- row.names(tax)
  dis_core_tax <- filter(tax, !tax$OTUID %in% core_list)
  dis_core_otu <- dis_core_tax$OTUID
  newphyseq <- prune_taxa(dis_core_otu, physeq)
  return(newphyseq)
}

physeq_ear_norm_non <- noncore_physeq(physeq_ear_norm, core_list_ear_tax)
physeq_eye_norm_non <- noncore_physeq(physeq_eye_norm, core_list_eye_tax)
physeq_rec_norm_non <- noncore_physeq(physeq_rec_norm, core_list_rec_tax)
physeq_mou_norm_non <- noncore_physeq(physeq_mou_norm, core_list_mou_tax)
physeq_nos_norm_non <- noncore_physeq(physeq_nos_norm, core_list_nos_tax)

aa <- divergence(physeq_ear_norm_Core, method = 'jaccard')
bb <- divergence(physeq_eye_norm_Core, method = 'jaccard')
cc <- divergence(physeq_rec_norm_Core, method = 'jaccard')
dd <- divergence(physeq_nos_norm_Core, method = 'jaccard')
ee <- divergence(physeq_mou_norm_Core, method = 'jaccard')

beta_div_box_plot_core <- data.frame(c(aa,bb,cc,dd,ee), 
                                     c(rep('Ears', times = length(aa)), rep('Eyes', times = length(bb)),
                                       rep('Rectum', times = length(cc)), rep('Nose', times = length(dd)),
                                       rep('Mouth', times = length(ee))))

colnames(beta_div_box_plot_core) <- c('Jaccard_distance', 'Body_Site')
beta_div_box_plot_core$Core <- 'Core'

aaa <- divergence(physeq_ear_norm_non, method = 'jaccard')
bbb <- divergence(physeq_eye_norm_non, method = 'jaccard')
ccc <- divergence(physeq_rec_norm_non, method = 'jaccard')
ddd <- divergence(physeq_nos_norm_non, method = 'jaccard')
eee <- divergence(physeq_mou_norm_non, method = 'jaccard')

beta_div_box_plot_non <- data.frame(c(aaa,bbb,ccc,ddd,eee), 
                                    c(rep('Ears', times = length(aaa)), rep('Eyes', times = length(bbb)),
                                      rep('Rectum', times = length(ccc)), rep('Nose', times = length(ddd)),
                                      rep('Mouth', times = length(eee))))

colnames(beta_div_box_plot_non) <- c('Jaccard_distance', 'Body_Site')
beta_div_box_plot_non$Core <- 'Non-core'

beta_div_taxa_plots <- rbind(beta_div_box_plot_core, beta_div_box_plot_non)

a <- ggplot(data = beta_div_taxa_plots, aes(x=Core, y = Jaccard_distance, fill=Body_Site)) +
  geom_boxplot() + facet_wrap(~Body_Site, nrow = 1) + ylab("Jaccard Distances") + 
  xlab("Core or Non-core Taxa") +
  scale_fill_manual("Body Site", values = c('Ears' = '#ECA72F', 'Eyes' = '#3DC5B6', 
                                            'Nose' = '#D45924', 'Mouth' = '#F65E5A', 'Rectum' = '#90839F')) +
  theme(axis.text.x = element_text(angle = 90))
a

ddply(beta_div_taxa_plots,~Core,summarise,mean=mean(Jaccard_distance), sd=sd(Jaccard_distance))
ddply(beta_div_taxa_plots,~Body_Site,summarise,mean=mean(Jaccard_distance), sd=sd(Jaccard_distance))

kruskal.test(Jaccard_distance ~ as.factor(Core), data = beta_div_taxa_plots)

kruskal.test(Jaccard_distance ~ as.factor(Body_Site), data = beta_div_taxa_plots)
out <- posthoc.kruskal.nemenyi.test(x=beta_div_taxa_plots$Jaccard_distance, 
                                    g=as.factor(beta_div_taxa_plots$Body_Site), 
                                    dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

beta_div_box_plot_core$SampleID <- rownames(beta_div_box_plot_core)
beta_div_box_plots_meta_c <- merge(metadata, beta_div_box_plot_core, by  = 'SampleID')
beta_div_box_plot_non$SampleID <- rownames(beta_div_box_plot_non)
beta_div_box_plots_meta_n <- merge(metadata, beta_div_box_plot_non, by  = 'SampleID')
beta_div_pmi_taxa <- rbind(beta_div_box_plots_meta_c, beta_div_box_plots_meta_n)

b <- ggplot(data = beta_div_pmi_taxa, aes(x=Core, y = Jaccard_distance, fill=FinePMI)) +
  geom_boxplot() + facet_wrap(~FinePMI, nrow = 1) + ylab("Jaccard Distances") + 
  xlab("Core or Non-core Taxa") +
  theme(axis.text.x = element_text(angle = 90))
b

ddply(beta_div_pmi_taxa,~FinePMI,summarise,mean=mean(Jaccard_distance), sd=sd(Jaccard_distance))

kruskal.test(Jaccard_distance ~ as.factor(FinePMI), data = beta_div_pmi_taxa)
out <- posthoc.kruskal.nemenyi.test(x=beta_div_pmi_taxa$Jaccard_distance, 
                                    g=as.factor(beta_div_pmi_taxa$FinePMI), 
                                    dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

#functions
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

GPfr_rec_24 <-subset_samples(piphy_rec_norm, FinePMI == 'Less24')
GPfr_eye_24 <-subset_samples(piphy_eye_norm, FinePMI == 'Less24')
GPfr_ear_24 <-subset_samples(piphy_ear_norm, FinePMI == 'Less24')
GPfr_nos_24 <-subset_samples(piphy_nos_norm, FinePMI == 'Less24')
GPfr_mou_24 <-subset_samples(piphy_mou_norm, FinePMI == 'Less24')

lets_find_the_core <- function(physeq) {
  otu <- data.frame(otu_table(physeq))
  tax <- data.frame(tax_table(physeq))
  sample_nums <- rowSums(otu == 0)
  otu$OTUID <- row.names(otu)
  tax$OTUID <- row.names(tax)
  len <- length(colnames(otu)) -1
  core_perc <- sample_nums / len
  new_otu <- cbind(otu, core_perc)
  new_core_taxa <- new_otu$OTUID[new_otu$core_perc <= 0.10]
  new_core_otu <- filter(new_otu, new_otu$OTUID %in% new_core_taxa)
  new_core_otu <- subset(new_core_otu, select=-c(OTUID, core_perc))
  new_core_tax <- filter(tax, tax$OTUID %in% new_core_taxa)
  return(new_core_taxa)
}

core_list_rec <- lets_find_the_core(GPfr_rec_24)
core_list_eye <- lets_find_the_core(GPfr_eye_24)
core_list_ear <- lets_find_the_core(GPfr_ear_24)
core_list_nos <- lets_find_the_core(GPfr_nos_24)
core_list_mou <- lets_find_the_core(GPfr_mou_24)


piphy_ear_norm_Core <- prune_taxa(core_list_ear, piphy_ear_norm)
piphy_eye_norm_Core <- prune_taxa(core_list_eye, piphy_eye_norm)
piphy_nos_norm_Core <- prune_taxa(core_list_nos, piphy_nos_norm)
piphy_mou_norm_Core <- prune_taxa(core_list_mou, piphy_mou_norm)
piphy_rec_norm_Core <- prune_taxa(core_list_rec, piphy_rec_norm)

noncore_physeq <- function(physeq, core_list) {
  otu <- data.frame(otu_table(physeq))
  tax <- data.frame(tax_table(physeq))
  otu$OTUID <- row.names(otu)
  tax$OTUID <- row.names(tax)
  dis_core_tax <- filter(tax, !tax$OTUID %in% core_list)
  dis_core_otu <- dis_core_tax$OTUID
  newphyseq <- prune_taxa(dis_core_otu, physeq)
  return(newphyseq)
}

piphy_ear_norm_non <- noncore_physeq(piphy_ear_norm, core_list_ear)
piphy_eye_norm_non <- noncore_physeq(piphy_eye_norm, core_list_eye)
piphy_rec_norm_non <- noncore_physeq(piphy_rec_norm, core_list_rec)
piphy_mou_norm_non <- noncore_physeq(piphy_mou_norm, core_list_mou)
piphy_nos_norm_non <- noncore_physeq(piphy_nos_norm, core_list_nos)

af <- divergence(piphy_ear_norm_Core, method = 'jaccard')
bf <- divergence(piphy_eye_norm_Core, method = 'jaccard')
cf <- divergence(piphy_rec_norm_Core, method = 'jaccard')
df <- divergence(piphy_nos_norm_Core, method = 'jaccard')
ef <- divergence(piphy_mou_norm_Core, method = 'jaccard')

beta_div_box_plot_core_func <- data.frame(c(af,bf,cf,df,ef), 
                                          c(rep('Ears', times = length(af)), rep('Eyes', times = length(bf)),
                                            rep('Rectum', times = length(cf)), rep('Nose', times = length(df)),
                                            rep('Mouth', times = length(ef))))

colnames(beta_div_box_plot_core_func) <- c('Jaccard_distance', 'Body_Site')
beta_div_box_plot_core_func$Core <- 'Core'

aaf <- divergence(piphy_ear_norm_non, method = 'jaccard')
bbf <- divergence(piphy_eye_norm_non, method = 'jaccard')
ccf <- divergence(piphy_rec_norm_non, method = 'jaccard')
ddf <- divergence(piphy_nos_norm_non, method = 'jaccard')
eef <- divergence(piphy_mou_norm_non, method = 'jaccard')

beta_div_box_plot_non_func <- data.frame(c(aaf,bbf,ccf,ddf,eef), 
                                         c(rep('Ears', times = length(aaf)), rep('Eyes', times = length(bbf)),
                                           rep('Rectum', times = length(ccf)), rep('Nose', times = length(ddf)),
                                           rep('Mouth', times = length(eef))))

colnames(beta_div_box_plot_non_func) <- c('Jaccard_distance', 'Body_Site')
beta_div_box_plot_non_func$Core <- 'Non-core'

beta_div_taxa_plots_func <- rbind(beta_div_box_plot_core_func, beta_div_box_plot_non_func)

c <- ggplot(data = beta_div_taxa_plots_func, aes(x=Core, y = Jaccard_distance, fill=Body_Site)) +
  geom_boxplot() + facet_wrap(~Body_Site, nrow = 1) + ylab("Jaccard Distances") + 
  xlab("Core or Non-core Functions") +
  scale_fill_manual("Body Site", values = c('Ears' = '#ECA72F', 'Eyes' = '#3DC5B6', 
                                            'Nose' = '#D45924', 'Mouth' = '#F65E5A', 'Rectum' = '#90839F')) +
  theme(axis.text.x = element_text(angle = 90))
c

ddply(beta_div_taxa_plots_func,~Core,summarise,mean=mean(Jaccard_distance), sd=sd(Jaccard_distance))
ddply(beta_div_taxa_plots_func,~Body_Site,summarise,mean=mean(Jaccard_distance), sd=sd(Jaccard_distance))

kruskal.test(Jaccard_distance ~ as.factor(Core), data = beta_div_taxa_plots_func)

kruskal.test(Jaccard_distance ~ as.factor(Body_Site), data = beta_div_taxa_plots_func)
out <- posthoc.kruskal.nemenyi.test(x=beta_div_taxa_plots_func$Jaccard_distance, 
                                    g=as.factor(beta_div_taxa_plots_func$Body_Site), 
                                    dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

beta_div_box_plot_core_func$SampleID <- rownames(beta_div_box_plot_core_func)
beta_div_box_plots_meta_cf <- merge(metadata, beta_div_box_plot_core_func, by  = 'SampleID')
beta_div_box_plot_non_func$SampleID <- rownames(beta_div_box_plot_non_func)
beta_div_box_plots_meta_nf <- merge(metadata, beta_div_box_plot_non_func, by  = 'SampleID')
beta_div_pmi_func <- rbind(beta_div_box_plots_meta_cf, beta_div_box_plots_meta_nf)

d <- ggplot(data = beta_div_pmi_func, aes(x=Core, y = Jaccard_distance, fill=FinePMI)) +
  geom_boxplot() + facet_wrap(~FinePMI, nrow = 1) + ylab("Jaccard Distances") + 
  xlab("Core or Non-core Taxa") +
  theme(axis.text.x = element_text(angle = 90))
d

ddply(beta_div_pmi_taxa,~FinePMI,summarise,mean=mean(Jaccard_distance), sd=sd(Jaccard_distance))

kruskal.test(Jaccard_distance ~ as.factor(FinePMI), data = beta_div_pmi_taxa)
out <- posthoc.kruskal.nemenyi.test(x=beta_div_pmi_taxa$Jaccard_distance, 
                                    g=as.factor(beta_div_pmi_taxa$FinePMI), 
                                    dist='Tukey', p.adjust.method = 'bonf' )
print(out$p.value)

#Final Figure
theme_set(theme_classic(base_size = 18))
tiff("Fig4.TIF", width = 6000, height = 4000, res=300)
ggarrange(a,b,c,d, 
          labels = c("A", "B", "C", "D"),
          nrow = 2, ncol = 2)
dev.off()