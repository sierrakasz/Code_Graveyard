rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/Forensic_Pig_Combined/")

#packages
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(multcompView)
library(phyloseq)
library(plyr)
library(tidyverse)

#make a table in excel for rare data
tax_group <- read.csv("tax_group_fp_combo.csv")
tax_group <- tax_group %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group <- tax_group %>%  ungroup()

regroup_physeq_object <-function(table) {
  tax <- table %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
  tax <- as.matrix(tax)
  otu <- table %>% select(contains("P"))
  otu <- otu[,-1]
  OTU=otu_table(otu, taxa_are_rows=TRUE)
  TAX=tax_table(tax)
  taxa_names(TAX)=row.names(OTU)
  physeq_all=phyloseq(OTU,TAX)
  return(physeq_all)
}

physeq <- regroup_physeq_object(tax_group)

metadata=(read.csv("AquaticPigMetadataCombined.csv",header=TRUE))
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
sampdat$Decomp_Stage <- factor(sampdat$Decomp_Stage, levels = c('Submerged_Fresh', 'Early_Floating',
                                                                'Floating_Decay', 'Advanced_Floating_Decay', 
                                                                'Sunken_Remains'))

physeq=merge_phyloseq(physeq, sampdat)
physeq_benbow <- subset_samples(physeq, Study == "Benbow")
physeq_wallace <- subset_samples(physeq, Study == "Wallace")

#stacked bar graph
relative_abundance_setup <- function(physeq) {
  GPrPhylum=tax_glom(physeq, "Phylum")
  PhylumLevel_rec = transform_sample_counts(GPrPhylum, function(x) x / sum(x))
  PhylumLevel_rec = filter_taxa(PhylumLevel_rec, function(x) mean(x) > 0.005, TRUE) 
  df <- psmelt(PhylumLevel_rec) 
  df$Abundance=df$Abundance*100
  Trtdata <- ddply(df, c("Phylum", 'Study', 'Decomp_Stage'), summarise,
                   N    = length(Abundance),
                   mean = mean(Abundance),
                   sd   = sd(Abundance),
                   se   = sd / sqrt(N)
  )
}

df_graph <- relative_abundance_setup(physeq)
study_names <- list('Benbow' = 'Freshwater', "Wallace"= "Estuarine")
study_labeller <- function(variable,value){
  return(study_names[value])
}
theme_set(theme_bw(base_size = 30))
tiff("rel_abund_both.TIF", width = 1200, height = 900)
cdataplot=ggplot(df_graph, aes(x=Decomp_Stage,y=mean))+
  geom_bar(aes(fill = Phylum),colour="black", stat="identity")+
  facet_wrap(~Study, labeller=study_labeller) +
  xlab("Decomposition Stage")+ylab("Relative Abundance (> 5%)") + 
  scale_x_discrete(labels=c('Submerged_Fresh' ='SF', 'Early_Floating'='EF',
                            'Floating_Decay'='FD', 'Advanced_Floating_Decay'='AFD', 
                            'Sunken_Remains'='SR')) +
  theme(panel.grid = element_blank())
cdataplot
dev.off()

#family level
#stacked bar graph
relative_abundance_setup <- function(physeq) {
  GPrFamily=tax_glom(physeq, "Family")
  FamilyLevel = transform_sample_counts(GPrFamily, function(x) x / sum(x))
  FamilyLevel = filter_taxa(FamilyLevel, function(x) mean(x) > 0.005, TRUE) 
  df <- psmelt(FamilyLevel) 
  df$Abundance=df$Abundance*100
  Trtdata <- ddply(df, c("Family", 'Study', 'Decomp_Stage'), summarise,
                   N    = length(Abundance),
                   mean = mean(Abundance),
                   sd   = sd(Abundance),
                   se   = sd / sqrt(N)
  )
}

df_graph_fam <- relative_abundance_setup(physeq)
write.csv(df_graph_fam, "family_rel_abund.csv")




relative_abundance_setup <- function(physeq) {
  GPrPhylum=tax_glom(physeq, "Phylum")
  PhylumLevel_rec = transform_sample_counts(GPrPhylum, function(x) x / sum(x))
  PhylumLevel_rec = filter_taxa(PhylumLevel_rec, function(x) mean(x) > 0.001, TRUE) 
  df <- psmelt(PhylumLevel_rec) 
  df$Abundance=df$Abundance*100
  Trtdata <- ddply(df, c("Phylum", 'Study'), summarise,
                   N    = length(Abundance),
                   mean = mean(Abundance),
                   sd   = sd(Abundance),
                   se   = sd / sqrt(N)
  )
}

df_graph <- relative_abundance_setup(physeq)
















GPrPhylum=tax_glom(physeq_benbow, "Phylum")
PhylumLevel = transform_sample_counts(GPrPhylum, function(x) x / sum(x))
PhylumLevel = filter_taxa(PhylumLevel, function(x) mean(x) > 0.005, TRUE) 

df <- psmelt(PhylumLevel)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum", "Decomp_Stage"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
write.csv(Trtdata, "benbow_rel_abund_total.csv")
theme_set(theme_bw(base_size = 24))
tiff("rel_abund_benbow.TIF", width = 800, height = 800)
cdataplot=ggplot(Trtdata, aes(x=Decomp_Stage,y=mean))+
  geom_bar(aes(fill = Phylum),colour="black", stat="identity")+
  xlab("Decomposition Stage")+ylab("Relative Abundance (> 1%)") + 
  scale_x_discrete(labels=c('Submerged_Fresh' ='SF', 'Early_Floating'='EF',
                              'Floating_Decay'='FD', 'Advanced_Floating_Decay'='AFD', 
                              'Sunken_Remains'='SR'))
cdataplot
dev.off()

GPrPhylumw=tax_glom(physeq_wallace, "Phylum")
PhylumLevelw = transform_sample_counts(GPrPhylumw, function(x) x / sum(x))
PhylumLevelw = filter_taxa(PhylumLevelw, function(x) mean(x) > 0.005, TRUE) 

dfw <- psmelt(PhylumLevelw)
dfw$Abundance=dfw$Abundance*100
Trtdataw <- ddply(dfw, c("Phylum", "Decomp_Stage"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
write.csv(Trtdata, "wallace_rel_abund_nodecomp.csv")


tiff("rel_abund_wallace.TIF", width = 1200, height = 800)
wdataplot=ggplot(Trtdataw, aes(x=Decomp_Stage,y=mean))+
  geom_bar(aes(fill = Phylum),colour="black", stat="identity")+
  xlab("Decomposition Stage")+ylab("Relative Abundance (> 1%)") + 
  scale_x_discrete(labels=c('Submerged_Fresh' ='SF', 'Early_Floating'='EF',
                            'Floating_Decay'='FD', 'Advanced_Floating_Decay'='AFD', 
                            'Sunken_Remains'='SR')) +
  scale_fill_manual(values=c("red", "orange", "yellow", "green", "blue",
                               "purple", "pink", "gray", "white", "black"))
wdataplot
dev.off()

bm <- compare_means(Abundance~Decomp_Stage, data=df,group.by="Phylum",method = "kruskal.test",p.adjust.method="fdr")
write.csv(bm, "benbow_phy_kw_stats.csv")
Means=compare_means(Abundance ~ Decomp_Stage, data = df, 
                    group.by = "Phylum", p.adjust.method = "fdr")
NComparisons<-length(unique(metadata$Decomp_Stage))*length(unique(Trtdata$Phylum))
SigList<-length(unique(Trtdata$Phylum))
SigLetters2<-vector(length=NComparisons)
#vec<-unlist(lst)
sink("benbow_phy_wil_stats.csv")
for (i in levels(Means$Phylum)){
  Tax<-i
  TaxAbundance<-subset(Means,Phylum==i )
  print(TaxAbundance)
  Hyphenated<-as.character(paste0(TaxAbundance$group1,"-",TaxAbundance$group2))
  difference<-TaxAbundance$p.adj
  names(difference)<-Hyphenated
  Letters<-multcompLetters(difference)
  return(Letters)
  
}
sink()

bm <- compare_means(Abundance~Decomp_Stage, data=df,group.by="Phylum",method = "kruskal.test",p.adjust.method="fdr")
write.csv(bm, "wallace_phy_kw_stats.csv")
Means=compare_means(Abundance ~ Decomp_Stage, data = df, 
                    group.by = "Phylum", p.adjust.method = "fdr")
NComparisons<-length(unique(metadata$Decomp_Stage))*length(unique(Trtdata$Phylum))
SigList<-length(unique(Trtdata$Phylum))
SigLetters2<-vector(length=NComparisons)
#vec<-unlist(lst)
sink("wallace_phy_wil_stats.csv")
for (i in levels(Means$Phylum)){
  Tax<-i
  TaxAbundance<-subset(Means,Phylum==i )
  print(TaxAbundance)
  Hyphenated<-as.character(paste0(TaxAbundance$group1,"-",TaxAbundance$group2))
  difference<-TaxAbundance$p.adj
  names(difference)<-Hyphenated
  Letters<-multcompLetters(difference)
  print(Letters)
  
}
sink()

GPrFamily=tax_glom(physeq_benbow, "Family")
FamilyLevel = transform_sample_counts(GPrFamily, function(x) x / sum(x))
FamilyLevel = filter_taxa(FamilyLevel, function(x) mean(x) > 0.001, TRUE) 

df <- psmelt(FamilyLevel)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum", "Family"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
d <- Trtdata

GPrFamily=tax_glom(physeq_wallace, "Family")
FamilyLevel = transform_sample_counts(GPrFamily, function(x) x / sum(x))
FamilyLevel = filter_taxa(FamilyLevel, function(x) mean(x) > 0.001, TRUE) 

df <- psmelt(FamilyLevel)
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum", "Family"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
d <- Trtdata





relative_abundance_setup <- function(physeq) {
  GPrPhylum=tax_glom(physeq, "Phylum")
  PhylumLevel_rec = transform_sample_counts(GPrPhylum, function(x) x / sum(x))
  PhylumLevel_rec = filter_taxa(PhylumLevel_rec, function(x) mean(x) > 0.001, TRUE) 
  df <- psmelt(PhylumLevel_rec) 
  df$Abundance=df$Abundance*100
  Trtdata <- ddply(df, c("Phylum", 'Study', 'Decomp_Stage'), summarise,
                   N    = length(Abundance),
                   mean = mean(Abundance),
                   sd   = sd(Abundance),
                   se   = sd / sqrt(N)
  )
}


GPrPhylum=tax_glom(physeq, "Phylum")
PhylumLevel = transform_sample_counts(GPrPhylum, function(x) x / sum(x))
PhylumLevel = filter_taxa(PhylumLevel, function(x) mean(x) > 0.001, TRUE) 
df <- psmelt(PhylumLevel)  
df_graph <- relative_abundance_setup(physeq)
write.csv(df_graph, "phylum_relative_abund_forensicpig.csv")
df_graph$Decomp_Stage <- factor(df_graph$Decomp_Stage, levels = c('Submerged_Fresh', 'Early_Floating',
                                                                  'Floating_Decay', 'Advanced_Floating_Decay', 
                                                                  'Sunken_Remains'))

bm <- compare_means(Abundance~Study, data=df, group.by="Phylum", method = "kruskal.test",p.adjust.method="bonferroni")
write.csv(bm, "study_phy_kw_stats.csv")
Means=compare_means(Abundance~Study, data = df, 
                    group.by = "Phylum", p.adjust.method = "bonferroni")
NComparisons<-length(unique(metadata$Study))*length(unique(df_graph$Phylum))
SigList<-length(unique(df_graph$Phylum))
SigLetters2<-vector(length=NComparisons)
#vec<-unlist(lst)
sink("study_phy_wil_stats.csv")
for (i in levels(Means$Phylum)){
  Tax<-i
  TaxAbundance<-subset(Means,Phylum==i )
  print(TaxAbundance)
  Hyphenated<-as.character(paste0(TaxAbundance$group1,"-",TaxAbundance$group2))
  difference<-TaxAbundance$p.adj
  names(difference)<-Hyphenated
  Letters<-multcompLetters(difference)
  return(Letters)
}
sink()


df_graph_benbow <- subset(df_graph, Study == "Benbow")
df_graph_benbow
df_graph_wallace <- subset(df_graph, Study == "Wallace")

theme_set(theme_classic(base_size = 14))
a=ggplot(df_graph_benbow, aes(x=Decomp_Stage, y=mean, fill=Decomp_Stage))+
  geom_bar(stat="identity")+
  ylab("Relative Abundance Phylum Level (> 1%)")+facet_wrap(~Phylum, scales="free")+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1.0)
a 
tiff("rel_abund_benbow.TIF", width = 1200, height = 600)
a + scale_fill_manual(values = c("#FF0000","#006699", "#FF6600", "#FFCC00")) + 
  theme(axis.text.x = element_blank()) + 
  xlab("")
dev.off()

b=ggplot(df_graph_wallace, aes(x=Decomp_Stage, y=mean, fill=Decomp_Stage))+
  geom_bar(stat="identity")+
  ylab("Relative Abundance Phylum Level (> 1%)")+facet_wrap(~Phylum, scales="free")+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1.0)
b
tiff("rel_abund_wallace.TIF", width = 1200, height = 600)
b + scale_fill_manual(values = c("#FF0000", "#006699", "#FF6600", "#FFCC00", "#00A08A")) + 
  theme(axis.text.x = element_blank()) + 
  xlab("")
dev.off()

relative_abundance_setup <- function(physeq) {
  GPrPhylum=tax_glom(physeq, "Phylum")
  PhylumLevel_rec = transform_sample_counts(GPrPhylum, function(x) x / sum(x))
  PhylumLevel_rec = filter_taxa(PhylumLevel_rec, function(x) mean(x) > 0.0001, TRUE) 
  df <- psmelt(PhylumLevel_rec) 
  df$Abundance=df$Abundance*100
  Trtdata <- ddply(df, c("Phylum", 'Study'), summarise,
                   N    = length(Abundance),
                   mean = mean(Abundance),
                   sd   = sd(Abundance),
                   se   = sd / sqrt(N)
  )
}

df_graph <- relative_abundance_setup(physeq)
df_graph$Day <- factor(df_graph$Day, levels= c('1', '4', '7', '8', '10', '12', 
                                               '13', '16', '19', '21'))
df_graph_benbow <- subset(df_graph, Study == "Benbow")
df_graph_wallace <- subset(df_graph, Study == "Wallace")

c=ggplot(df_graph_benbow, aes(x=Day, y=mean, fill=Day))+
  geom_bar(stat="identity")+
  ylab("Relative Abundance Phylum Level (> 1%)")+facet_wrap(~Phylum, scales="free")+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1.0)
c 
tiff("rel_abund_benbow_day.TIF", width = 1200, height = 600)
c +  theme(axis.text.x = element_blank()) + 
  scale_fill_manual(values = c("#66023C", "#CB2314", "#273046", "#354823", "#E1BD6D", "#A9A9A9")) +
  xlab("")
dev.off()

d=ggplot(df_graph_wallace, aes(x=Day, y=mean, fill=Day))+
  geom_bar(stat="identity")+
  ylab("Relative Abundance Phylum Level (> 1%)")+facet_wrap(~Phylum, scales="free")+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1.0)
d
tiff("rel_abund_wallace_day.TIF", width = 1200, height = 600)
d + scale_fill_manual(values = c("#66023C", "#CB2314", "#273046", 
                                 "#CC5500","#354823", "#E1BD6D", "#A9A9A9")) +
  theme(axis.text.x = element_blank()) + 
  xlab("")
dev.off()

