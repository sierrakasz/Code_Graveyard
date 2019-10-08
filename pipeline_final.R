rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/ThesisProject_Pipelines/")

#packages
library(FSA)
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(microbiome)
library(phyloseq)
library(plyr)
library(PMCMR)
library(randomForest)
library(tidyverse)
library(UpSetR)
library(vegan)

#import mothur
import_mothur_files <- function(shared, tax) {
  physeq_mothur <- import_mothur(mothur_shared_file = shared,
                                 mothur_constaxonomy_file = tax)
  colnames(tax_table(physeq_mothur)) <- c("Kingdom", "Phylum", "Class", 
                                             "Order", "Family",  "Genus")
  return(physeq_mothur)
}

physeq_mothur <- import_mothur_files("mothur biom files/mothur.all.shared", 
                    "mothur biom files/mothur.all.taxonomy")
sample_names(physeq_mothur) <- paste("M_188_", sample_names(physeq_mothur), sep="")

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

physeq_qiime <- import_qiime_files("qiime biom files/table-with-taxonomy.biom")
sample_names(physeq_qiime) <- paste("Q_188_", sample_names(physeq_qiime), sep="")

#import MG-RAST
import_MGRAST <- function(file) {
  file_name <- read.csv(file)
  otu <- select(file_name, contains('WCME'))
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
physeq_MG <- import_MGRAST("MG-Rast biom files\\HPMMSAll.csv")
sample_names(physeq_MG) <- paste("G_188_", sample_names(physeq_MG), sep="")

#merging phyloseq objects together
merge = merge_phyloseq(physeq_qiime,physeq_mothur,physeq_MG)

#rarefaction
physeq_rare <- rarefy_even_depth(merge, 
                            rngseed = 711, sample.size = 1000)

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
tax_group <- tidying_data(physeq_rare)
write.csv(tax_group, "tax_group.csv")
tax_group <- read.csv("tax_group.csv")
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

#import metadata and combine
metadata=(read.csv("Metadata/HPMMMeta_r_merge.csv",header=TRUE))
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
physeq=merge_phyloseq(physeq_all, sampdat)

#Removing unwanted taxa
tax_table(PhylumLevel)
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}
badTaxa <- c("sp351", "sp325")
PhylumLevel <- pop_taxa(PhylumLevel, badTaxa)

#subsample data
physeq_rec <- subset_samples(physeq, Sample_Area == "Rectum")
physeq_mou <- subset_samples(physeq, Sample_Area == "Mouth")

physeq_Q_r <- subset_samples(physeq_rec, Pipeline == "QIIME")
physeq_M_r <- subset_samples(physeq_rec, Pipeline == "mothur")
physeq_G_r <- subset_samples(physeq_rec, Pipeline == "MG-RAST")

physeq_QM_r <- merge_phyloseq(physeq_Q_r, physeq_M_r)
physeq_QG_r <- merge_phyloseq(physeq_Q_r, physeq_G_r)
physeq_GM_r <- merge_phyloseq(physeq_M_r, physeq_G_r)

physeq_Q_m <- subset_samples(physeq_mou, Pipeline == "QIIME")
physeq_M_m <- subset_samples(physeq_mou, Pipeline == "mothur")
physeq_G_m <- subset_samples(physeq_mou, Pipeline == "MG-RAST")

physeq_QM_m <- merge_phyloseq(physeq_Q_m, physeq_M_m)
physeq_QG_m <- merge_phyloseq(physeq_Q_m, physeq_G_m)
physeq_GM_m <- merge_phyloseq(physeq_M_m, physeq_G_m)


#Tax glom
GPrPhylum=tax_glom(physeq_rec, "Phylum")
PhylumLevel = transform_sample_counts(GPrPhylum, function(x) x / sum(x))
PhylumLevel = filter_taxa(PhylumLevel, function(x) mean(x) > 0.0001, TRUE) 

#Relative abundance
#realtive abundance by variable
df <- psmelt(PhylumLevel) 
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Phylum", "Pipeline"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
write.csv(Trtdata, "rel_abundance_phy_rec.csv")

theme_set(theme_classic(base_size = 11))
a=ggplot(Trtdata, aes(x=Pipeline, y=mean, fill=Pipeline))+
  geom_bar(stat="identity")+
  ylab("Relative Abundance Phylum Level (> 1%)")+facet_wrap(~Phylum, scales="free")+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1.0)
a
tiff("rel_abund_rec_phy.TIF", width = 800, height = 600)
a + scale_fill_manual(values = c("#787878", "#ffb31a", "#5c5c8a")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab("")
dev.off()

#relative abundance family level
#Tax glom
GPrFamily=tax_glom(physeq_rec, "Family")
FamilyLevel = transform_sample_counts(GPrFamily, function(x) x / sum(x))
FamilyLevel = filter_taxa(FamilyLevel, function(x) mean(x) > 0.0001, TRUE) 

#Relative abundance
#realtive abundance by variable
df <- psmelt(FamilyLevel) 
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Family", "Pipeline"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
write.csv(Trtdata, "rel_abundance_rec_fam.csv")

theme_set(theme_classic(base_size = 11))
a=ggplot(Trtdata, aes(x=Pipeline, y=mean, fill=Pipeline))+
  geom_bar(stat="identity")+
  ylab("Relative Abundance Family Level (> 1%)")+facet_wrap(~Family, scales="free")+
  geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1.0)
tiff("rel_abund_rec_fam.TIF", width = 1200, height = 1000)
a + scale_fill_manual(values = c("#787878", "#ffb31a", "#5c5c8a")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab("")
dev.off()

#simulated data
df <- psmelt(physeq) 
df$Abundance=df$Abundance*100
Trtdata <- ddply(df, c("Kingdom", "Phylum", "Class", 
                       "Order", "Family",  "Genus", 
                       "Sample_Area"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
write.csv(Trtdata, "rel_abundance_insilico.csv")

#ANCOM
#open ANCOM.R and run code to make functions
otu_ancom_make <- function(physeq) {
  otu_ancom <- data.frame(otu_table(physeq))
  otu_ancom <- data.frame(t(otu_ancom))
  Sample.ID <- rownames(otu_ancom)
  rownames(otu_ancom) <- NULL
  otu_ancom <- cbind(Sample.ID, otu_ancom)
  return(otu_ancom)
}
otu_ancom <- otu_ancom_make(physeq_QG_r)

metadata_ancom <- metadata
colnames(metadata_ancom)[1] <- "Sample.ID"

comparison_test <- ANCOM.main(OTUdat = otu_ancom,
                              Vardat = metadata_ancom,
                              adjusted = FALSE,
                              repeated = F,
                              main.var = "Pipeline",
                              adj.formula = NULL,
                              repeat.var=NULL,
                              longitudinal=FALSE,
                              random.formula=NULL,
                              multcorr=2,
                              sig=0.05,
                              prev.cut=0.90)
w_values <- data.frame(comparison_test$W.taxa)
tax <- data.frame(tax_table(physeq_QG_r))
tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
tax$otu.names <- rownames(tax)
ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
ancom_sign_taxa <- ancom_sign_taxa[,-1]
ancom_sign_taxa <- ancom_sign_taxa %>% filter(detected_0.9 == TRUE)
write.table(ancom_sign_taxa, sep=",", "ancom_sign_taxa_pipeline_QG_r.csv", row.names=FALSE)

#Alpha Diversity
#Fix this
#remove indexing 
#figure out stats
#put plots together

erich <- estimate_richness(physeq, measures = c("Observed", "Chao1", "ACE", "Shannon", "InvSimpson"))
write.table(erich,"rich.txt",sep="\t",row.names=TRUE)
#fix in excel, add qiime values
View(erich)

#make data tidy 
erich2 <- erich %>%
  gather(Index, Observation, c("Observed", "Chao1", "ACE", "Shannon", "InvSimpson"), na.rm = TRUE)
erich2 <- erich2[,3:5]
View(erich2)

#add metadata
rich = merge(erich2, metadata)
View(rich)

rich <- erich[,1:4]
div <- erich[,-2:-4]

rich2 <- rich %>%
  gather(Index, Observation, c("Observed", "Chao1", "ACE"), na.rm = TRUE)
div2 <- div %>%
  gather(Index, Observation, c("Shannon", "InvSimpson"), na.rm = TRUE)

#add metadata
rich = merge(rich2, metadata)
div = merge(div2, metadata)

p <- ggplot(rich, aes(x=Index, y=Observation, fill=Index)) +
  geom_boxplot() +
  labs(title="Alpha Diversity Metrics")
p + facet_grid(Subsample~Pipeline, scales="free")

t <- ggplot(div, aes(x=Index, y=Observation, fill=Index)) +
  geom_boxplot() +
  labs(title="Alpha Diversity Metrics")
t + facet_grid(Subsample~Pipeline, scales="free")

p + facet_grid(Subsample~Pipeline, scales="fixed")
t + facet_grid(Subsample~Pipeline, scales="fixed")

#kruskal wallis test
#nemenyi  test

#kruskal.test(alpha_div~pipeline, data = alpha)
#out <- posthoc.kruskal.nemenyi.test(x=alpha_div, g=pipeline, dist="Tukey", p.adjust.method = 'bonf')
#print(otu$statistic)
#https://cran.r-project.org/web/packages/PMCMR/vignettes/PMCMR.pdf

ktest <- t(sapply(erich, function(x) unlist(kruskal.test(x~sample_data(physeq)$Sample_Area)[c("estimate","p.value","statistic","conf.int")])))
ktest
posthoc.kruskal.nemenyi.test()

met <- metadata[,c(1,3,4)]
dunn <- merge(met, erich)
dunnTest(Chao1~Pipeline, data=dunn, method = 'bonferroni')

#beta diversity
#distance= wunifrac, jaccard 
theme_set(theme_classic(base_size = 12))
ord = ordinate(physeq_mou, method="PCoA", distance="(A+B-2*J)/(A+B-J)")
ordplot=plot_ordination(physeq_mou, ord, color="Pipeline")
ordplot
tiff("total_beta_jaccard.TIF", width = 800, height = 600)
ordplot+ stat_ellipse(type= "norm",geom = "polygon", alpha = 0, size=1.1, aes(color = Sample_Area))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_point(aes(color = sample_data(physeq)$Pipeline)) + scale_color_manual(values = c("#787878", "#ffb31a", "#34495E", "#5c5c8a", "#D2B48C"))
dev.off()

tiff("mou_beta_jaccard.TIF", width = 800, height = 600)
ordplot+ geom_point(aes(color = sample_data(physeq_mou)$Pipeline)) + scale_color_manual(values = c("#787878", "#ffb31a", "#5c5c8a"))
dev.off()

#change distance
GPdist=phyloseq::distance(physeq_QG_m, "(A+B-2*J)/(A+B-J)")
sampledf <- data.frame(sample_data(physeq_QG_m))

# Beta diversity
adonis(GPdist ~ Pipeline, data = sampledf)

#Homogeniy of beta dispersion
beta <- betadisper(GPdist, sampledf$Pipeline)
permutest(beta)

#core
#UpsetR - need to fix this code
#Add core microbiome idea
#transpose otu data, Make sample ID the row names

## fix this ##
otus <- data.frame(otu_table(physeq_188))
head(otus)
totus <- data.frame(t(otus))
head(totus)
rownames(totus)
totus$SampleID <- rownames(totus)

#merge metadata, and convert sample IDs to just MoD
met <- metadata[,c(1,4)]
mtotus <- merge(totus, met)
head(mtotus)
mtotus <- mtotus[,-1]
total <- as.vector(colSums(Filter(is.numeric, mtotus)))
View(total)

#bring MoD to first column
new_df <- mtotus %>% group_by(Pipeline) %>% summarise_all(funs(sum))
head(new_df)
new_df <- data.frame(t(new_df))
colnames(new_df) <- as.character(unlist(new_df[1,]))
new_df = new_df[-1, ]
new_df$OTU <- rownames(new_df)
rownames(new_df) <- NULL

sapply(new_df, function(x) if(is.numeric(x)) replace(x,x>0,1) else x)

Upset <- cbind(new_df, total)
upset <- Upset[,c(4,1,2,3,5)]

#export to excel to change format
#Change to binary per OTU
#OTU number,(binary), Total

write.csv(upset,"pipeline_upset.csv", row.names=FALSE)
pipeCore=read.csv("pipeline_upset.csv",header=TRUE)
View(pipeCore)
pipeCore[,2:4]=sapply(pipeCore[,2:4],as.character)
sapply(pipeCore, function(x) if(is.numeric(x)) replace(x,x>0,1) else x)

upset(pipeCore, point.size = 3, line.size = 1.5, order.by = c("degree", 
                                                              "freq"), 
      empty.intersections = "on",
      mainbar.y.label = "Number of Taxa", sets.x.label = "Core Taxa", 
      queries = list(list(query = intersects, 
                          params = list("MG.RAST","mothur", "QIIME"), color= 'orange',
                          active = T), 
                     list(query = intersects,
                          params = list("MG.RAST"), color = 'red', active = T),
                     list(query = intersects, 
                          params = list("mothur"), color = "green", active = T),
                     list(query = intersects, 
                          params = list("QIIME"), color = "blue", active = T)))

#random forest
#functions
random_foresting_pipeline <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$Pipeline)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data, ntree=1000)
  return(Forest)
}
random_foresting_sample_area <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$Sample_Area)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data, ntree=1000)
  return(Forest)
}
random_foresting_MoD <- function(data_phy) {
  fd = data_phy
  predictors=t(otu_table(fd))
  dim(predictors)
  resp <- as.factor(sample_data(fd)$MoD)
  rf.data <- data.frame(resp, predictors)
  Forest <- randomForest(resp~., data=rf.data, ntree=1000)
  return(Forest)
}

forest_predictors <- function(forest) {
  imp <- importance(forest)
  imp <- data.frame(predictors = rownames(imp), imp)
  imp.sort <- arrange(imp, desc(MeanDecreaseGini))
  imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
  imp.20 <- imp.sort[1:20, ]
  tax <- data.frame(tax_table(physeq))
  head(tax)
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
  tax$predictors <- rownames(tax)
  imp.20 <- merge(imp.20, tax)
  write.table(imp.20, sep=",", "random_forest_predictors.csv", append = TRUE)
  return(imp.20)
}
                              
forest_predictors(Forest)

ggplot(imp.20, aes(x=predictors, y=MeanDecreaseGini, color=MeanDecreaseGini)) +
  geom_point(size = 3) +
  geom_segment(aes(x=predictors, 
                   xend=predictors, 
                   y=0, 
                   yend=MeanDecreaseGini)) +
  labs(title="Important Taxa Predictors of Manner of Death for MG-RAST")