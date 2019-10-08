setwd("C:/Users/sierr/Documents/ThesisProject_Pipelines/")

#packages
#library(ape)
library(DESeq2)
library(dplyr)
library(FSA)
library(ggforce)
library(ggplot2)
library(ggpubr)
library(gplots)
library(ggthemes)
library(knitr)
library(lme4)
library(microbiome)
library(phyloseq)
library(plyr)
library(plotly)
library(PMCMR)
library(randomForest)
library(tidyr)
library(tidyverse)
library(UpSetR)
library(vegan)

#import files
metadata=(read.csv("Metadata/HPMMMeta_r_merge.csv",header=TRUE))

sampdat=sample_data(metadata)
row.names(sampdat)
sample_names(sampdat)=metadata$SampleID

physeq=merge_phyloseq(physeq_all, sampdat)

physeq <- rarefy_even_depth(physeq, 
                            rngseed = 711, replace = FALSE, trimOTUs = FALSE)
physeq

#abundance charts
#by sample
theme_set(theme_bw(base_size = 8)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

physeq = transform_sample_counts(physeq, function(x) x / sum(x))

physeq_30 <- subset_samples(physeq, Subsample == "30")
physeq_60 <- subset_samples(physeq, Subsample == "60")
physeq_120 <- subset_samples(physeq, Subsample == "120")
physeq_188 <- subset_samples(physeq, Subsample == "188")
merge36 <- merge_phyloseq(physeq_30, physeq_60)
merge312 <- merge_phyloseq(physeq_30, physeq_120)
merge318 <- merge_phyloseq(physeq_30, physeq_188)
merge612 <- merge_phyloseq(physeq_120, physeq_60)
merge618 <- merge_phyloseq(physeq_188, physeq_60)
merge1218 <- merge_phyloseq(physeq_188, physeq_120)

physeq_Q <- subset_samples(physeq, Pipeline == "QIIME")
physeq_M <- subset_samples(physeq, Pipeline == "mothur")
physeq_G <- subset_samples(physeq, Pipeline == "MG-RAST")
physeq_QM <- merge_phyloseq(physeq_Q, physeq_M)
physeq_QG <- merge_phyloseq(physeq_Q, physeq_G)
physeq_GM <- merge_phyloseq(physeq_M, physeq_G)

physeq_rec <- subset_samples(physeq_188, Sample_Area == "Rectum")
physeq_mou <- subset_samples(physeq_188, Sample_Area == "Mouth")



GPrPhylum=tax_glom(physeq, "Phylum")
GPrPhylum
PhylumLevel = filter_taxa(GPrPhylum, function(x) mean(x) > 0.00001, TRUE) 
PhylumLevel

#remove unclassified
tax_table(PhylumLevel)

pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

badTaxa <- c("sp351", "sp325")
PhylumLevel <- pop_taxa(PhylumLevel, badTaxa)
PhylumLevel

PhylumLevel = transform_sample_counts(PhylumLevel, function(x) x / sum(x))

q <-plot_bar(PhylumLevel, x="SampleID", fill="Phylum", title="Relative Abundance Phylum Level\nMG-Rast Subsample 30")
q + facet_wrap(~MoD)
q + facet_grid(rows = vars(Sample_Area), cols=vars(MoD))

#realtive abundance by variable
df <- psmelt(PhylumLevel) 
df$Abundance=df$Abundance*100
length(df)
Trtdata <- ddply(df, c("Phylum", "Pipeline"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)

View(Trtdata)
write.table(Trtdata,"phylum.txt",sep="\t",row.names=FALSE)

a=ggplot(Trtdata, aes(x=MoD, y=mean,fill=MoD))+geom_bar(stat="identity")+ylab("Relative Abundance Phylum Level (> 1%)")+facet_wrap(~Phylum)+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1.0)
a + theme(axis.text.x = element_text(angle = 45, hjust = 1))
a + ggtitle("Relative Abundance of Microbiome Phlya\nMG-RAST Subsample 30") + xlab("Manner of Death") + labs(fill = "Manner of Death")

#reset physeq object
diagdds = phyloseq_to_deseq2(physeq_QM, ~ Pipeline)
#change for which one is on top or bottom
diagdds$Pipeline <- relevel(diagdds$Pipeline , ref = "mothur")
gm_mean = function(x, na.rm=TRUE){  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")

res = results(diagdds, cooksCutoff = FALSE)
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq)[rownames(sigtab), ], "matrix"))
head(sigtab)

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))

log_plot <- ggplot(sigtab, aes(x=Phylum, y=log2FoldChange, color=Phylum)) +
  geom_boxplot( aes(x=Phylum, y=log2FoldChange),outlier.shape=1) +
  stat_boxplot( aes(x=Phylum, y=log2FoldChange), 
                                     geom='errorbar', linetype=1, width=0.5) +
  geom_point(size=2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("mothur vs. QIIME Differential Abundance\nPhylum Level")
log_plot

log_plot + geom_errorbar(aes(ymin=log2FoldChange-lfcSE,ymax=log2FoldChange+lfcSE),color="black",width=1.0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("mothur vs. QIIME Differential Abundance\nPhylum Level")


#tree
plot_tree(PhylumLevel, color = "Sample_Area", 
          label.tips = "Phylum", nodelabf = nodeplotblank, 
          plot.margin = 0.5, ladderize = TRUE,
          title = "")

#Alpha diversity
alph <- prune_species(speciesSums(physeq) > 0, physeq)
#Chao1 and ACE are species richness, Simpson and shannon are diversity
plot_richness(alph, x = "Sample_Area", color="Sample_Area", 
              measures=c("Chao1", "ACE", "Shannon", "InvSimpson"),
              scales = "free",
              title = "Alpha Diversity Metrics\nMG-RAST Subsample 30")

#alphst = merge_samples(alph, "Sample_Area")
#sample_data(alphst)$Broad_PMI <- factor(sample_names(alphst))
#sample_data(alphst)$Sample_Area <- factor(sample_names(alphst))
#a = plot_richness(alphst, x="Sample_Area", color="Sample_Area", 
#                  measures=c("Chao1", "ACE", "Shannon", "InvSimpson"),
#                  scales = "fixed")
#a + geom_point(size=5, alpha=0.7)

erich <- estimate_richness(physeq_188, measures = c("Observed", "Chao1", "ACE", "Shannon", "InvSimpson"))
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

#test metric for normality, use T test or kruskal results
shapiro.test(erich[,6])
#not normal: 2, 5, 6
ttest <- t(sapply(erich, function(x) unlist(t.test(x~sample_data(physeq)$Sample_Area)[c("estimate","p.value","statistic","conf.int")])))
ttest

ktest <- t(sapply(erich, function(x) unlist(kruskal.test(x~sample_data(physeq)$Sample_Area)[c("estimate","p.value","statistic","conf.int")])))
ktest

met <- metadata[,c(1,3,4)]
dunn <- merge(met, erich)
dunnTest(Chao1~Pipeline, data=dunn, method = 'bonferroni')

#beta diversity
#change method and distance
#method= DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"
#distance= bray, wunifrac, uunifrac, jaccard 

ord = ordinate(physeq_188, method="PCoA", distance="bray")
ordplot=plot_ordination(physeq_188, ord,"SampleID", color="Pipeline",
                        title ="PCoA of Bray Curtis Index\nMG-RAST Subsample 30") +
  stat_ellipse(geom = "polygon", alpha = 1/8, aes(color = Pipeline))
ordplot 

ordplot+ stat_ellipse(type= "norm",geom = "polygon", alpha = 1/4, aes(fill = Sample_Area))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#plot_heatmap(physeq, "PCoA", "bray", sample.label="SampleID", 
#             sample.order="Sample_Area", 
#             low="red", high="blue", na.value="black", 
#             title="PCoA of Bray Curtis Index\nMG-RAST Subsample 30")

#change distance
GPdist=phyloseq::distance(physeq_188, "bray")
sampledf <- data.frame(sample_data(physeq_188))

# Beta diversity
adonis(GPdist ~ Pipeline, data = sampledf)

#Homogeniy of beta dispersion
beta <- betadisper(GPdist, sampledf$Pipeline)
permutest(beta)

#core
#heatmap
physeq.rel <- microbiome::transform(physeq, "compositional")
physeq.core <- core(physeq.rel, detection = 0, prevalence = .5)
core.taxa <- taxa(physeq.core)

det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)
plot_core(physeq.rel, prevalences = prevalences, detections = det, plot.type = "lineplot") + xlab("Relative Abundance (%)") +
  ggtitle("Core Microbiome Size\nMG-RAST Subsample 30")

detections <- 10^seq(log10(1), log10(max(abundances(physeq))/10), length = 10)

library(RColorBrewer)
p <- plot_core(physeq, plot.type = "heatmap", 
               prevalences = seq(0.1, 1, .1),
               detections = detections,
               colours = rev(brewer.pal(5, "Spectral")),
               min.prevalence = .2, horizontal = TRUE)
p + theme_gray(base_size = 14
) + ggtitle("Core Microbiome\nMG-RAST Subsample 30")

#core venn diagram
#mg-rast 
otu.data <- otufull
View(otufull)
metadata <- read_csv("HPMMMeta_r.csv")

#qiime and mothur
OTU1 = as(otu_table(physeq), "matrix")
OTUdf = as.data.frame(OTU1)
View(OTUdf)
rownames(OTUdf) <- c()
otu.data <- OTUdf
View(otu.data)

#core by sample area
df <- otu.data %>% mutate(otu.num = 1:nrow(otu.data)) %>% 
  select(otu.num, everything()) %>% 
  gather(contains("WCME"), key = "SampleID", value = otu.val) %>% 
  left_join(metadata, by = "SampleID") %>%
  select(SampleID, otu.num, otu.val, Sample_Area) %>% 
  filter(otu.val > 0) %>%
  group_by(Sample_Area, otu.num) %>%
  tally()

mouth.df <- df %>% filter(Sample_Area == "Mouth")
rectum.df <- df %>% filter(Sample_Area == "Rectum")

intersect.val <- nrow(semi_join(mouth.df, rectum.df, by = "otu.num"))
mouth.val <- nrow(anti_join(mouth.df, rectum.df, by = "otu.num"))
rectum.val <- nrow(anti_join(rectum.df, mouth.df, by = "otu.num"))

df.venn <- tibble(x = c(0,0),
                  y = c(1,-0.5),
                  labels = c('Mouth', 'Rectum'))

df.labels <- tibble(x = c(0,0,0), y = c(1.5,0.3,-1), counts = c(mouth.val, intersect.val, rectum.val))

ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .3, size = 1, colour = 'grey') +
  coord_fixed() +
  theme_void() +
  labs(fill = NULL) +
  annotate("text", x = df.labels$x, y = df.labels$y, label = df.labels$counts, size = 10) +
  ggtitle("Core Microbiome Between Sample Area\nMG-RAST Subsample 30")

#core by MoD
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
ForestData=physeq_G
predictors=t(otu_table(ForestData))
dim(predictors)
response <- as.factor(sample_data(ForestData)$MoD)
rf.data <- data.frame(response, predictors)
Forest <- randomForest(response~., data = rf.data, ntree = 1000)
print(Forest)#returns overall Random Forest results

imp <- importance(Forest)#all the steps that are imp or imp. are building a dataframe that contains info about the taxa used by the Random Forest test to classify treatment 
imp <- data.frame(predictors = rownames(imp), imp)
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
imp.20 <- imp.sort[1:20, ]
tax <- data.frame(tax_table(physeq))
head(tax)
tax <- tax[,1:6]
tax$predictors <- rownames(tax)
imp.20 <- merge(imp.20, tax)
write.table(imp.20, sep=",", "random_forest_predictors.csv", append = TRUE)

ggplot(imp.20, aes(x=predictors, y=MeanDecreaseGini, color=MeanDecreaseGini)) +
  geom_point(size = 3) +
  geom_segment(aes(x=predictors, 
                   xend=predictors, 
                   y=0, 
                   yend=MeanDecreaseGini)) +
  labs(title="Important Taxa Predictors of Manner of Death for MG-RAST")

