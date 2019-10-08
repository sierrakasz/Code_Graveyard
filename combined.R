setwd("C:/Users/sierr/Documents/Thesis Project/")

library(ape)
library(DAtest)
library(dplyr)
library(dunn.test)
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
library(randomForest)
library(tidyr)
library(tidyverse)
library(vegan)

#import files

#mg-rast
otufull=read.csv("MG-Rast biom files\\otu_mgrast_thirty.csv",header=TRUE)

metadata=(read.csv("MG-Rast biom files\\HPMMMeta_r.csv",header=TRUE))

taxa=as.matrix(read.csv("MG-Rast biom files\\taxa_thirty.csv"))

#Formating OTU table
OTU=otu_table(otufull, taxa_are_rows=TRUE)
row.names(OTU)

TAX=tax_table(taxa)
taxa_names(TAX)

sampdat=sample_data(metadata)
row.names(sampdat)
sample_names(sampdat)=metadata$SampleID

taxa_names(TAX)=row.names(OTU)

physeq=phyloseq(OTU,TAX, sampdat)
sample_variables(physeq)
colnames(tax_table(physeq)) <- c("Domain",
                                 "Phylum",
                                 "Class",
                                 "Order",
                                 "Family",
                                 "Genus"
)

#qiime
biom <- import_biom("qiime biom files/qiime.sixty.biom")

metadata=(read.csv("Metadata/HPMMMeta_r.csv",header=TRUE))

tree=read_tree("/tree.nwk", errorIfNULL=TRUE)

sampdat=sample_data(metadata)
row.names(sampdat)
sample_names(sampdat)=metadata$SampleID


physeq=merge_phyloseq(biom, sampdat)
rank_names(physeq)
colnames(tax_table(physeq)) <- c("Domain",
                                 "Phylum",
                                 "Class",
                                 "Order",
                                 "Family",
                                 "Genus",
                                 "Species")
#mothur
biom <- import_biom("/mothur biom files/mothur.sixty.biom")

metadata=(read.csv("Metadata/HPMMMeta_r.csv",header=TRUE))

sampdat=sample_data(metadata)
row.names(sampdat)
sample_names(sampdat)=metadata$SampleID


physeq=merge_phyloseq(biom, sampdat)
rank_names(physeq)
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", 
                                 "Order", "Family",  "Genus")
#abundance charts
theme_set(theme_bw(base_size = 9)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

physeq = transform_sample_counts(physeq, function(x) x / sum(x))

GPrPhylum=tax_glom(physeq, "Phylum")
GPrPhylum
PhylumLevel = filter_taxa(GPrPhylum, function(x) mean(x) > 0.00001, TRUE) 
PhylumLevel

GPrFamily=tax_glom(physeq,"Family")
GPrFamily
FamilyLevel = filter_taxa(GPrFamily, function(x) mean(x) > 0.00001, TRUE)
FamilyLevel
f <- plot_bar(FamilyLevel, x="SampleID", fill="Family",  title="Relative Abundance Family Level")
f + facet_wrap(Sample_Area~Broad_PMI)

#remove unclassified
tax_table(PhylumLevel)

pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

badTaxa <- c("sp273", "sp272", "sp251", "sp157", "sp209", "sp295")
PhylumLevel <- pop_taxa(PhylumLevel, badTaxa)
PhylumLevel

q <-plot_bar(PhylumLevel, x="SampleID", fill="Phylum", title="Relative Abundance Phylum Level")
q + facet_wrap(~Sample_Area)

df <- psmelt(PhylumLevel) 
df$Abundance=df$Abundance*100
length(df)
Trtdata <- ddply(df, c("Phylum", "Sample_Area"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)

View(Trtdata)
write.table(Trtdata,"phylum.txt",sep="\t",row.names=FALSE)

a=ggplot(Trtdata, aes(x=Sample_Area, y=mean,fill=Sample_Area))+geom_bar(stat="identity")+ylab("Relative Abundance Phylum Level (> 1%)")+facet_wrap(~Phylum)+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1.0)
a + theme(axis.text.x = element_text(angle = 45, hjust = 1))
a + ggtitle("Relative Abundance of Microbiome Phlya") + xlab("Sample Area") + labs(fill = "Sample Area")

FamilyLevel
df <- psmelt(FamilyLevel) 
df$Abundance=df$Abundance*100
length(df)
Trtdata <- ddply(df, c("Family", "Sample_Area"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)

View(Trtdata)                
write.table(Trtdata,"family.txt",sep="\t",row.names=FALSE)

theme_set(theme_bw(base_size = 8)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

b=ggplot(Trtdata, aes(x=Sample_Area, y=mean,fill=Sample_Area))+geom_bar(stat="identity")+ylab("Relative Abundance Family Level (> 1%)")+facet_wrap(~Family, nrow=8)+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=1.0)
b + theme(axis.text.x = element_text(angle = 45, hjust = 1))#+scale_fill_manual(values=cbPalette) + guides(fill=FALSE)+ scale_x_discrete(labels=c("Original", "PS/Bran", "PS Only"))

#ANCOM
DA.anc(physeq, predictor = "Sample_Area")

#tree
plot_tree(PhylumLevel, color = "Sample_Area", 
          label.tips = "Phylum", nodelabf = nodeplotblank, 
          plot.margin = 0.5, ladderize = TRUE,
          title = "")


#Alpha diversity
alph <- prune_species(speciesSums(physeq) > 0, physeq)
#Chao1 and ACE are species richness, Simpson and shannon are diversity
plot_richness(alph, x = "Broad_PMI", color="Sample_Area", 
              measures=c("Chao1", "ACE", "Shannon", "InvSimpson"),
              scales = "fixed")

alphst = merge_samples(alph, "Sample_Area")
sample_data(alphst)$Broad_PMI <- factor(sample_names(alphst))
sample_data(alphst)$Sample_Area <- factor(sample_names(alphst))
a = plot_richness(alphst, x="Sample_Area", color="Sample_Area", 
                  measures=c("Chao1", "ACE", "Shannon", "InvSimpson"),
                  scales = "fixed")
a + geom_point(size=5, alpha=0.7)

erich <- estimate_richness(physeq, measures = c("Chao1", "ACE", "Shannon", "InvSimpson"))
View(erich)
write.table(erich,"rich.txt",sep="\t",row.names=FALSE)

#test metric for normality, use T test or kruskal results
shapiro.test(erich[,6])
#not normal: 2, 5, 6
ttest <- t(sapply(erich, function(x) unlist(t.test(x~sample_data(physeq)$Sample_Area)[c("estimate","p.value","statistic","conf.int")])))
ttest

ktest <- t(sapply(erich, function(x) unlist(kruskal.test(x~sample_data(physeq)$Sample_Area)[c("estimate","p.value","statistic","conf.int")])))
ktest

#beta diversity
#change method and distance
#method= DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"
#distance= bray, wunifrac, uunifrac, jaccard 

ord=ordinate(physeq, method="NMDS", distance="bray")
ordplot=plot_ordination(physeq, ord,"SampleID", color="Sample_Area",
                        title ="")
ordplot+ stat_ellipse(type= "norm",geom = "polygon", alpha = 1/4, aes(fill = Sample_Area))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot_heatmap(physeq, "PCoA", "bray", sample.label="SampleID", 
             sample.order="Sample_Area", 
             low="red", high="blue", na.value="black", 
             title="")

#change distance
GPdist=phyloseq::distance(physeq, "bray")
sampledf <- data.frame(sample_data(physeq))

# Adonis test
adonis(GPdist ~ Sample_Area, data = sampledf)
adonis(GPdist ~ Broad_PMI, data = sampledf)
#Homogeniy of beta dispersion
beta <- betadisper(GPdist, sampledf$Sample_Area)
beta <- betadisper(GPdist, sampledf$Broad_PMI)
permutest(beta)

#core
physeq.rel <- microbiome::transform(physeq, "compositional")
physeq.core <- core(physeq.rel, detection = 0, prevalence = .5)
core.taxa <- taxa(physeq.core)

det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)
plot_core(physeq.rel, prevalences = prevalences, detections = det, plot.type = "lineplot") + xlab("Relative Abundance (%)")

#format otu tables like they are in mg-rast 
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
  annotate("text", x = df.labels$x, y = df.labels$y, label = df.labels$counts, size = 10)

#random forest
ForestData=physeq
predictors=t(otu_table(ForestData))
dim(predictors)
response <- as.factor(sample_data(ForestData)$Sample_Number)
rf.data <- data.frame(response, predictors)
MozzieForest <- randomForest(response~., data = rf.data, ntree = 1000)
print(MozzieForest)#returns overall Random Forest results

imp <- importance(MozzieForest)#all the steps that are imp or imp. are building a dataframe that contains info about the taxa used by the Random Forest test to classify treatment 
imp <- data.frame(predictors = rownames(imp), imp)
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
imp.20 <- imp.sort[1:20, ]

ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "purple") +
  coord_flip() +
  ggtitle("Important OTUs for Classifying Samples by Person")
imp.20$MeanDecreaseGini
otunames <- imp.20$predictors
r <- rownames(tax_table(ForestData)) %in% otunames
kable(tax_table(ForestData)[r, ])#returns a list of the most important predictors for Random Forest Classification
