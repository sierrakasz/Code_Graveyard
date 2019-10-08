library("plyr")
library("dplyr")
library("ggplot2")

otufull=read.table("C:\\Users\\sierr\\Documents\\Thesis project\\MG-Rast biom files\\otu_mgrast_thirty.txt",header=TRUE)#For other column add Other to end of name
head(otufull)
metadata=read.table("C:\\Users\\sierr\\Documents\\Thesis project\\HPMMMeta_r.txt",header=TRUE)
head(metadata)
taxa=as.matrix(read.csv("C:\\Users\\sierr\\Documents\\Thesis project\\taxa_thirty.csv"))#For other column add Other to end of name
head(taxa)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#000000","#CC79A7")
theme_set(theme_bw(base_size = 20)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))#sets the plotting theme

OTU=otu_table(otufull, taxa_are_rows=TRUE)
row.names(OTU)#head(OTU) 
(OTU)
TAX=tax_table(taxa)
taxa_names(TAX)
TAX
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
taxa_names(TAX)=row.names(OTU)
taxa_names(TAX)
physeq=phyloseq(OTU,TAX,sampdat)
sample_data(physeq)
sample_variables(physeq)
physeq

plot_bar(physeq, "Sample_Area", "Abundance", "phylum")
physeq = transform_sample_counts(physeq, function(x) x / sum(x))

GPrPhylum=tax_glom(physeq, "phylum")
GPrPhylum
PhylumLevel = filter_taxa(GPrPhylum, function(x) mean(x) > 1e-5, TRUE) #filter out any taxa lower tha 0.1%
PhylumLevel
GPrFamily=tax_glom(GPr,"Family")
FamilyLevel = filter_taxa(GPrFamily, function(x) mean(x) > 2.3e-2, TRUE) #filter out any taxa lower tha 1%
GPrGenus=tax_glom(GPr,"Genus")
GenusLevel = filter_taxa(GPrGenus, function(x) mean(x) > 2e-2, TRUE) #filter out any taxa lower tha 1%

#SubsetGut=subset_samples(FamilyLevel, Sample_Type!="Gut")
GenusLevel
df <- psmelt(PhylumLevel)#Change to phylum family or genus as needed
df$Abundance=df$Abundance*100
length(df)
Trtdata <- ddply(df, c("phylum", "Sample_Area"), summarise,
                 N    = length(Abundance),
                 mean = mean(Abundance),
                 sd   = sd(Abundance),
                 se   = sd / sqrt(N)
)
Trtdata
theme_set(theme_bw(base_size = 18)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

p=ggplot(Trtdata, aes(x=Sample_Area, y=mean,fill=Sample_Area))+geom_bar(stat="identity")+ylab("Relative Abundance (> 1%)")+facet_wrap(~phylum)+geom_errorbar(aes(ymin=mean-se,ymax=mean+se),color="black",width=0.5)
p + theme(axis.text.x = element_text(angle = 45, hjust = 1))#+scale_fill_manual(values=cbPalette) + guides(fill=FALSE)+ scale_x_discrete(labels=c("Original", "PS/Bran", "PS Only"))

ord=ordinate(physeq, "unifrac")
ordplot=plot_ordination(physeq, ord,"samples", color="Diet")#+geom_point(size=4)+scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)
ordplot+ stat_ellipse(type= "norm",geom = "polygon", alpha = 1/4, aes(fill = Diet))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())#+ theme(legend.justification=c(1,0), legend.position=c(1,0))

##Ordination by sample type
ordplot=plot_ordination(physeq, ord,"samples", color="Sample_Type")#+geom_point(size=4)+scale_colour_manual(values=cbPalette)+scale_fill_manual(values=cbPalette)
ordplot+ stat_ellipse(type= "norm",geom = "polygon", alpha = 1/4, aes(fill = Sample_Type))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())#+ theme(legend.justification=c(1,0), legend.position=c(1,0))




#PERMANOVAs

##Diet

GPdist=phyloseq::distance(physeq, "bray")
#MONMDS= ordinate(physeq, "NMDS",GPdist)


GPdist=phyloseq::distance(physeq, "wunifrac")

adonis(GPdist ~ Diet, as(sample_data(physeq), "data.frame"))
#By Diet and Sample type
#adonis(GPdist ~ Diet*Sample_Type, as(sample_data(physeq), "data.frame"))
#By Sample type
#adonis(GPdist ~ Sample_Type, as(sample_data(physeq), "data.frame"))


#Random Forest
##By Diet

ForestData=physeq#Change this one so you dont have to rewrite all variables
predictors=t(otu_table(ForestData))
dim(predictors)
response <- as.factor(sample_data(ForestData)$Sample_Area)
rf.data <- data.frame(response, predictors)
MozzieForest <- randomForest(response~., data = rf.data, ntree = 1000)
print(MozzieForest)#returns overall Random Forest results
imp <- importance(MozzieForest)#all the steps that are imp or imp. are building a dataframe that contains info about the taxa used by the Random Forest testto classify treatment 
imp <- data.frame(predictors = rownames(imp), imp)
imp.sort <- arrange(imp, desc(MeanDecreaseGini))
imp.sort$predictors <- factor(imp.sort$predictors, levels = imp.sort$predictors)
imp.20 <- imp.sort[1:20, ]
ggplot(imp.20, aes(x = predictors, y = MeanDecreaseGini)) +
  geom_bar(stat = "identity", fill = "indianred") +
  coord_flip() +
  ggtitle("Most important OTUs for classifying  samples\n by Diet")#\n in a string tells it to start a new line
imp.20$MeanDecreaseGini
otunames <- imp.20$predictors
r <- rownames(tax_table(ForestData)) %in% otunames
kable(tax_table(ForestData)[r, ])#returns a list of the most important predictors for Random Forest Classification