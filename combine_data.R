#import qiime files
biom <- import_biom("C:/Users/sierr/Documents/Thesis project/qiime biom files/qiime.thirty.biom")
(biom)
tax_table(biom)
otu_table(biom)
sample_names(biom) <- paste("Q_", sample_names(biom), sep="")
sample_names(biom)

physeq_qiime = merge_phyloseq(biom)
physeq_qiime

colnames(tax_table(physeq_qiime)) <- c("Kingdom",
                                 "Phylum",
                                 "Class",
                                 "Order",
                                 "Family",
                                 "Genus")
#import mothur
biomm <- import_biom("C:/Users/sierr/Documents/Thesis project/mothur biom files/mothur.thirty.biom")
(biomm)
tax_table(biomm)
otu_table(biomm)
sample_names(biomm) <- paste("M_", sample_names(biomm), sep="")
sample_names(biomm)

physeq_mothur = merge_phyloseq(biomm)
physeq_mothur
colnames(tax_table(physeq_mothur)) <- c("Kingdom", "Phylum", "Class", 
                                 "Order", "Family",  "Genus")

#import MG-RAST
otufull=read.csv("MG-Rast biom files\\otu_mgrast_thirty.csv",header=TRUE)
taxa=as.matrix(read.csv("MG-Rast biom files\\taxa_thirty.csv"))

OTU=otu_table(otufull, taxa_are_rows=TRUE)
row.names(OTU)

TAX=tax_table(taxa)
taxa_names(TAX)

taxa_names(TAX)=row.names(OTU)

sample_names(OTU) <- paste("G_", sample_names(OTU), sep="")
sample_names(OTU)

physeq_MG=phyloseq(OTU,TAX)
colnames(tax_table(physeq_MG)) <- c("Kingdom",
                                 "Phylum",
                                 "Class",
                                 "Order",
                                 "Family",
                                 "Genus"
)

#merge physeq
merge = merge_phyloseq(physeq_qiime, physeq_mothur, physeq_MG)
merge

#add metadata and re-merge
metadata=(read.csv("C:/Users/sierr/Documents/Thesis project/Metadata/HPMMMeta_r_merge.csv",header=TRUE))
View(metadata)

sampdat=sample_data(metadata)
row.names(sampdat)
(sampdat)
sample_names(sampdat)=metadata$SampleID


merge=merge_phyloseq(merge, sampdat)

#abundance
theme_set(theme_bw(base_size = 9)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

merge = transform_sample_counts(merge, function(x) x / sum(x))
merge = filter_taxa(merge, function(x) mean(x) > 0.00001, TRUE) 

GPrPhylum=tax_glom(merge, "Phylum")
GPrPhylum
PhylumLevel = filter_taxa(GPrPhylum, function(x) mean(x) > 0.00001, TRUE) 
PhylumLevel
tax_table(PhylumLevel)
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
q + facet_wrap(~Pipeline)
q + facet_grid(rows = vars(Sample_Area), cols=vars(MoD))

#ordinate
ord=ordinate(merge, method="PCoA", distance="jaccard")
ordplot=plot_ordination(merge, ord,"SampleID", color="Pipeline",
                        title ="PCoA of Bray Curtis Index\nMG-RAST Subsample 30")# +
  #stat_ellipse(geom = "polygon", alpha = 1/8, aes(color = Sample_Area))
ordplot 

GPdist=phyloseq::distance(merge, "bray")
sampledf <- data.frame(sample_data(merge))

# Beta diversity
adonis(GPdist ~ Pipeline*Sample_Area, data = sampledf)

#Homogeniy of beta dispersion
beta <- betadisper(GPdist, sampledf$Pipeline)
permutest(beta)


x <- data.frame(tax_table(merge))
write.csv(y, "otu_all.csv")
y <- data.frame(otu_table(merge))
total <- merge(x,y)
View(total)