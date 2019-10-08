#change method and distance
#method= DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA"
#distance= bray, wunifrac, uunifrac, jaccard 

ord=ordinate(physeq, method="NMDS", distance="bray")
ordplot=plot_ordination(physeq, ord,"SampleID", color="Sample_Area")
ordplot+ stat_ellipse(type= "norm",geom = "polygon", alpha = 1/4, aes(fill = Sample_Area))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#change distance
GPdist=phyloseq::distance(physeq, "bray")

sampledf <- data.frame(sample_data(physeq))

# Adonis test
adonis(GPdist ~ Sample_Area, data = sampledf)

#Homogeniy of beta dispersion

beta <- betadisper(GPdist, sampledf$Sample_Area)

permutest(beta)

