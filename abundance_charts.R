theme_set(theme_bw(base_size = 9)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))

physeq = transform_sample_counts(physeq, function(x) x / sum(x))
p <-plot_bar(physeq, "Sample_Area", "Abundance") + scale_fill_manual(values=cbPalette)
plot(p)

#bar plots with taxa and samples all together 
#remember to go to import and remove unclassified
#go through phylum first and then family
rank_names(physeq)
GPrPhylum=tax_glom(physeq, "Phylum")
GPrPhylum
PhylumLevel = filter_taxa(GPrPhylum, function(x) mean(x) > 0.00001, TRUE) 
PhylumLevel
#remove unclassified
tax_table(PhylumLevel)
q <-plot_bar(PhylumLevel, x="SampleID", fill="Phylum", title="Relative Abundance Phylum Level")
q + facet_wrap(~Sample_Area)

GPrFamily=tax_glom(physeq,"Family")
GPrFamily
FamilyLevel = filter_taxa(GPrFamily, function(x) mean(x) > 0.00001, TRUE)
FamilyLevel
f <- plot_bar(FamilyLevel, x="SampleID", fill="Family",  title="Relative Abundance Family Level")
f + facet_wrap(Sample_Area~Broad_PMI)


#bar plots with taxa spearated out
#think about comparing abundances
PhylumLevel
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
a + theme(axis.text.x = element_text(angle = 45, hjust = 1))#+scale_fill_manual(values=cbPalette) + guides(fill=FALSE)+ scale_x_discrete(labels=c("Original", "PS/Bran", "PS Only"))
plot(a)

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



