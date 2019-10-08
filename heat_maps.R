#heat map
#use reduced phylum and family levels from abundance charts

(PhylumLevel)
PhylumLevelAv = transform_sample_counts(PhylumLevel, function(x) x / sum(x))
(PhylumLevelAv)
rank_names(PhylumLevelAv)

(FamilyLevel)
FamilyLevelAv = transform_sample_counts(FamilyLevel, function(x) x / sum(x))
plot_heatmap(FamilyLevelAv, taxa.label = "family")

plot_heatmap(PhylumLevelAv, taxa.label = "Phylum", sample.label="SampleID", sample.order="Sample_Area", low="red", high="blue", na.value="black")

plot_heatmap(FamilyLevelAv, taxa.label = "Family", sample.label="SampleID", sample.order="Sample_Area", low="red", high="blue", na.value="black")

plot_heatmap(physeq,  "NMDS", "bray", sample.label="SampleID", sample.order="Sample_Area", low="red", high="blue", na.value="black", title="Bray-Curtis")

heatmap.2(as.matrix(GPdist))
