theme_set(theme_bw(base_size = 12)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))#sets the plotting theme

#Alpha diversity
alph <- prune_species(speciesSums(physeq) > 0, physeq)
alph
#Chao1 and ACE are species richness, Simpson and shannon are diversity
plot_richness(alph, x = "Broad_PMI", color="Sample_Area", measures=c("Chao1", "ACE", "Shannon", "InvSimpson"))

alphst = merge_samples(alph, "Sample_Area")
sample_data(alphst)$Broad_PMI <- factor(sample_names(alphst))
sample_data(alphst)$Sample_Area <- factor(sample_names(alphst))
a = plot_richness(alphst, x="Sample_Area", color="Sample_Area", measures=c("Chao1", "ACE", "Shannon", "InvSimpson"))
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



#richness and eveness
nsamp= nsamples(physeq)
for (i in 1:100) {
  # Calculate richness
  rich <- as.numeric(as.matrix(estimate_richness(physeq, measures = "Observed")))
  richness[ ,i] <- rich
  # Calculate evenness
  even <- as.numeric(as.matrix(estimate_richness(physeq, measures = "InvSimpson")))
  evenness[ ,i] <- even
}

SampleID <- row.names(richness)
mean <- apply(richness, 1, mean)
sd <- apply(richness, 1, sd)
measure <- rep("Richness", nsamp)
rich_stats <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of evenness estimates
SampleID <- row.names(evenness)
mean <- apply(evenness, 1, mean)
sd <- apply(evenness, 1, sd)
measure <- rep("Inverse Simpson", nsamp)
even_stats <- data.frame(SampleID, mean, sd, measure)


alpha <- rbind(rich_stats, even_stats)
View(alpha)
View(s)
s <- data.frame(sample_data(physeq))
alphadiv <- merge(alpha, s, by = 'SampleID') 
#alphadiv <- order_by(alphadiv)

ggplot(alphadiv, aes(x = BMI, y = mean, color = Sample_Area, group = Sample_Area, shape = Sample_Area)) +
  geom_point(size = 2) + 
  geom_line(size = 0.8) +
  facet_wrap(~measure, ncol = 1, scales = "free") +
  scale_color_manual(values = c("#E96446", "#302F3D", "#87CEFA")) +
  scale_x_discrete(
    breaks = c("7/8", "8/4", "9/2", "10/6"),
    labels = c("Jul", "Aug", "Sep", "Oct"), 
    drop = FALSE
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
