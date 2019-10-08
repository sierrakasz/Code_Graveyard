
read_report <- read.csv("C:/Users/sierr/Documents/Thesis project/R Script/unclassified_reads.csv")
View(read_report)

family.df <- read_report %>% filter(Taxonomic.Level == "Family")
phylum.df <- read_report %>% filter(Taxonomic.Level == "Phylum")

family_mg <- family.df %>% filter(pipeline == "MG-RAST")
family_mot <- family.df %>% filter(pipeline == "mothur")
family <- rbind(family_mg, family_mot)

kruskal.test(phylum.df$reads, phylum.df$pipeline, data=phylum.df, p.adjust.method = "bonf")

dunnTest(reads~pipeline, data=phylum.df, method = 'bonferroni')

posthoc.kruskal.nemenyi.test(phylum.df$reads, phylum.df$pipeline, data=phylum.df, p.adjust.method = "bonf")

pairwise.wilcox.test(phylum.df$reads, phylum.df$pipeline, data=phylum.df, p.adjust.method = "bonf")


levels(read_report$pipeline)

group_by(read_report, pipeline) %>%
  summarise(
    count = n(),
    mean = mean(reads, na.rm = TRUE),
    sd = sd(reads, na.rm = TRUE)
  )

#boxplot
#change the labels
mycompare <- list(c("mothur", "QIIME"), c("mothur", "MG-RAST"), c("QIIME", "MG-RAST"))
p <- ggbarplot(read_report, x = "pipeline", y = "reads", 
          order = c("mothur", "QIIME", "MG-RAST", "Unfiltered"),
          fill = "Subsample", palette = c("#00AFBB", "#E7B800", "#FC4E07", "#00CC33"),
          ylab = "Number of Reads", xlab = "Pipeline",
          title= "Total Reads Post Filtering Between Pipelines",
          position= position_dodge()) + stat_compare_means(label = 'p.format', comparisons = mycompare)
#facet(p, facet.by = "Taxonomic.Level")
plot(p)

class(read_report)
read_reportb <- read_report[1:12,1:3]
View(read_reportb)




#Is the data normal? Signifcant = non normal
shapiro.test(read_report$reads)

#Decides hetero or homoscedastic t test
#Normal data, variance, signficant = equal variance
bartlett.test(reads ~ pipeline, data = read_report)

#Non normal data variance, signficant = equal variance
fligner.test(reads ~ pipeline, data = read_report)

#Anova
res.aov <- aov(reads ~ pipeline, data = read_report)
summary(res.aov)

#comparing means between groups 
TukeyHSD(res.aov)

#Kruskal test, non para
kruskal.test(read_report$reads, read_report$pipeline)

#Pairwise, non para
dunn.test(read_report$reads, read_report$pipeline)

#Pairwise, non para bigger sample
pairwise.wilcox.test(read_report$reads, read_report$pipeline, p.adjust.method = "BH")
