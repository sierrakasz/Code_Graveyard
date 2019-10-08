alpha <- read.csv("C:/Users/sierr/Documents/Thesis project/R Script/alphdiv.csv")
View(alpha)

shapiro.test(alpha$Shannon)
bartlett.test(Shannon ~ pipeline, data = alpha)
res.aov <- aov(Shannon ~ pipeline, data = alpha)
summary(res.aov)
TukeyHSD(res.aov)

shapiro.test(alpha$InvSimpson)
fligner.test(InvSimpson ~ pipeline, data = alpha)
kruskal.test(alpha$InvSimpson, alpha$pipeline)
dunn.test(alpha$InvSimpson, alpha$pipeline)
pairwise.wilcox.test(alpha$InvSimpson, alpha$pipeline, p.adjust.method = "BH")


phy <- read.csv("R Script/relabundance_phylum.csv")
View(phy)
shapiro.test(phy$Means)
fligner.test(Means ~ Pipeline, data = phy)
kruskal.test(phy$Means, phy$Pipeline)
dunn.test(phy$Means, phy$Pipeline)

fam <- read.csv("R Script/relabundance_family.csv")
View(fam)
shapiro.test(fam$Means)
fligner.test(Means ~ Pipeline, data = fam)
kruskal.test(fam$Means, fam$Pipeline)
dunn.test(fam$Means, fam$Pipeline)
