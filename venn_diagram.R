
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
