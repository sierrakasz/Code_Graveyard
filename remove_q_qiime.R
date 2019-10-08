#fixing taxa names from QIIME
a <- data.frame(tax_table(biom))

b <- a %>% 
  mutate(Rank1 = str_replace(Rank1, "D_0__", ""))
head(b)

c <- b %>% 
  mutate(Rank2 = str_replace(Rank2, "D_1__", ""))
head(c)

d <- c %>% 
  mutate(Rank3 = str_replace(Rank3, "D_2__", ""))

e <- d %>% 
  mutate(Rank4 = str_replace(Rank4, "D_3__", ""))

f <- e %>% 
  mutate(Rank5 = str_replace(Rank5, "D_4__", ""))

g <- f %>% 
  mutate(Rank6 = str_replace(Rank6, "D_5__", ""))

h <- g[,-7]
head(h)

h <- as.matrix(h)

i <- data.frame(otu_table(biom))

OTUq=otu_table(i, taxa_are_rows=TRUE)
row.names(OTUq)

TAXq=tax_table(h)
taxa_names(TAXq)

taxa_names(TAXq)=row.names(OTUq)

sample_names(OTUq)

physeq_qiime = merge_phyloseq(OTUq, TAXq)
