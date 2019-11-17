library(ggplot2)
library(tidyverse)

df <- data.frame(MOD = c('Natural', 'Accident', 'Homicide', 'Suicide'), 
                 cases = c(57,71,37,23))

df

theme_set(theme_classic(base_size = 22))
pie <- ggplot(df, aes(x="", y=cases, fill=MOD))+
  geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0) +
  geom_text(aes(label = paste0(cases)), position = position_stack(vjust = 0.5), size = 6) +
  xlab('') + ylab('Number of Cases') + 
  scale_fill_manual(values=c("#F26419", "#55DDE0", "#2F4858", "#F6AE2D")) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

pie
