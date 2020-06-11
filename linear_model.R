rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/ThesisProject_AKP/")

#packages
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(lme4)
library(phyloseq)
library(piecewiseSEM)
library(plyr)
library(PMCMR)
library(tidyverse)
library(vegan)


core_data <- read.csv('core_percentage_norare.csv', header= T)

hist(core_data$Core_percent)


core_data$FinePMI <- factor(core_data$FinePMI, levels = c("Less24", "25-48", "49-72", "Great73"))

ggplot(data = core_data, aes(x=FinePMI, y = Core_percent, color = Sample_Area)) +
  geom_boxplot()

#add random effects of family/genus 

m0 <- lmer(Core_percent ~ 1 + (1|Pack_ID), data = core_data, REML = F)
m1 <- lmer(Core_percent ~ FinePMI + (1|Pack_ID), data = core_data, REML = F)
m1.1 <- lmer(Core_percent ~ Sample_Area + (1|Pack_ID), data = core_data, REML = F)
m2 <- lmer(Core_percent ~ FinePMI + Sample_Area + (1|Pack_ID), data = core_data, REML = F)
m3 <- lmer(Core_percent ~ FinePMI + Sample_Area + CoD_Simple + (1|Pack_ID), data = core_data, REML = F)
m4 <- lmer(Core_percent ~ FinePMI + Sample_Area + Violent + (1|Pack_ID), data = core_data, REML = F)
m5 <- lmer(Core_percent ~ FinePMI + Sample_Area + Violent + CoD_Simple +(1|Pack_ID), data = core_data, REML = F)
m6 <- lmer(Core_percent ~ FinePMI + Sample_Area + Violent * CoD_Simple + (1|Pack_ID), data = core_data, REML = F)

anova(m0, m1, m1.1, m2, m3, m4, m5, m6, test = 'Chisq')

rsquared(m0)
rsquared(m1)
rsquared(m1.1)
rsquared(m2)
rsquared(m3)
rsquared(m4)
rsquared(m5)
rsquared(m6)

summary(m1.1)
