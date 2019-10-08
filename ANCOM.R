rm(list=ls())
getwd()
setwd("C:/Users/sierr/Documents/Forensic_Pig_Combined/")

#packages
library(ggforce)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(PMCMR)
library(phyloseq)
library(plyr)
library(tidyverse)
library(exactRankTests)
library(nlme)

#make a table in excel for rare data
tax_group <- read.csv("tax_group_fp_combo.csv")
tax_group <- tax_group %>% group_by(Kingdom, Phylum, Class, Order, Family, Genus) %>% summarise_all(funs(sum))
tax_group <- tax_group %>%  ungroup()

regroup_physeq_object <-function(table) {
  tax <- table %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
  tax <- as.matrix(tax)
  otu <- table %>% select(contains("P"))
  otu <- otu[,-1]
  OTU=otu_table(otu, taxa_are_rows=TRUE)
  TAX=tax_table(tax)
  taxa_names(TAX)=row.names(OTU)
  physeq_all=phyloseq(OTU,TAX)
  return(physeq_all)
}

physeq <- regroup_physeq_object(tax_group)

metadata=(read.csv("AquaticPigMetadataCombined.csv",header=TRUE))
sampdat=sample_data(metadata)
sample_names(sampdat)=metadata$SampleID
sampdat$Day <- factor(sampdat$Day, levels = c('1', '4', '7', '8', '10', '12', 
                                              '13', '16', '19', '21'))
sampdat$Decomp_Stage <- factor(sampdat$Decomp_Stage, levels = c('Submerged_Fresh', 'Early_Floating',
                                                                'Floating_Decay', 'Advanced_Floating_Decay', 
                                                                'Sunken_Remains'))
physeq=merge_phyloseq(physeq, sampdat)

physeq_benbow <- subset_samples(physeq, Study == "Benbow")
physeq_wallace <- subset_samples(physeq, Study == "Wallace")


pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[allTaxa %in% badTaxa]
  return(prune_taxa(allTaxa, physeq))
}

badTaxa <- c("sp1000", "sp1004", "sp1039", "sp1041", "sp1083", "sp1089", "sp1124",
             "sp1133", "sp19", "sp2045", "sp34", "sp934", "sp948", "sp953",
             "sp960", "sp967", "sp981", "sp971", "sp977", "sp997")
physeq_benbow_ancom <- pop_taxa(physeq_benbow, badTaxa)

badTaxa2 <- c("sp1225", "sp1669", "sp1749", "sp1761", "sp1944", "sp1951", "sp1972",
              "sp1986", "sp1991", "sp2057", "sp2060", "sp2086", "sp2119",
              "sp2332", "sp320", "sp373", "sp431", "sp458", "sp501", "sp561")
physeq_wallace_ancom <- pop_taxa(physeq_wallace, badTaxa2)

#code from https://sites.google.com/site/siddharthamandal1985/research
ancom.W = function(otu_data,var_data,
                   adjusted,repeated,
                   main.var,adj.formula,
                   repeat.var,long,rand.formula,
                   multcorr,sig){
  
  n_otu=dim(otu_data)[2]-1
  
  otu_ids=colnames(otu_data)[-1]
  
  if(repeated==F){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID",all.y=T),row.names=NULL)
    #data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var)],by="Sample.ID",all.y=T),row.names=NULL)
  }else if(repeated==T){
    data_comp=data.frame(merge(otu_data,var_data,by="Sample.ID"),row.names=NULL)
    # data_comp=data.frame(merge(otu_data,var_data[,c("Sample.ID",main.var,repeat.var)],by="Sample.ID"),row.names=NULL)
  }
  
  base.formula = paste0("lr ~ ",main.var)
  if(repeated==T){
    repeat.formula = paste0(base.formula," | ", repeat.var)
  }
  if(adjusted==T){
    adjusted.formula = paste0(base.formula," + ", adj.formula)
  }
  
  if( adjusted == F & repeated == F ){
    fformula  <- formula(base.formula)
  } else if( adjusted == F & repeated == T & long == T ){
    fformula  <- formula(base.formula)   
  }else if( adjusted == F & repeated == T & long == F ){
    fformula  <- formula(repeat.formula)   
  }else if( adjusted == T & repeated == F  ){
    fformula  <- formula(adjusted.formula)   
  }else if( adjusted == T & repeated == T  ){
    fformula  <- formula(adjusted.formula)   
  }else{
    stop("Problem with data. Dataset should contain OTU abundances, groups, 
         and optionally an ID for repeated measures.")
  }
  
  
  
  if( repeated==FALSE & adjusted == FALSE){
    if( length(unique(data_comp[,which(colnames(data_comp)==main.var)]))==2 ){
      tfun <- exactRankTests::wilcox.exact
    } else{
      tfun <- stats::kruskal.test
    }
  }else if( repeated==FALSE & adjusted == TRUE){
    tfun <- stats::aov
  }else if( repeated== TRUE & adjusted == FALSE & long == FALSE){
    tfun <- stats::friedman.test
  }else if( repeated== TRUE & adjusted == FALSE & long == TRUE){
    tfun <- nlme::lme
  }else if( repeated== TRUE & adjusted == TRUE){
    tfun <- nlme::lme
  }
  
  logratio.mat <- matrix(NA, nrow=n_otu, ncol=n_otu)
  for(ii in 1:(n_otu-1)){
    for(jj in (ii+1):n_otu){
      data.pair <- data_comp[,which(colnames(data_comp)%in%otu_ids[c(ii,jj)])]
      lr <- log((1+as.numeric(data.pair[,1]))/(1+as.numeric(data.pair[,2])))
      
      lr_dat <- data.frame( lr=lr, data_comp,row.names=NULL )
      
      if(adjusted==FALSE&repeated==FALSE){  ## Wilcox, Kruskal Wallis
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==FALSE&repeated==TRUE&long==FALSE){ ## Friedman's 
        logratio.mat[ii,jj] <- tfun( formula=fformula, data = lr_dat)$p.value
      }else if(adjusted==TRUE&repeated==FALSE){ ## ANOVA
        model=tfun(formula=fformula, data = lr_dat,na.action=na.omit)   
        picker=which(gsub(" ","",row.names(summary(model)[[1]]))==main.var)  
        logratio.mat[ii,jj] <- summary(model)[[1]][["Pr(>F)"]][picker]
      }else if(repeated==TRUE&long==TRUE){ ## GEE
        model=tfun(fixed=fformula,data = lr_dat,
                   random = formula(rand.formula),
                   correlation=corAR1(),
                   na.action=na.omit)   
        picker=which(gsub(" ","",row.names(anova(model)))==main.var)
        logratio.mat[ii,jj] <- anova(model)[["p-value"]][picker]
      }
      
    }
  } 
  
  ind <- lower.tri(logratio.mat)
  logratio.mat[ind] <- t(logratio.mat)[ind]
  
  
  logratio.mat[which(is.finite(logratio.mat)==FALSE)] <- 1
  
  mc.pval <- t(apply(logratio.mat,1,function(x){
    s <- p.adjust(x, method = "BH")
    return(s)
  }))
  
  a <- logratio.mat[upper.tri(logratio.mat,diag=FALSE)==TRUE]
  
  b <- matrix(0,ncol=n_otu,nrow=n_otu)
  b[upper.tri(b)==T] <- p.adjust(a, method = "BH")
  diag(b)  <- NA
  ind.1    <- lower.tri(b)
  b[ind.1] <- t(b)[ind.1]
  
  #########################################
  ### Code to extract surrogate p-value
  surr.pval <- apply(mc.pval,1,function(x){
    s0=quantile(x[which(as.numeric(as.character(x))<sig)],0.95)
    # s0=max(x[which(as.numeric(as.character(x))<alpha)])
    return(s0)
  })
  #########################################
  ### Conservative
  if(multcorr==1){
    W <- apply(b,1,function(x){
      subp <- length(which(x<sig))
    })
    ### Moderate
  } else if(multcorr==2){
    W <- apply(mc.pval,1,function(x){
      subp <- length(which(x<sig))
    })
    ### No correction
  } else if(multcorr==3){
    W <- apply(logratio.mat,1,function(x){
      subp <- length(which(x<sig))
    })
  }
  
  return(W)
  }



ANCOM.main = function(OTUdat,Vardat,
                      adjusted,repeated,
                      main.var,adj.formula,
                      repeat.var,longitudinal,
                      random.formula,
                      multcorr,sig,
                      prev.cut){
  
  p.zeroes=apply(OTUdat[,-1],2,function(x){
    s=length(which(x==0))/length(x)
  })
  
  zeroes.dist=data.frame(colnames(OTUdat)[-1],p.zeroes,row.names=NULL)
  colnames(zeroes.dist)=c("Taxon","Proportion_zero")
  
  zero.plot = ggplot(zeroes.dist, aes(x=Proportion_zero)) + 
    geom_histogram(binwidth=0.1,colour="black",fill="white") + 
    xlab("Proportion of zeroes") + ylab("Number of taxa") +
    theme_bw()
  
  #print(zero.plot)
  
  OTUdat.thinned=OTUdat
  OTUdat.thinned=OTUdat.thinned[,c(1,1+which(p.zeroes<prev.cut))]
  
  otu.names=colnames(OTUdat.thinned)[-1]
  
  W.detected   <- ancom.W(OTUdat.thinned,Vardat,
                          adjusted,repeated,
                          main.var,adj.formula,
                          repeat.var,longitudinal,random.formula,
                          multcorr,sig)
  
  W_stat       <- W.detected
  
  
  ### Bubble plot
  
  W_frame = data.frame(otu.names,W_stat,row.names=NULL)
  W_frame = W_frame[order(-W_frame$W_stat),]
  
  W_frame$detected_0.9=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.8=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.7=rep(FALSE,dim(W_frame)[1])
  W_frame$detected_0.6=rep(FALSE,dim(W_frame)[1])
  
  W_frame$detected_0.9[which(W_frame$W_stat>0.9*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.8[which(W_frame$W_stat>0.8*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.7[which(W_frame$W_stat>0.7*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  W_frame$detected_0.6[which(W_frame$W_stat>0.6*(dim(OTUdat.thinned[,-1])[2]-1))]=TRUE
  
  final_results=list(W_frame,zero.plot)
  names(final_results)=c("W.taxa","PLot.zeroes")
  return(final_results)
}

decomp_pairs_benbow <- function(physeq) {
  phy_SF <- subset_samples(physeq, Decomp_Stage == 'Submerged_Fresh')
  phy_EF <- subset_samples(physeq, Decomp_Stage == 'Early_Floating')
  phy_FD <- subset_samples(physeq, Decomp_Stage == 'Floating_Decay')
  phy_AFD <- subset_samples(physeq, Decomp_Stage == 'Advanced_Floating_Decay')
  return(list(phy_EF.FD <- merge_phyloseq(phy_EF, phy_FD), 
              phy_FD.AFD <- merge_phyloseq(phy_AFD, phy_FD),
              phy_EF.AFD <- merge_phyloseq(phy_EF, phy_AFD),
              phy_SF.EF <- merge_phyloseq(phy_SF, phy_EF),
              phy_SF.FD <- merge_phyloseq(phy_SF, phy_FD),
              phy_SF.AFD <- merge_phyloseq(phy_SF, phy_AFD)))
}

decomp_list_benbow <- decomp_pairs_benbow(physeq_benbow_ancom)

decomp_pairs_wallace <- function(physeq) {
  phy_SF <- subset_samples(physeq, Decomp_Stage == 'Submerged_Fresh')
  phy_EF <- subset_samples(physeq, Decomp_Stage == 'Early_Floating')
  phy_FD <- subset_samples(physeq, Decomp_Stage == 'Floating_Decay')
  phy_AFD <- subset_samples(physeq, Decomp_Stage == 'Advanced_Floating_Decay')
  phy_SR <- subset_samples(physeq, Decomp_Stage == 'Sunken_Remains')
  return(list(phy_EF.FD <- merge_phyloseq(phy_EF, phy_FD), 
              phy_FD.AFD <- merge_phyloseq(phy_AFD, phy_FD),
              phy_EF.AFD <- merge_phyloseq(phy_EF, phy_AFD),
              phy_SF.EF <- merge_phyloseq(phy_SF, phy_EF),
              phy_SF.FD <- merge_phyloseq(phy_SF, phy_FD),
              phy_SF.AFD <- merge_phyloseq(phy_SF, phy_AFD),
              phy_SF.SR <- merge_phyloseq(phy_SF, phy_SR),
              phy_EF.SR <- merge_phyloseq(phy_EF, phy_SR),
              phy_AFD.SR <- merge_phyloseq(phy_AFD, phy_SR),
              phy_FD.SR <- merge_phyloseq(phy_FD, phy_SR)))
}

decomp_list_wallace <- decomp_pairs_wallace(physeq_wallace_ancom)

day_pairs_benbow <- function(physeq) {
  physeq_1 <- subset_samples(physeq, Day == '1')
  physeq_4 <- subset_samples(physeq, Day == '4')
  physeq_8 <- subset_samples(physeq, Day == '8')
  physeq_12 <- subset_samples(physeq, Day == '12')
  physeq_16 <- subset_samples(physeq, Day == '16')
  physeq_21 <- subset_samples(physeq, Day == '21')
  return(list(physeq_1.4 <- merge_phyloseq(physeq_1, physeq_4),
              physeq_1.8 <- merge_phyloseq(physeq_1, physeq_8),
              physeq_1.12 <- merge_phyloseq(physeq_1, physeq_12),
              physeq_1.16 <- merge_phyloseq(physeq_1, physeq_16),
              physeq_1.21 <- merge_phyloseq(physeq_1, physeq_21),
              physeq_4.8 <- merge_phyloseq(physeq_8, physeq_4),
              physeq_4.12 <- merge_phyloseq(physeq_12, physeq_4),
              physeq_4.16 <- merge_phyloseq(physeq_16, physeq_4),
              physeq_4.21 <- merge_phyloseq(physeq_21, physeq_4),
              physeq_8.12 <- merge_phyloseq(physeq_8, physeq_12),
              physeq_8.16 <- merge_phyloseq(physeq_8, physeq_16),
              physeq_8.21 <- merge_phyloseq(physeq_8, physeq_21),
              physeq_12.16 <- merge_phyloseq(physeq_16, physeq_12),
              physeq_12.21 <- merge_phyloseq(physeq_21, physeq_12),
              physeq_16.21 <- merge_phyloseq(physeq_16, physeq_21)))
}

day_list_benbow <- day_pairs_benbow(physeq_benbow)

day_pairs_wallace <- function(physeq) {
  physeq_1 <- subset_samples(physeq, Day == '1')
  physeq_4 <- subset_samples(physeq, Day == '4')
  physeq_7 <- subset_samples(physeq, Day == '7')
  physeq_10 <- subset_samples(physeq, Day == '10')
  physeq_13 <- subset_samples(physeq, Day == '13')
  physeq_16 <- subset_samples(physeq, Day == '16')
  physeq_19 <- subset_samples(physeq, Day == '19')
  return(list(physeq_1.4 <- merge_phyloseq(physeq_1, physeq_4),
              physeq_1.7 <- merge_phyloseq(physeq_1, physeq_7),
              physeq_1.10 <- merge_phyloseq(physeq_1, physeq_10),
              physeq_1.13 <- merge_phyloseq(physeq_1, physeq_13),
              physeq_1.16 <- merge_phyloseq(physeq_1, physeq_16),
              physeq_1.19 <- merge_phyloseq(physeq_1, physeq_19),
              physeq_4.7 <- merge_phyloseq(physeq_7, physeq_4),
              physeq_4.10 <- merge_phyloseq(physeq_4, physeq_10),
              physeq_4.13 <- merge_phyloseq(physeq_13, physeq_4),
              physeq_4.16 <- merge_phyloseq(physeq_16, physeq_4),
              physeq_4.19 <- merge_phyloseq(physeq_19, physeq_4),
              physeq_7.10 <- merge_phyloseq(physeq_7, physeq_10),
              physeq_7.13 <- merge_phyloseq(physeq_7, physeq_13),
              physeq_7.16 <- merge_phyloseq(physeq_7, physeq_16),
              physeq_7.19 <- merge_phyloseq(physeq_7, physeq_19),
              physeq_10.13 <- merge_phyloseq(physeq_10, physeq_13),
              physeq_10.16 <- merge_phyloseq(physeq_10, physeq_16),
              physeq_10.19 <- merge_phyloseq(physeq_10, physeq_19),
              physeq_13.16 <- merge_phyloseq(physeq_16, physeq_13),
              physeq_13.19 <- merge_phyloseq(physeq_19, physeq_13),
              physeq_16.19 <- merge_phyloseq(physeq_16, physeq_19)))
}

day_list_wallace <- day_pairs_wallace(physeq_wallace)


otu_ancom_make <- function(physeq) {
  otu_ancom <- data.frame(otu_table(physeq))
  otu_ancom <- data.frame(t(otu_ancom))
  Sample.ID <- rownames(otu_ancom)
  rownames(otu_ancom) <- NULL
  otu_ancom <- cbind(Sample.ID, otu_ancom)
  return(otu_ancom)
}

metadata_ancom <- metadata
colnames(metadata_ancom)[1] <- 'Sample.ID'

ancom_results_decomp_benbow <- list()
ancom_results_day_benbow <- list()
ancom_results_decomp_wallace <- list()
ancom_results_day_wallace <- list()

for(i in 1:length(decomp_list_benbow)) {
  otu_ancom <- otu_ancom_make(decomp_list_benbow[[i]])
  comparison_test <- ANCOM.main(OTUdat = otu_ancom,
                                Vardat = metadata_ancom,
                                adjusted = FALSE,
                                repeated = F,
                                main.var = "Decomp_Stage",
                                adj.formula = NULL,
                                repeat.var=NULL,
                                longitudinal=FALSE,
                                random.formula=NULL,
                                multcorr=2,
                                sig=0.05,
                                prev.cut=0.90)
  w_values <- data.frame(comparison_test$W.taxa)
  tax <- data.frame(tax_table(physeq_benbow))
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
  tax$otu.names <- rownames(tax)
  ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
  ancom_sign_taxa <- ancom_sign_taxa[,-1]
  ancom_results_decomp_benbow[[i]] <- ancom_sign_taxa
}

for(i in 1:length(day_list_benbow)) {
  otu_ancom <- otu_ancom_make(day_list_benbow[[i]])
  comparison_test <- ANCOM.main(OTUdat = otu_ancom,
                                Vardat = metadata_ancom,
                                adjusted = FALSE,
                                repeated = F,
                                main.var = "Day",
                                adj.formula = NULL,
                                repeat.var=NULL,
                                longitudinal=FALSE,
                                random.formula=NULL,
                                multcorr=2,
                                sig=0.05,
                                prev.cut=0.90)
  w_values <- data.frame(comparison_test$W.taxa)
  tax <- data.frame(tax_table(physeq_benbow))
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
  tax$otu.names <- rownames(tax)
  ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
  ancom_sign_taxa <- ancom_sign_taxa[,-1]
  ancom_results_day_benbow[[i]] <- ancom_sign_taxa
}

for(i in 1:length(decomp_list_wallace)) {
  otu_ancom <- otu_ancom_make(decomp_list_wallace[[1]])
  comparison_test <- ANCOM.main(OTUdat = otu_ancom,
                                Vardat = metadata_ancom,
                                adjusted = FALSE,
                                repeated = F,
                                main.var = "Decomp_Stage",
                                adj.formula = NULL,
                                repeat.var=NULL,
                                longitudinal=FALSE,
                                random.formula=NULL,
                                multcorr=2,
                                sig=0.05,
                                prev.cut=0.90)
  w_values <- data.frame(comparison_test$W.taxa)
  tax <- data.frame(tax_table(physeq_wallace))
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
  tax$otu.names <- rownames(tax)
  ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
  ancom_sign_taxa <- ancom_sign_taxa[,-1]
  ancom_results_decomp_wallace[[i]] <- ancom_sign_taxa
}

for(i in 1:length(day_list_wallace)) {
  otu_ancom <- otu_ancom_make(day_list_wallace[[i]])
  comparison_test <- ANCOM.main(OTUdat = otu_ancom,
                                Vardat = metadata_ancom,
                                adjusted = FALSE,
                                repeated = F,
                                main.var = "Day",
                                adj.formula = NULL,
                                repeat.var=NULL,
                                longitudinal=FALSE,
                                random.formula=NULL,
                                multcorr=2,
                                sig=0.05,
                                prev.cut=0.90)
  w_values <- data.frame(comparison_test$W.taxa)
  tax <- data.frame(tax_table(physeq_wallace))
  tax <- tax %>% select(Kingdom, Phylum, Class, Order, Family, Genus)
  tax$otu.names <- rownames(tax)
  ancom_sign_taxa <- merge(w_values, tax, by='otu.names')
  ancom_sign_taxa <- ancom_sign_taxa[,-1]
  ancom_results_day_wallace[[i]] <- ancom_sign_taxa
}

sink("decomp_benbow_ANCOM.csv")
for(i in 1:length(ancom_results_decomp_benbow)) {
  print(ancom_results_decomp_benbow[[i]])
}
sink()

sink("day_benbow_ANCOM.csv")
for(i in 1:length(ancom_results_day_benbow)) {
  print(ancom_results_day_benbow[[i]]%>% filter(detected_0.9 == 'TRUE'))
}
sink()

sink("decomp_wallace_ANCOM.csv")
for(i in 1:length(ancom_results_decomp_wallace)) {
  print(ancom_results_decomp_wallace[[i]])
}
sink()


sink("day_wallace_ANCOM.csv")
for(i in 1:length(ancom_results_day_wallace)) {
  print(ancom_results_day_wallace[[i]]%>% filter(detected_0.9 == 'TRUE'))
}
sink()
