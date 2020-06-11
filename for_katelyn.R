#for katelyn 

#change this depending on the MOD/COD
preparing_data_for_core <- function(physeq) {
  otus <- data.frame(otu_table(physeq))
  totus <- data.frame(t(otus))
  totus$SampleID <- rownames(totus)
  #change MoD/ CoD_Simple2 when applicable
  met <- metadata[,c('SampleID', 'Suicide')]
  mtotus <- merge(totus, met)
  mtotus <- mtotus[,-1]
  total <- as.vector(colSums(Filter(is.numeric, mtotus)))
  new_df <- mtotus %>% group_by(Suicide) %>% summarise_all(funs(sum))
  new_df <- data.frame(t(new_df))
  colnames(new_df) <- as.character(unlist(new_df[1,]))
  new_df = new_df[-1, ]
  new_df$OTU <- rownames(new_df)
  rownames(new_df) <- NULL
  Upset <- cbind(new_df, total)
  return(Upset)
}

all_otus_mod_rec <- preparing_data_for_core(physeq_by_bodysites[['Rectum']])
all_otus_mod_mou <- preparing_data_for_core(physeq_by_bodysites[['Mouth']])

#clean up the data
#change variable name based on bodysite
all_otus_mod_mou <- all_otus_mod_mou %>% filter(total != 0)
all_otus_mod_mou$Suicide <- as.numeric(as.character(all_otus_mod_mou$Suicide))
all_otus_mod_mou$Non <- as.numeric(as.character(all_otus_mod_mou$Non))

#otus withins MODs only
#change variable names based on bodysite
non_otus_mod_mou <- all_otus_mod_mou %>% filter(Non != 0) %>% 
  filter(Suicide == 0) 
sui_otus_mod_mou <- all_otus_mod_mou %>% filter(Non == 0) %>% 
  filter(Suicide != 0) 

#subset physeq object using the important taxa
pop_taxa = function(physeq, impTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[(allTaxa %in% impTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

#pull out just taxa names, not boruta decision
impTaxa <- sui_otus_mod_mou$OTU
#additional X in names from some reason, remove it
impTaxa <- gsub('X', '', impTaxa)
#new phyloseq obj.
#change based on bodysite
core_mod_mou_bd <- pop_taxa(physeq_by_bodysites[['Mouth']], impTaxa)
core_mod_mou_bd

#clean up the data
#change variable name based on bodysite
all_otus_mod_rec <- all_otus_mod_rec %>% filter(total != 0)
all_otus_mod_rec$Suicide <- as.numeric(as.character(all_otus_mod_rec$Suicide))
all_otus_mod_rec$Non <- as.numeric(as.character(all_otus_mod_rec$Non))

#otus withins MODs only
#change variable names based on bodysite
non_otus_mod_rec <- all_otus_mod_rec %>% filter(Non != 0) %>% 
  filter(Suicide == 0) 
sui_otus_mod_rec <- all_otus_mod_rec %>% filter(Non == 0) %>% 
  filter(Suicide != 0) 

#subset physeq object using the important taxa
pop_taxa = function(physeq, impTaxa){
  allTaxa = taxa_names(physeq)
  allTaxa <- allTaxa[(allTaxa %in% impTaxa)]
  return(prune_taxa(allTaxa, physeq))
}

#pull out just taxa names, not boruta decision
impTaxa <- sui_otus_mod_rec$OTU
#additional X in names from some reason, remove it
impTaxa <- gsub('X', '', impTaxa)
#new phyloseq obj.
#change based on bodysite
core_mod_rec_bd <- pop_taxa(physeq_by_bodysites[['Rectum']], impTaxa)
core_mod_rec_bd

case_study_mou <- read.csv('sui_core_taxa_mou.csv')
case_study_rec <- read.csv('sui_core_taxa_rec.csv')

tax_mou <- data.frame(tax_table(core_mod_mou_bd))

tax_mou <- read.csv('sui_Full_taxa_mou.csv')

tax_mou_dataf <- rbind(tax_mou, case_study_mou)

tax_mou_dd <- tax_mou_dataf[duplicated(tax_mou_dataf),]

write.csv(tax_mou_dd, 'sui_robust_taxa_mou.csv')

tax_rec <- data.frame(tax_table(core_mod_mou_bd))

write.csv(tax_rec, 'sui_Full_taxa_rec.csv')

tax_rec <- read.csv('sui_Full_taxa_rec.csv')

tax_rec_dataf <- rbind(tax_rec, case_study_rec)

tax_rec_df <- tax_rec_dataf[duplicated(tax_rec_dataf),]

write.csv(tax_rec_df, 'sui_robust_taxa_rec.csv')
