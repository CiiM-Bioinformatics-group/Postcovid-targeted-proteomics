rm(list = ls())
try(dev.off())

# Longitudinal analysis for each of the biomarkers found in step 1

library(dplyr)
library(ggplot2)
library(magrittr)
library(openxlsx)
library(ggpubr)
library(reshape2)
library(limma)

setwd("/Users/martijnzoodsma/Documents/PhD/corona/Postcovid-targeted-proteomics/")
load('data/data.RData')

PVALUE_CUTOFF <- 0.05

#### Longitudinal analysis with the biomarkers from step 1
biomarkers <- unique(c(
  read.xlsx('output/biomarkers.xlsx', sheet= 'ICU', colNames = F)$X1, 
  read.xlsx('output/biomarkers.xlsx', sheet = 'nonICU', colNames = F)$X1
))

biomarkers

# annot.radboud %<>% filter(time %in% c('W1T1', 'W1T2', 'W1T3', 'W2T1', 'W2T2', 'W2T3'))
annot.radboud %<>% filter(time %in% c('W1T1', 'W1T2', 'W1T3'))
annot.radboud <- na.omit(annot.radboud)
radboud <- radboud[which(rownames(radboud) %in% rownames(annot.radboud)), ]

means <- colMeans(radboud, na.rm = T)

df <- cbind(radboud, annot.radboud %>% dplyr::select(time, condition, sampleID)) %>% 
  reshape2::melt(., id.vars = c('sampleID', 'time', 'condition'))

# In this way we select a subset of samples that are at least somewhat consistently sampled.
samples <- df %>% filter(time == 'W1T2') %>% pull(sampleID) %>% unique()
df <- df %>% filter(sampleID %in% samples)

pvals <- data.frame('biomarker' = biomarkers, pval.time = NA, pval.condition = NA, pval.interaction = NA)

for (biomarker in biomarkers) {
  df2 <- df %>% filter(variable == biomarker)
  
  # make wide df again so we can quickly replace all the NA
  melt.sub2 <- reshape2::dcast(df2, sampleID + condition ~ time)
  
  # Look for the rows where the NA count is too high. These are patients that switched ICU - nonICU and will have 2/3 timepoints in each.
  # We dont want to do imputation for these. Instead, choose the 
  na_count <- apply(melt.sub2, 1, function(x) sum(is.na(x)))
  melt.sub3 <- melt.sub2[!na_count >= 2, ]
  
  melt.sub3[is.na(melt.sub3)] <- means[biomarker]    
  long <- melt(melt.sub3, variable.name = 'time', id.vars = c('sampleID', 'condition'))
  
  # Anova
  x <- aov(data=long, formula = value ~ condition + time + condition * time + Error(sampleID))
  y <- unlist(summary(x))  
  summary(x)
  pvals[which(pvals$biomarker == biomarker), 'pval.condition'] <- y['Error: sampleID.Pr(>F)1']
  pvals[which(pvals$biomarker == biomarker), 'pval.time'] <- y['Error: Within.Pr(>F)1']
  pvals[which(pvals$biomarker == biomarker), 'pval.interaction'] <- y['Error: Within.Pr(>F)2']  
}

# Nr of sigs
pvals %>%  filter(pval.time < 0.05) %>% nrow()
pvals %>% filter(pval.condition < 0.05) %>% nrow()
pvals %>% filter(pval.interaction < 0.05) %>% nrow()

#### Replication in Breda cohort
# How many of these can we replicate in the Breda cohort? First we check how many present in Breda then attempt to replicate
sig.condition <- pvals %>% filter(pval.condition < 0.05) %>% pull(biomarker)
sig.time <- pvals %>% filter(pval.time < 0.05) %>% pull(biomarker)
sig.interaction <- pvals %>% filter(pval.interaction < 0.05) %>% pull(biomarker)

# Take only the ones present in Breda cohort
condition.breda <- sig.condition[which(sig.condition %in% colnames(breda))]
time.breda <- sig.time[which(sig.time %in% colnames(breda))]
interaction.breda <-sig.interaction[which(sig.interaction %in% colnames(breda))]

df <- cbind(breda, annot.breda %>% dplyr::select(timepoint, condition, patient_nr))
melt.df <- melt(df, id.vars = c('timepoint', 'condition', 'patient_nr'))
melt.df <- na.omit(melt.df)
samples <- melt.df %>% filter(timepoint == 'T2') %>% pull(patient_nr)
means <- colMeans(breda)

melt.df %<>% 
  na.omit() %>% 
  filter(timepoint != 'T6') %>% 
  filter(timepoint %in% c('T1', 'T2', 'T3')) %>% 
  filter(patient_nr %in% samples)

condition.breda.sig <- c()
time.breda.sig <- c()
interaction.breda.sig <- c()


# Condition 
for (biomarker in condition.breda) {
  sub <- melt.df %>% filter(variable == biomarker)
  
  # There are a few duplicates. Remove them for now?
  sub2 <- sub %>% group_by(timepoint, patient_nr) %>% mutate(dupl = n() > 1)
  sub2 %<>% filter(dupl == FALSE) %>% dplyr::select(-dupl)
  # make wide df again so we can quickly replace all the NA
  melt.sub2 <- reshape2::dcast(sub2, patient_nr + condition ~ timepoint, value.var = 'value')
  
  melt.sub2[is.na(melt.sub2)] <- means[biomarker]    
  long <- melt(melt.sub2, variable.name = 'time', id.vars = c('patient_nr', 'condition'))
  x <- aov(data=long, formula = value ~ condition + time + condition * time + Error(patient_nr))
  
  y <- unlist(summary(x))  
  
  pval.condition <- y['Error: patient_nr.Pr(>F)1']
  # pval.time <- y['Error: Within.Pr(>F)1']
  # pval.interaction <- y['Error: Within.Pr(>F)2']
  
  if( pval.condition < PVALUE_CUTOFF) { condition.breda.sig <- c(condition.breda.sig, biomarker)}
  # if( pval.time < PVALUE_CUTOFF) { time.breda.sig <- c(time.breda.sig, biomarker)}
  # if( pval.interaction < PVALUE_CUTOFF) { interaction.breda.sig <- c(interaction.breda.sig, biomarker)}
  
}

# Time
for (biomarker in time.breda) {
  sub <- melt.df %>% filter(variable == biomarker)
  
  # There are a few duplicates. Remove them for now?
  sub2 <- sub %>% group_by(timepoint, patient_nr) %>% mutate(dupl = n() > 1)
  sub2 %<>% filter(dupl == FALSE) %>% dplyr::select(-dupl)
  # make wide df again so we can quickly replace all the NA
  melt.sub2 <- reshape2::dcast(sub2, patient_nr + condition ~ timepoint, value.var = 'value')
  
  melt.sub2[is.na(melt.sub2)] <- means[biomarker]    
  long <- melt(melt.sub2, variable.name = 'time', id.vars = c('patient_nr', 'condition'))
  x <- aov(data=long, formula = value ~ condition + time + condition * time + Error(patient_nr))
  
  y <- unlist(summary(x))  
  
  # pval.condition <- y['Error: patient_nr.Pr(>F)1']
  pval.time <- y['Error: Within.Pr(>F)1']
  # pval.interaction <- y['Error: Within.Pr(>F)2']
  
  # if( pval.condition < PVALUE_CUTOFF) { condition.breda.sig <- c(condition.breda.sig, biomarker)}
  if( pval.time < PVALUE_CUTOFF) { time.breda.sig <- c(time.breda.sig, biomarker)}
  # if( pval.interaction < PVALUE_CUTOFF) { interaction.breda.sig <- c(interaction.breda.sig, biomarker)}
  
}

# Interaction
for (biomarker in interaction.breda) {
  sub <- melt.df %>% filter(variable == biomarker)
  
  # There are a few duplicates. Remove them for now?
  sub2 <- sub %>% group_by(timepoint, patient_nr) %>% mutate(dupl = n() > 1)
  sub2 %<>% filter(dupl == FALSE) %>% dplyr::select(-dupl)
  # make wide df again so we can quickly replace all the NA
  melt.sub2 <- reshape2::dcast(sub2, patient_nr + condition ~ timepoint, value.var = 'value')
  
  melt.sub2[is.na(melt.sub2)] <- means[biomarker]    
  long <- melt(melt.sub2, variable.name = 'time', id.vars = c('patient_nr', 'condition'))
  x <- aov(data=long, formula = value ~ condition + time + condition * time + Error(patient_nr))
  
  y <- unlist(summary(x))  
  
  # pval.condition <- y['Error: patient_nr.Pr(>F)1']
  # pval.time <- y['Error: Within.Pr(>F)1']
  pval.interaction <- y['Error: Within.Pr(>F)2']
  
  # if( pval.condition < PVALUE_CUTOFF) { condition.breda.sig <- c(condition.breda.sig, biomarker)}
  # if( pval.time < PVALUE_CUTOFF) { time.breda.sig <- c(time.breda.sig, biomarker)}
  if( pval.interaction < PVALUE_CUTOFF) { interaction.breda.sig <- c(interaction.breda.sig, biomarker)}
  
}

# How many?
table(time.breda %in% time.breda.sig)
table(condition.breda %in% condition.breda.sig)
table(interaction.breda %in% interaction.breda.sig)

# Replication at the nominal pvalue levels.
# Afterwards, adjust Pvalues using BH
pvals$pval.time <- p.adjust(p = pvals$pval.time, method = 'BH')
pvals$pval.condition <- p.adjust(p = pvals$pval.condition, method=  'BH')
pvals$pval.interaction <- p.adjust(p = pvals$pval.interaction, method = 'BH')

# Export all proteins and their pvalues
write.csv(pvals, file = 'output/pvalues_anova.csv')
