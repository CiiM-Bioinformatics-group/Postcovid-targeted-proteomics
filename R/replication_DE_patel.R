rm(list = ls())
try(dev.off())

library(openxlsx)
library(dplyr)
library(magrittr)
library(ggplot2)

setwd('/Users/martijnzoodsma/Documents/PhD/corona/Postcovid-targeted-proteomics/')
conv <- read.xlsx('data/convert_olink_names.xlsx')

# Patel data
patel <- read.csv('data/Patel et al//tabula-supplementary information Patel et al 2021 _ subset.csv', header=T, sep = ',')
severe <- patel %>% filter(Comparison == "Control_v_Severe") 
critical <- patel %>% filter(Comparison == "Control_v_Critical") 

# Own DE results
ICU <- openxlsx::read.xlsx('output/DE.xlsx', sheet = 'healthy_vs_ICU')
nonICU <- openxlsx::read.xlsx('output/DE.xlsx', sheet = 'healthy_vs_nonICU')


# ICU vs critical
# How many of the significant proteins in ICU can we replicate significantly and with the same direction in Patel data?
ICU %<>% filter(OlinkID %in% critical$OLINK.ID) %>% arrange(match(OlinkID, critical$OLINK.ID))
critical %<>% filter(OLINK.ID %in% ICU$OlinkID) %>% arrange(match(OLINK.ID, ICU$OlinkID))
all(ICU$OlinkID == critical$OLINK.ID)

# How many same direction? Without significance threshold
table(ICU$logFC < 0, critical$logFC < 0)

ICU %<>% filter(significance == TRUE)
critical %<>% filter(OLINK.ID %in% ICU$OlinkID) %>% arrange(match(OLINK.ID, ICU$OlinkID))

table(critical$adj.P.Val < 0.05)
table(critical$P.Value < 0.05)

# Make df for export
ICU_replication <- data.frame(
  'protein' = ICU$Assay, 
  'OlinkID' = ICU$OlinkID, 
  'logFC.radboud' = ICU$logFC, 
  'adj.pval.radboud' = ICU$adj.P.Val, 
  'logFC.patel' = critical$logFC, 
  'adj.pval.patel' = critical$adj.P.Val
)

ICU_replication$same.dir <- ifelse(sign(ICU_replication$logFC.radboud) == sign(ICU_replication$logFC.patel), TRUE, FALSE)
ICU_replication$both_significant <- ifelse(ICU_replication$adj.pval.radboud < 0.05 & ICU_replication$adj.pval.patel < 0.05, TRUE, FALSE)

#################
# non-ICU vs severe
# How many of the significant proteins in ICU can we replicate significantly and with the same direction in Patel data?
nonICU %<>% filter(OlinkID %in% severe$OLINK.ID) %>% arrange(match(OlinkID, severe$OLINK.ID))
severe %<>% filter(OLINK.ID %in% nonICU$OlinkID) %>% arrange(match(OLINK.ID, severe$OLINK.ID))
all(nonICU$OlinkID == severe$OLINK.ID)

# How many same direction without significance threshold?
table(nonICU$logFC < 0, severe$logFC < 0)

nonICU %<>% filter(significance == TRUE)
severe %<>% filter(OLINK.ID %in% nonICU$OlinkID) %>% arrange(match(OLINK.ID, nonICU$OlinkID))

table(severe$adj.P.Val < 0.05)
table(severe$P.Value < 0.05)

# Make df for export
nonICU_replication <- data.frame(
  'protein' = nonICU$Assay, 
  'OlinkID' = nonICU$OlinkID, 
  'logFC.radboud' = nonICU$logFC, 
  'adj.pval.radboud' = nonICU$adj.P.Val, 
  'logFC.patel' = severe$logFC, 
  'adj.pval.patel' = severe$adj.P.Val
)

nonICU_replication$same.dir <- ifelse(sign(nonICU_replication$logFC.radboud) == sign(nonICU_replication$logFC.patel), TRUE, FALSE)
nonICU_replication$both_significant <- ifelse(nonICU_replication$adj.pval.radboud < 0.05 & nonICU_replication$adj.pval.patel < 0.05, TRUE, FALSE)

exp <- list(
  'ICU_replication' = ICU_replication, 
  'nonICU_replication' = nonICU_replication
)

write.xlsx('output/replication_DE_patel_results.xlsx', x = exp)
