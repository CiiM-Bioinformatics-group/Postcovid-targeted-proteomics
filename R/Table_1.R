# Table 1 celltypes calculations
rm(list = ls())
dev.off()
library(dplyr)
library(magrittr)
# MHH cohort
setwd('/Users/martijnzoodsma/Documents/PhD/corona/Postcovid-targeted-proteomics/')
load('data/data.RData')

cellcounts <- annot.mhh %>% select(condition, CD45Lymphocytes, CD3Tcells, CD19Bcells, plasmablasts)
cellcounts <- cellcounts[complete.cases(cellcounts), ]
table(cellcounts$condition)

# Measurements are in cells / uL = cells / 10^-6 L.
# To be in cells / 10^-9 cells, divide by 10^3

# Lymphocytes
cellcounts %>% filter(condition == 'healthy') %>% pull(CD45Lymphocytes) %>% mean() / 1000
cellcounts %>% filter(condition == 'healthy') %>% pull(CD45Lymphocytes) %>% sd() / 1000

cellcounts %>% filter(condition == 'postcovid') %>% pull(CD45Lymphocytes)  %>% mean() / 1000
cellcounts %>% filter(condition == 'postcovid') %>% pull(CD45Lymphocytes)  %>% sd() / 1000

# T
cellcounts %>% filter(condition == 'healthy') %>% pull(CD3Tcells) %>% mean() / 1000
cellcounts %>% filter(condition == 'healthy') %>% pull(CD3Tcells)%>% sd() / 1000

cellcounts %>% filter(condition == 'postcovid') %>% pull(CD3Tcells) %>% mean() / 1000
cellcounts %>% filter(condition == 'postcovid') %>% pull(CD3Tcells)%>% sd() / 1000

# B
cellcounts %>% filter(condition == 'healthy') %>% pull(CD19Bcells) %>% mean() / 1000
cellcounts %>% filter(condition == 'healthy') %>% pull(CD19Bcells) %>% sd() / 1000

cellcounts %>% filter(condition == 'postcovid') %>% pull(CD19Bcells) %>% mean() / 1000
cellcounts %>% filter(condition == 'postcovid') %>% pull(CD19Bcells) %>% sd() / 1000

# Plasmablasts
cellcounts %>% filter(condition == 'healthy') %>% pull(plasmablasts) %>% mean() / 1000
cellcounts %>% filter(condition == 'healthy') %>% pull(plasmablasts) %>% sd() / 1000

cellcounts %>% filter(condition == 'postcovid') %>% pull(plasmablasts) %>% mean() / 1000
cellcounts %>% filter(condition == 'postcovid') %>% pull(plasmablasts) %>% sd() / 1000



# Radboud

library(openxlsx)
library(dplyr)
library(magrittr)

rm(list = ls())
load('data/data.RData')

annot <- read.xlsx('/Users/martijnzoodsma/Documents/PhD/corona/data/Nijmegen cohort/pat_info/Clinical data ALL RUMC patients (incl transfers) Castor export 220720_Correct timepoints RUMC plasma v3.xlsx')
annot <- cbind(reshape2::colsplit(annot$record_id, 
                                  pattern = '=', 
                                  names = c('RUMC', 'BERN')) %>% 
                 select(RUMC) %>% mutate(RUMC = gsub(' ', '', RUMC)), annot)

samples <- unique(annot.radboud$pat_id)
annot %<>% filter(RUMC %in% samples)

annot.radboud <- annot.radboud %>% 
  filter(time == 'W1T1') %>% 
  arrange(match(pat_id, samples))

annot %<>% filter(RUMC %in% annot.radboud$pat_id) %>% arrange(match(RUMC, annot.radboud$pat_id)) 

# How many patients with presumed vs proven COVID-19 split by condition?
table(annot$w1_tp1_location, annot$COVID_confirmed)

cellcounts <- annot %>% select(colnames(annot)[grepl('adm_diff', colnames(annot))])
cellcounts <- cbind(cellcounts, annot %>% select(RUMC))
colnames(cellcounts) <- c('neutrophils', 'lymphocytes', 'monocytes', 'eosinophils', 'basophils', 'sample')
cellcounts[cellcounts == 'Not done'] <- NA

apply(cellcounts, 2, is.na) %>% colSums()
cellcounts <- cellcounts[complete.cases(cellcounts), ]
apply(cellcounts, 2, is.na) %>% colSums()

cellcounts$neutrophils <- as.double(cellcounts$neutrophils)
cellcounts$lymphocytes <- as.double(cellcounts$lymphocytes)
cellcounts$monocytes <- as.double(cellcounts$monocytes)
cellcounts$eosinophils <- as.double(cellcounts$eosinophils)
cellcounts$basophils <- as.double(cellcounts$basophils) 

df <- cbind(cellcounts, 
            annot.radboud %>% filter(pat_id %in% cellcounts$sample) %>% arrange(match(pat_id, cellcounts$sample)) %>% select(condition))
table(df$condition)
# Lymphocytes
df %>% filter(condition == 'ICU') %>% pull(lymphocytes) %>% na.omit() %>% mean()
df %>% filter(condition == 'ICU') %>% pull(lymphocytes) %>% na.omit() %>% sd()

df %>% filter(condition == 'non-ICU') %>% pull(lymphocytes) %>% na.omit() %>% mean()
df %>% filter(condition == 'non-ICU') %>% pull(lymphocytes) %>% na.omit() %>% sd()

# monocytes
df %>% filter(condition == 'ICU') %>% pull(monocytes) %>% na.omit() %>% mean()
df %>% filter(condition == 'ICU') %>% pull(monocytes) %>% na.omit() %>% sd()

df %>% filter(condition == 'non-ICU') %>% pull(monocytes) %>% na.omit() %>% mean()
df %>% filter(condition == 'non-ICU') %>% pull(monocytes) %>% na.omit() %>% sd()

# neutrophils
df %>% filter(condition == 'ICU') %>% pull(neutrophils) %>% na.omit() %>% mean()
df %>% filter(condition == 'ICU') %>% pull(neutrophils) %>% na.omit() %>% sd()

df %>% filter(condition == 'non-ICU') %>% pull(neutrophils) %>% na.omit() %>% mean()
df %>% filter(condition == 'non-ICU') %>% pull(neutrophils) %>% na.omit() %>% sd()

# eosinophils
df %>% filter(condition == 'ICU') %>% pull(eosinophils) %>% na.omit() %>% mean()
df %>% filter(condition == 'ICU') %>% pull(eosinophils) %>% na.omit() %>% sd()

df %>% filter(condition == 'non-ICU') %>% pull(eosinophils) %>% na.omit() %>% mean()
df %>% filter(condition == 'non-ICU') %>% pull(eosinophils) %>% na.omit() %>% sd()

# basophils
df %>% filter(condition == 'ICU') %>% pull(basophils) %>% na.omit() %>% mean()
df %>% filter(condition == 'ICU') %>% pull(basophils) %>% na.omit() %>% sd()

df %>% filter(condition == 'non-ICU') %>% pull(basophils) %>% na.omit() %>% mean()
df %>% filter(condition == 'non-ICU') %>% pull(basophils) %>% na.omit() %>% sd()

## Tests
# Divide the MHH counts by 1000 so we arrive at the same unit
wilcox.test(x = annot.mhh %>% filter(condition == 'healthy') %>% pull(CD45Lymphocytes) %>% na.omit() / 1000, 
            y = df %>% filter(condition == 'ICU') %>% pull(lymphocytes) %>% na.omit())

wilcox.test(x = annot.mhh %>% filter(condition == 'healthy') %>% pull(CD45Lymphocytes) %>% na.omit() / 1000, 
            y = df %>% filter(condition == 'non-ICU') %>% pull(lymphocytes) %>% na.omit())

wilcox.test(x = annot.mhh %>% filter(condition == 'healthy') %>% pull(CD45Lymphocytes) %>% na.omit() / 1000, 
            y = annot.mhh %>% filter(condition == 'postcovid') %>% pull(CD45Lymphocytes) %>% na.omit() / 1000)

