# Load data, load annotation.
# merge and clean. Save as .RData object for further use

rm(list = ls())
dev.off()

suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(magrittr))

# Set to the correct directory
setwd('/Users/martijnzoodsma/Documents/PhD/corona/Postcovid-targeted-proteomics/')
conv <- read.xlsx('data/convert_olink_names.xlsx')

# Load the data per cohort. 
# MHH
inf <- read.xlsx('data/MHH/inflammation_panel.xlsx', rowNames = T, colNames = T, sep.names = ' ')
cm <- read.xlsx('data/MHH/cardiometabolic_panel.xlsx', rowNames = T, colNames = T, sep.names = ' ')
cvd <- read.xlsx('data/MHH/cardiovascular_panel.xlsx', rowNames = T, colNames = T, sep.names = ' ')
neu <- read.xlsx('data/MHH/neurology_panel.xlsx', rowNames = T, colNames = T, sep.names = ' ')

incl <- unique(c(rownames(inf), rownames(cm), rownames(cvd), rownames(neu)))

# Map the names back to Olink IDs. These have the advantage to be unique which is easier
# and allows us to keep track of the duplicates without too much trouble
colnames(inf) <- conv %>% 
  filter(Olink.panel == 'Olink INFLAMMATION') %>% 
  filter(Assay %in% colnames(inf)) %>% 
  arrange(match(Assay, colnames(inf))) %>% 
  pull(OlinkID)

colnames(cm) <- conv %>% 
  filter(Olink.panel == 'Olink CARDIOMETABOLIC') %>% 
  filter(Assay %in% colnames(cm)) %>% 
  arrange(match(Assay, colnames(cm))) %>% 
  pull(OlinkID)

colnames(cvd) <- conv %>% 
  filter(Olink.panel == 'Olink CARDIOVASCULAR II') %>% 
  filter(Assay %in% colnames(cvd)) %>% 
  arrange(match(Assay, colnames(cvd))) %>% 
  pull(OlinkID)

colnames(neu) <- conv %>% 
  filter(Olink.panel == 'Olink NEUROLOGY') %>% 
  filter(Assay %in% colnames(neu)) %>% 
  arrange(match(Assay, colnames(neu))) %>% 
  pull(OlinkID)

# Merge in two steps
tmp <- merge(x = inf, y = neu, by = 'row.names', all = T)
rownames(tmp) <- tmp$Row.names
tmp$Row.names <- NULL

tmp <- merge(x = tmp, y = cm, by = 'row.names', all = T)
rownames(tmp) <- tmp$Row.names
tmp$Row.names <- NULL

tmp <- merge(x = tmp, y = cvd, by = 'row.names', all = T)
rownames(tmp) <- tmp$Row.names
tmp$Row.names <- NULL

mhh <- tmp

# Annotation
annot.mhh <- read.xlsx('data/MHH/patient_information_complete.xlsx', rowNames = T, colNames = T)
annot.mhh <- annot.mhh[which(rownames(annot.mhh) %in% rownames(mhh)), ]
annot.mhh <- annot.mhh[rownames(mhh), ]


# Radboud
inf <- read.xlsx('data/Radboud/merged_data/inflammation.xlsx', rowNames = T, colNames = T, sep.names = ' ')
cm <- read.xlsx('data/Radboud/merged_data/cardiometabolic.xlsx', rowNames = T, colNames = T, sep.names = ' ')
cvd <- read.xlsx('data/Radboud/merged_data/cardiovascular.xlsx', rowNames = T, colNames = T, sep.names = ' ')

colnames(inf) <- conv %>% 
  filter(Olink.panel == 'Olink INFLAMMATION') %>% 
  filter(Assay %in% colnames(inf)) %>% 
  arrange(match(Assay, colnames(inf))) %>% 
  pull(OlinkID)

colnames(cm) <- conv %>% 
  filter(Olink.panel == 'Olink CARDIOMETABOLIC') %>% 
  filter(Assay %in% colnames(cm)) %>% 
  arrange(match(Assay, colnames(cm))) %>% 
  pull(OlinkID)

colnames(cvd) <- conv %>% 
  filter(Olink.panel == 'Olink CARDIOVASCULAR II') %>% 
  filter(Assay %in% colnames(cvd)) %>% 
  arrange(match(Assay, colnames(cvd))) %>% 
  pull(OlinkID)

# merge in two steps
tmp <- merge(inf, cm, all=T, by = 'row.names')
rownames(tmp) <- tmp$Row.names
tmp$Row.names <- NULL

tmp <- merge(tmp, cvd, all=T, by = 'row.names')
rownames(tmp) <- tmp$Row.names
tmp$Row.names <- NULL

radboud <- tmp

# Annotation
annot.radboud <- read.xlsx('data/Radboud/patient_information.xlsx', rowNames = T, colNames = T)

info <- data.frame('full_id' = rownames(radboud))

info <- cbind(info,
              reshape2::colsplit(info$full_id, " ", c('pat_id', 'time'))) %>% select(-time)

annot.radboud <- merge(info, annot.radboud, by.x = 'full_id', by.y = 'row.names', all=T, sort=F)

annot.radboud <- annot.radboud[which(!is.na(annot.radboud$time)), ]
rownames(annot.radboud) <- annot.radboud$full_id

annot.radboud <- annot.radboud[which(rownames(annot.radboud) %in% rownames(radboud)), ]
annot.radboud %<>% arrange(match(rownames(annot.radboud), rownames(radboud)))

radboud <- radboud[which(rownames(radboud) %in% rownames(annot.radboud)), ]
radboud %<>% arrange(match(rownames(radboud), rownames(annot.radboud)))

# Breda
breda <- read.xlsx('data/Breda/integration_breda_MHH_inf.xlsx', sheet = 'Breda final', rowNames = T, colNames = T, sep.names = ' ')
colnames(breda) <- conv %>% 
  filter(Olink.panel == 'Olink INFLAMMATION') %>% 
  filter(Assay %in% colnames(breda)) %>% 
  arrange(match(Assay, colnames(cvd))) %>% 
  pull(OlinkID)

annot.breda <- read.xlsx('data/Breda/patient_info_breda.xlsx')
rownames(annot.breda) <- annot.breda$Olink.ID
annot.breda <-subset(annot.breda, select = -Olink.ID)

common <- intersect(rownames(breda), rownames(annot.breda))
breda <- breda[which(rownames(breda) %in% common), ]
annot.breda <- annot.breda[which(rownames(annot.breda) %in% common), ]
annot.breda %<>% arrange(match(rownames(breda), rownames(annot.breda)))

# Checks
stopifnot(all(rownames(mhh) == rownames(annot.mhh)))
stopifnot(all(rownames(radboud) == rownames(annot.radboud)))
stopifnot(all(rownames(breda) == rownames(annot.breda)))

# Normalise annotation
annot.breda %<>% 
  mutate(condition = ifelse(ICU == 'yes', 'ICU', 'non-ICU')) 

annot.mhh$cohort <- 'MHH'
annot.breda$cohort <- 'Breda'
annot.radboud$cohort <- 'Radboud'

rm(cm, cvd, inf, info, neu, tmp, common, incl)

# Remove the bridging samples from the MHH cohort. Not necessary anymore
annot.mhh %<>% filter(condition != 'bridge')
mhh <- mhh[which(rownames(mhh) %in% rownames(annot.mhh)), ]
all(rownames(mhh) == rownames(annot.mhh))

save.image(file = 'data/data.RData')
