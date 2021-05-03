rm(list=ls())
dev.off()

setwd("/Users/martijnzoodsma/Documents/PhD/corona/Postcovid-targeted-proteomics/")

suppressPackageStartupMessages(library(OlinkAnalyze))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))

THRESHOLD_MISSING_DATA <- 0.20

## Load in the data and combine
## NPX function from Olink produces dataframes in long format. Exclude samples that did not pass QC and proteins with too much missingness.
## Cast to wide df and export
cm.df <- read_NPX(filename = 'data/MHH/raw/MHH cohort analysis_NPX_CardioMetabolic_v1.xlsx')
cvd.df <- read_NPX(filename = 'data/MHH/raw/MHH cohort analysis_NPX_CVD II_v1.xlsx')
inf.df <- read_NPX(filename = 'data/MHH/raw/MHH cohort analysis_NPX_INF_v1.xlsx')
neu.df <- read_NPX(filename = 'data/MHH/raw/MHH cohort analysis_NPX_Neur_v1.xlsx')

cm.df %<>%
  mutate(MissingFreq = as.numeric(MissingFreq)) %>% 
  filter(QC_Warning == 'Pass') %>%
  filter(MissingFreq <= THRESHOLD_MISSING_DATA) %>% 
  reshape2::dcast(data = ., formula = SampleID ~ Assay, value.var = 'NPX')

cvd.df %<>%
  mutate(MissingFreq = as.numeric(MissingFreq)) %>% 
  filter(QC_Warning == 'Pass') %>% 
  filter(MissingFreq <= THRESHOLD_MISSING_DATA) %>% 
  reshape2::dcast(data = ., formula = SampleID ~ Assay, value.var = 'NPX')

inf.df %<>%
  mutate(MissingFreq = as.numeric(MissingFreq)) %>% 
  filter(QC_Warning == 'Pass') %>% 
  filter(MissingFreq <= THRESHOLD_MISSING_DATA) %>% 
  reshape2::dcast(data = ., formula = SampleID ~ Assay, value.var = 'NPX')

neu.df %<>%
  mutate(MissingFreq = as.numeric(MissingFreq)) %>% 
  filter(QC_Warning == 'Pass') %>% 
  filter(MissingFreq <= THRESHOLD_MISSING_DATA) %>% 
  reshape2::dcast(data = ., formula = SampleID ~ Assay, value.var = 'NPX')


# Export
write.xlsx(x = inf.df, file = 'data/MHH/inflammation_panel.xlsx')
write.xlsx(x = cm.df, file = 'data/MHH/cardiometabolic_panel.xlsx')
write.xlsx(x = cvd.df, file = 'data/MHH/cardiovascular_panel.xlsx')
write.xlsx(x = neu.df, file = 'data/MHH/neurology_panel.xlsx')
