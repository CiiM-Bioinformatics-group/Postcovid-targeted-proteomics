# Table 1 celltypes calculations
rm(list = ls())
dev.off()
library(dplyr)
library(magrittr)
# MHH cohort
setwd('/Users/martijnzoodsma/Documents/PhD/corona/Postcovid-targeted-proteomics/')
load('data/data.RData')

annot.mhh %>% pull(time.convalescence) %>% na.omit() %>% summary()
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

annot.radboud %<>% filter(time == 'W1T1')
annot.radboud %>% group_by(condition) %>% summarise(mean = mean(age), sd = sd(age))
table(annot.radboud$gender, annot.radboud$condition)
annot.radboud %>% group_by(condition) %>% summarise(mean = mean(na.omit(BMI)), sd = sd(na.omit(BMI)))
table(annot.radboud$condition, annot.radboud$`PCR_proven_COVID-19`)

table(annot.radboud$condition, is.na(annot.radboud$BMI))
# Blood cell counts. Remove NAs
counts <- annot.radboud[]

x <- c('condition', colnames(annot.radboud)[grepl('bloodcell', colnames(annot.radboud))])
tmp <- annot.radboud[, x]
tmp <- na.omit(tmp) #106 / 143
colnames(tmp)

tmp %>% group_by(condition) %>% 
  summarise(mean = mean(White.bloodcell.differentiation.neutrophils.at.admission), 
            sd = sd(White.bloodcell.differentiation.neutrophils.at.admission))

tmp %>% group_by(condition) %>% 
  summarise(mean = mean(White.bloodcell.differentiation.monocytes.at.admission), 
            sd = sd(White.bloodcell.differentiation.monocytes.at.admission))

tmp %>% group_by(condition) %>% 
  summarise(mean = mean(White.bloodcell.differentiation.lymphocytes.at.admission), 
            sd = sd(White.bloodcell.differentiation.lymphocytes.at.admission))

tmp %>% group_by(condition) %>% 
  summarise(mean = mean(White.bloodcell.differentiation.eosinophils.at.admission), 
            sd = sd(White.bloodcell.differentiation.eosinophils.at.admission))

tmp %>% group_by(condition) %>% 
  summarise(mean = mean(White.bloodcell.differentiation.basophils.at.admission), 
            sd = sd(White.bloodcell.differentiation.basophils.at.admission))

# Breda
annot.breda %<>% filter(timepoint == 'T1')
annot.breda %<>% na.omit()
str(annot.breda)
annot.breda %>% group_by(condition) %>% summarise(mean = mean(as.numeric(age)), sd = sd(as.numeric(age)))
table(annot.breda$gender, annot.breda$condition)

# annot.breda %>% group_by(condition) %>% summarise(mean = mean(as.numeric(BMI)), sd = sd(na.omit(BMI)))
table(annot.breda$condition, annot.breda$`PCR_proven_COVID-19`)


# Compare demographic data
rm(list = ls())
load('data/data.RData')

breda <- read.xlsx('/Users/martijnzoodsma/Documents/PhD/corona/data/Breda/clean_patinfo_breda.xlsx')

annot.radboud %<>% filter(time == 'W1T1')
annot.breda %<>% filter(timepoint == 'T1')
annot.breda %<>% na.omit(annot.breda)

table(annot.radboud$condition)
table(annot.breda$condition)
table(annot.mhh$condition)

annot.breda <- merge(x = annot.breda, y = breda, by.x = 'patient_nr', by.y = 'number')
annot.breda$age <- annot.breda$age.x
annot.breda$gender <- annot.breda$gender.x

# Age, gender, BMI,
df <- rbind(
  annot.radboud %>% filter(time == 'W1T1') %>% mutate(cohort = 'radboud') %>% select(age, gender, condition, cohort), 
  annot.mhh %>% mutate(cohort = 'mhh') %>% select(age, gender, condition, cohort), 
  annot.breda %>% mutate(cohort = 'breda') %>% filter(timepoint=='T1') %>% select(age, gender, condition, cohort)
)
str(df)
df$age <- as.numeric(df$age)
df %>% group_by(condition) %>% summarise(mean = mean(age), sd = sd(age))

res.aov <- aov(age ~ condition, data = df)
summary(res.aov)

# Ages vs healthy
wilcox.test(
  x = df %>% filter(condition == 'ICU') %>% pull(age) %>% as.numeric(), 
  y = df %>% filter(condition == 'healthy') %>% pull(age) %>% as.numeric()
)

wilcox.test(
  x = df %>% filter(condition == 'non-ICU') %>% pull(age) %>% as.numeric(), 
  y = df %>% filter(condition == 'healthy') %>% pull(age) %>% as.numeric()
)

chisq.test(df$gender, df$condition)

df <- rbind(
  annot.radboud %>% filter(time == 'W1T1') %>% mutate(cohort = 'radboud') %>% select(condition, cohort, BMI), 
  annot.breda %>% mutate(cohort = 'breda') %>% filter(timepoint=='T1') %>% select(condition, cohort, BMI)
)

df %>% group_by(condition, cohort) %>% summarise(mean = mean(na.omit(BMI)), sd = sd(na.omit(BMI)), n = n())

res.aov <- aov(BMI ~ condition, data=df)                                    
summary(res.aov)





## Tests for cell counts
load('data/data.RData')
x = annot.mhh %>% filter(condition == 'healthy') %>% pull(CD45Lymphocytes) %>% na.omit() / 1000
y = annot.radboud %>% filter(condition == 'ICU') %>% pull(White.bloodcell.differentiation.lymphocytes.at.admission)
wilcox.test(x,y)
summary(x)
summary(y)

x = annot.mhh %>% filter(condition == 'healthy') %>% pull(CD45Lymphocytes) %>% na.omit() / 1000
y = annot.radboud %>% filter(condition == 'non-ICU') %>% pull(White.bloodcell.differentiation.lymphocytes.at.admission)
sd(na.omit(y))
wilcox.test(x,y)
summary(x)
summary(y)
head(tmp)
head(annot.radboud)
annot.radboud %>% pull(White.bloodcell.differentiation.lymphocytes.at.admission) %>% na.omit()
tmp %>% pull(White.bloodcell.differentiation.lymphocytes.at.admission) %>% na.omit()
head(tmp)
# Divide the MHH counts by 1000 so we arrive at the same unit
wilcox.test(x = tmp %>% filter(condition == 'ICU') %>% pull(White.bloodcell.differentiation.lymphocytes.at.admission),
            y = annot.mhh %>% filter(condition == 'healthy') %>% pull(CD45Lymphocytes) %>% na.omit() / 1000, alternative = 'less')
# p = 7e-5 that MHH counts are higher than RAD

wilcox.test(x = tmp %>% filter(condition == 'non-ICU') %>% pull(White.bloodcell.differentiation.lymphocytes.at.admission),
            y = annot.mhh %>% filter(condition == 'healthy') %>% pull(CD45Lymphocytes) %>% na.omit() / 1000, alternative = 'less')
# p = 2e-9 that MHH non-ICU counts are higher than RAD

wilcox.test(x = annot.mhh %>% filter(condition == 'postcovid') %>% pull(CD45Lymphocytes) %>% na.omit() / 1000,
            y = annot.mhh %>% filter(condition == 'healthy') %>% pull(CD45Lymphocytes) %>% na.omit() / 1000, alternative = 'greater')
