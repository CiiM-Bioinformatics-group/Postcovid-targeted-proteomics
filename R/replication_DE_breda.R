# DE for breda vs MHH cohorts, Venn diagram and heatmaps for different genes
# ICU / non-ICU / postcovid
rm(list = ls())
setwd('/Users/martijnzoodsma/Documents/PhD/corona/Postcovid-targeted-proteomics/')
load('data/data.RData')

library(limma)
library(openxlsx)
library(dplyr)
library(magrittr)

# DE limma corrected for age and gender
# Only T1 for each patient
annot <- rbind(annot.mhh %>% select(age, gender, condition, cohort), 
               annot.breda %>% filter(timepoint == 'T1') %>% select(age, gender, condition, cohort))
annot$age <- as.numeric(annot$age)

df <- rbind(
  mhh %>% select(intersect(colnames(mhh), colnames(breda))), 
  breda %>% select(intersect(colnames(mhh), colnames(breda)))
)

df %<>% filter(rownames(df) %in% intersect(rownames(df), rownames(annot)))
annot %<>% filter(rownames(annot) %in% intersect(rownames(annot), rownames(df)))
stopifnot(all(rownames(df) == rownames(annot)))

# DE MHH healthy vs Breda ICU
sub.annot <- annot %>% filter(condition %in% c('healthy', 'ICU'))
sub.df <- df %>% filter(row.names(df) %in% rownames(sub.annot))
sub.annot$condition <- factor(sub.annot$condition, levels = c('healthy', 'ICU'))

# Limma model
design <- model.matrix(~ sub.annot$condition + sub.annot$age + sub.annot$gender)
colnames(design) <- c('Intercept', 'condition', 'age', 'gender_male')
fit <- lmFit(as.matrix(t(sub.df)), design = design)
fit <- eBayes(fit)
healthy_vs_ICU <- topTable(fit = fit, adjust.method = "BH", number = Inf, coef = 'condition')
healthy_vs_ICU$OlinkID <- rownames(healthy_vs_ICU)
healthy_vs_ICU$significance <- ifelse(healthy_vs_ICU$adj.P.Val < 0.05, TRUE, FALSE)
healthy_vs_ICU$direction <- ifelse(healthy_vs_ICU$logFC < 0.0, 'Downregulated', 'Upregulated')
healthy_vs_ICU <- cbind(healthy_vs_ICU, 
                        conv %>% 
                          filter(OlinkID %in% healthy_vs_ICU$OlinkID) %>% 
                          arrange(match(OlinkID, healthy_vs_ICU$OlinkID)) %>% 
                          select(Assay, Uniprot.ID, Olink.panel))


# DE MHH healthy vs Breda non-ICU
sub.annot <- annot %>% filter(condition %in% c('healthy', 'non-ICU'))
sub.df <- df %>% filter(row.names(df) %in% rownames(sub.annot))
sub.annot$condition <- factor(sub.annot$condition, levels = c('healthy', 'non-ICU'))

# Limma model
design <- model.matrix(~ sub.annot$condition + sub.annot$age + sub.annot$gender)
colnames(design) <- c('Intercept', 'condition', 'age', 'gender_male')
fit <- lmFit(as.matrix(t(sub.df)), design = design)
fit <- eBayes(fit)

healthy_vs_nonICU <- topTable(fit = fit, adjust.method = "BH", number = Inf, coef = 'condition')
healthy_vs_nonICU$OlinkID <- rownames(healthy_vs_nonICU)
healthy_vs_nonICU$significance <- ifelse(healthy_vs_nonICU$adj.P.Val < 0.05, TRUE, FALSE)
healthy_vs_nonICU$direction <- ifelse(healthy_vs_nonICU$logFC < 0.0, 'Downregulated', 'Upregulated')
healthy_vs_nonICU <- cbind(healthy_vs_nonICU, 
                           conv %>% 
                             filter(OlinkID %in% healthy_vs_nonICU$OlinkID) %>% 
                             arrange(match(OlinkID, healthy_vs_nonICU$OlinkID)) %>% 
                             select(Assay, Uniprot.ID, Olink.panel))


# DE MHH healthy vs MHH postcovid
# Already done this. Check results to ensure they are the same.
sub.annot <- annot %>% filter(condition %in% c('healthy', 'postcovid'))
sub.df <- df %>% filter(row.names(df) %in% rownames(sub.annot))
sub.annot$condition <- factor(sub.annot$condition, levels = c('healthy', 'postcovid'))

# Limma model
design <- model.matrix(~ sub.annot$condition + sub.annot$age + sub.annot$gender)
colnames(design) <- c('Intercept', 'condition', 'age', 'gender_male')
fit <- lmFit(as.matrix(t(sub.df)), design = design)
fit <- eBayes(fit)
healthy_vs_postcovid <- topTable(fit = fit, adjust.method = "BH", number = Inf, coef = 'condition')
healthy_vs_postcovid$OlinkID <- rownames(healthy_vs_postcovid)
healthy_vs_postcovid$significance <- ifelse(healthy_vs_postcovid$adj.P.Val < 0.05, TRUE, FALSE)
healthy_vs_postcovid$direction <- ifelse(healthy_vs_postcovid$logFC < 0.0, 'Downregulated', 'Upregulated')
healthy_vs_postcovid <- cbind(healthy_vs_postcovid, 
                              conv %>% 
                                filter(OlinkID %in% healthy_vs_postcovid$OlinkID) %>% 
                                arrange(match(OlinkID, healthy_vs_postcovid$OlinkID)) %>% 
                                select(Assay, Uniprot.ID, Olink.panel))

res <- list(healthy_vs_ICU = healthy_vs_ICU, 
            healthy_vs_nonICU = healthy_vs_nonICU, 
            healthy_vs_postcovid = healthy_vs_postcovid)

openxlsx::write.xlsx(x = res, file = 'output/DE_breda_MHH.xlsx')


# How comparable are MHH healthy vs Radboud ICU / MHH healthy vs Breda ICU?
# How comaprable are MHH healthy vs Radboud nonICU / MHH healthy vs Breda nonICU?

# ICU first
rm(list = ls())
radboud <- read.xlsx('output/DE.xlsx', sheet = 'healthy_vs_ICU')
breda <- read.xlsx('output/DE_breda_MHH.xlsx', sheet = 'healthy_vs_ICU')

radboud %<>% filter(OlinkID %in% breda$OlinkID) %>% arrange(match(OlinkID, breda$OlinkID))
breda %<>% filter(OlinkID %in% radboud$OlinkID) %>% arrange(match(OlinkID, radboud$OlinkID))
all(radboud$OlinkID == breda$OlinkID)

# How many of the significant proteins can we replicate significantly and with the same direction in Breda?
radboud %<>% filter(significance == TRUE)
breda %<>% filter(OlinkID %in% radboud$OlinkID) %>% arrange(match(OlinkID, radboud$OlinkID))

table(breda$significance)

# Make df for export
df <- data.frame(
  'protein' = radboud$Assay, 
  'OlinkID' = radboud$OlinkID, 
  'logFC.radboud' = radboud$logFC, 
  'adj.pval.radboud' = radboud$adj.P.Val, 
  'logFC.breda' = breda$logFC, 
  'adj.pval.breda' = breda$adj.P.Val
)

write.xlsx(x = df, file = 'output/replication_DE_breda_results.xlsx', sheet = 'ICU')



##
radboud <- read.xlsx('output/DE.xlsx', sheet = 'healthy_vs_nonICU')
breda <- read.xlsx('output/DE_breda_MHH.xlsx', sheet = 'healthy_vs_nonICU')

radboud %<>% filter(OlinkID %in% breda$OlinkID) %>% arrange(match(OlinkID, breda$OlinkID))
breda %<>% filter(OlinkID %in% radboud$OlinkID) %>% arrange(match(OlinkID, radboud$OlinkID))
all(radboud$OlinkID == breda$OlinkID)

# How many of the significant proteins can we replicate significantly and with the same direction in Breda?
radboud %<>% filter(significance == TRUE)
breda %<>% filter(OlinkID %in% radboud$OlinkID) %>% arrange(match(OlinkID, radboud$OlinkID))

table(breda$significance)

# Make df for export
df <- data.frame(
  'protein' = radboud$Assay, 
  'OlinkID' = radboud$OlinkID, 
  'logFC.radboud' = radboud$logFC, 
  'adj.pval.radboud' = radboud$adj.P.Val, 
  'logFC.breda' = breda$logFC, 
  'adj.pval.breda' = breda$adj.P.Val
)

write.xlsx(x = df, file = 'output/replication_DE_breda_results.xlsx', sheet = 'nonICU')



