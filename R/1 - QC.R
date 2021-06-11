# PCA plots for QC. Do any of the samples stand out?

library(dplyr)
library(magrittr)
library(ggplot2)
library(ggsci)
library(ggpubr)

rm(list = ls())
load('data/data.RData')

# PCA on each of the cohorts separately
# MHH
mhh <- na.omit(mhh)
annot.mhh <- annot.mhh[which(rownames(annot.mhh) %in% rownames(mhh)), ]

all(rownames(mhh) == rownames(annot.mhh))

pca <- prcomp(x = mhh, center = T, scale. = T)
pcs <- data.frame(pca$x)

pdf('output/pca_MHH_cohort.pdf', width = 5, height = 3)
ggplot() +
  geom_point(data = pcs, aes(PC1, PC2, color = annot.mhh$condition)) +
  theme_classic() +
  scale_color_nejm() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = 'PC1', y = 'PC2', color='Condition', title = 'MHH cohort')
dev.off()

# Radboud. Only T1
radboud <- na.omit(radboud)
annot.radboud %<>% filter(time == 'W1T1')
annot.radboud <- annot.radboud[which(rownames(annot.radboud) %in% rownames(radboud)), ]
radboud <- radboud[which(rownames(radboud) %in% rownames(annot.radboud)), ]

all(rownames(radboud) == rownames(annot.radboud))

pca <- prcomp(x = radboud, center = T, scale. = T)
pcs <- data.frame(pca$x)

pdf('output/pca_radboud_cohort.pdf', width = 4, height = 3)
ggplot() +
  geom_point(data = pcs, aes(PC1, PC2, color = annot.radboud$condition)) +
  theme_classic() +
  scale_color_nejm() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = 'PC1', y = 'PC2', color='Condition', title = 'Radboud cohort')
dev.off()


# Breda
breda <- na.omit(breda)
annot.breda <- na.omit(annot.breda)
annot.breda %<>% filter(timepoint == 'T1')
breda <- breda[which(rownames(breda) %in% rownames(annot.breda)), ]

all(rownames(breda) == rownames(annot.breda))

pca <- prcomp(x = breda, center = T, scale. = T)
pcs <- data.frame(pca$x)

pdf('output/pca_breda_cohort.pdf', width = 4, height = 3)
ggplot() +
  geom_point(data = pcs, aes(PC1, PC2, color = annot.breda$condition)) +
  theme_classic() +
  scale_color_nejm() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = 'PC1', y = 'PC2', color='Condition', title = 'Breda cohort')
dev.off()


## Remove apparent outliers
# Radboud: 1 sample
# MHH: 1 sample
# Breda: 0 sample
rm(list = ls())
load('data/data.RData')
rem <- c('TK_41', 'RUMC_0112')

annot.mhh <- annot.mhh[which(!rownames(annot.mhh) %in% rem), ]
mhh <- mhh[which(!rownames(mhh) %in% rem), ]

annot.radboud <- annot.radboud[which(!annot.radboud$sampleID %in% rem), ]
radboud <- radboud[which(rownames(radboud) %in% rownames(annot.radboud)), ]

rm(rem)
save.image(file = 'data/data.RData')




# PCA on the inflammatory panel across all three cohorts
rm(list = ls())
load('data/data.RData')

breda <- breda[, which(colnames(breda) %in% colnames(radboud)), ]


mhh %<>% select(all_of(colnames(breda)))
radboud %<>% select(all_of(colnames(breda)))

df <- rbind(mhh, radboud, breda) %>% na.omit()
annot <- rbind(annot.mhh %>% select(cohort, condition), 
               annot.radboud %>% filter(time == 'W1T1') %>% select(cohort, condition), 
               annot.breda %>% filter(timepoint == 'T1') %>% select(cohort, condition))

annot <- na.omit(annot)
df <- df[which(rownames(df) %in% rownames(annot)), ]
annot <- annot[which(rownames(annot) %in% rownames(df)), ]

all(rownames(df) == rownames(annot))

pca <- prcomp(x = df, center = T, scale. = T)
pcs <- data.frame(pca$x)

pdf('output/pca_inf_panel.pdf', width = 4, height = 3)
ggplot() +
  geom_point(data = pcs, aes(PC1, PC2, color = annot$condition)) +
  theme_classic() +
  scale_color_nejm() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = 'PC1', y = 'PC2', color='Condition', title = 'Inflammation panel\nBreda + MHH + Nijmegen')
dev.off()




# PCA separated by the panels. Only for Radboud and MHH cohorts
# Inflammation panel
rm(list = ls())
load('data/data.RData')

annot.radboud %<>% filter(time == 'W1T1')

annot <- rbind(annot.mhh %>% select(condition, cohort), 
               annot.radboud %>% select(condition, cohort))

df <- rbind(mhh %>% select(intersect(colnames(mhh), colnames(radboud))), 
            radboud %>% select(intersect(colnames(mhh), colnames(radboud))))

df <- na.omit(df)
annot <- na.omit(annot)

annot %<>% filter(rownames(annot) %in% intersect(rownames(df), rownames(annot)))
df %<>% filter(rownames(df) %in% intersect(rownames(df), rownames(annot)))

inf.proteins <- conv %>% filter(Olink.panel == 'Olink INFLAMMATION') %>% pull(OlinkID)

inf <- df[, which(colnames(df) %in% inf.proteins)]

pca <- prcomp(x = inf, center = T, scale. = T)
summary(pca)
plot.inf <- ggplot() +
  geom_point(data = data.frame(pca$x), aes(PC1, PC2, color = annot$condition)) +
  theme_classic() +
  scale_color_nejm() +
  labs(x = 'PC1', y = 'PC2', color='Condition', title = 'Inflammation panel') +
  theme(plot.title = element_text(hjust = 0.5))

# Cardiometabolic panel
cm.proteins <- conv %>% filter(Olink.panel == 'Olink CARDIOMETABOLIC') %>% pull(OlinkID)

cm <- df[, which(colnames(df) %in% cm.proteins)]

pca <- prcomp(x = cm, center = T, scale. = T)
summary(pca)
plot.cm <- ggplot() +
  geom_point(data = data.frame(pca$x), aes(PC1, PC2, color = annot$condition)) +
  theme_classic() +
  scale_color_nejm() +
  labs(x = 'PC1', y = 'PC2', color='Condition', title = 'Cardiometabolic panel') +
  theme(plot.title = element_text(hjust = 0.5))

# Cardiovascular II panel
cv.proteins <- conv %>% filter(Olink.panel == 'Olink CARDIOVASCULAR II') %>% pull(OlinkID)

cv <- df[, which(colnames(df) %in% cv.proteins)]

pca <- prcomp(x = cv, center = T, scale. = T)
summary(pca)
plot.cv <- ggplot() +
  geom_point(data = data.frame(pca$x), aes(PC1, PC2, color = annot$condition)) +
  theme_classic() +
  scale_color_nejm() +
  labs(x = 'PC1', y = 'PC2', color='Condition', title = 'Cardiovascular panel') +
  theme(plot.title = element_text(hjust = 0.5))


pdf('output/pca_panels_sep.pdf', width = 15, height = 5, onefile = F)
ggarrange(plot.inf, plot.cm, plot.cv, nrow = 1, common.legend = T, legend = 'top')
dev.off()


