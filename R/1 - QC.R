# PCA plots for QC. Do any of the samples stand out?

library(dplyr)
library(magrittr)
library(ggplot2)
library(ggsci)
library(ggpubr)

rm(list = ls())
setwd('/Users/martijnzoodsma/Documents/PhD/corona/final')
load('data/data.RData')

# PCA on the Radboud and MHH cohorts combined
annot.radboud %<>% filter(time == 'W1T1')

annot <- rbind(annot.mhh %>% select(condition, cohort), 
               annot.radboud %>% select(condition, cohort))

df <- rbind(mhh %>% select(intersect(colnames(mhh), colnames(radboud))), 
            radboud %>% select(intersect(colnames(mhh), colnames(radboud))))

df <- na.omit(df)
annot <- na.omit(annot)

annot %<>% filter(rownames(annot) %in% intersect(rownames(df), rownames(annot)))
df %<>% filter(rownames(df) %in% intersect(rownames(df), rownames(annot)))

pca <- prcomp(x = df, center = T, scale. = T)
pcs <- data.frame(pca$x)

annot$condition <- factor(annot$condition, levels = c('ICU', 'non-ICU', 'postcovid', 'healthy'))
pdf('output/pca.pdf', width = 5, height = 5)
ggplot() +
  geom_point(data = pcs, aes(PC1, PC2, color = annot$condition)) +
  theme_classic() +
  scale_color_nejm() +
  labs(x = 'PC1 [27.74%]', y = 'PC2 [11.66%]', color='Condition')
dev.off()

pcs[which(pcs$PC2 > 20), ]


## Remove two samples because outliers in the first plot
rm(list = ls())
load('data/data.RData')
rem <- c('TK_41', 'RUMC_0112')

annot.mhh <- annot.mhh[which(!rownames(annot.mhh) %in% rem), ]
mhh <- mhh[which(!rownames(mhh) %in% rem), ]

annot.radboud <- annot.radboud[which(!annot.radboud$sampleID %in% rem), ]
radboud <- radboud[which(rownames(radboud) %in% rownames(annot.radboud)), ]

rm(rem)
save.image(file = 'data/data.RData')



# PCA separated by the panels
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
  labs(x = 'PC1 [33.23%]', y = 'PC2 [13.73%]', color='Condition', title = 'Inflammation panel') +
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
  labs(x = 'PC1 [27.42%]', y = 'PC2 [20.52%]', color='Condition', title = 'Cardiometabolic panel') +
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
  labs(x = 'PC1 [29.13%]', y = 'PC2 [12.42%]', color='Condition', title = 'Cardiovascular panel') +
  theme(plot.title = element_text(hjust = 0.5))


pdf('output/pca_panels_sep.pdf', width = 15, height = 5, onefile = F)
ggarrange(plot.inf, plot.cm, plot.cv, nrow = 1, common.legend = T, legend = 'top')
dev.off()


