rm(list = ls())
try(dev.off())

setwd('/Users/martijnzoodsma/Documents/PhD/corona/Postcovid-targeted-proteomics/')

library(dplyr)
library(magrittr)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(pheatmap)

load('data/data.RData')


annot.radboud %<>% filter(time %in% c('W1T1', 'W1T2', 'W1T3'))
radboud <- radboud[which(rownames(radboud) %in% rownames(annot.radboud)), ]

radboud <- na.omit(radboud)
annot.radboud <- annot.radboud[which(rownames(annot.radboud) %in% rownames(radboud)), ]

# 61 patients that were consistently samples over the first three timepoints
samples <- annot.radboud %>% group_by(sampleID) %>% summarise(n= n()) %>% filter(n >= 3) %>% pull(sampleID)
annot.radboud %<>% filter(sampleID %in% samples)
radboud <- radboud[which(rownames(radboud) %in% rownames(annot.radboud)), ]

all(rownames(radboud) == rownames(annot.radboud))

# Identify genes with highest variance
vars <- sapply(radboud, var)
genes <- names(vars[which(vars > 1)])
radboud %<>% select(all_of(genes))


radboud <- as.matrix(radboud)
radboud.mean.remove <- matrix(0, nrow = nrow(radboud), ncol = ncol(radboud))

# Do we need to scale in the PCA?
x <- reshape2::melt(radboud)
ggplot() +
  geom_density(data=x, aes(x = value, color = Var2)) +
  theme(legend.position = 'none')

x <- reshape2::melt(scale(radboud))
ggplot() +
  geom_density(data=x, aes(x = value, color = Var2)) +
  theme(legend.position = 'none')

# yes, we need to scale


# Remove mean for the both groups separately
for (i in seq(1, nrow(radboud), 3)) {
  beta_group_mean<- colMeans(radboud[i:(i+2),])
  radboud.mean.remove[i,]<- radboud[i,]- beta_group_mean
  radboud.mean.remove[i+1,]<- radboud[i+1,]- beta_group_mean
  radboud.mean.remove[i+2,]<- radboud[i+2,]- beta_group_mean
}

pca <- prcomp(x = radboud.mean.remove, center=TRUE, scale. = TRUE)
pcs <- data.frame(pca$x, row.names = rownames(radboud))

pcs <- cbind(pcs, annot.radboud)

# Identify the samples that go down from T1 to T2 and flip them
flip <- c()
for (sample in unique(pcs$sampleID)) {
  
  sub <- pcs %>% filter(sampleID == sample)
  
  # Is T2 smaller than T1 in PC1?
  T1 <- sub %>% filter(time == 'W1T1') %>% pull(PC1)
  T2 <- sub %>% filter(time == 'W1T2') %>% pull(PC1)
  
  if (T1 > T2) { flip <- c(flip, sample) }
}

# Flip the samples that need flipping
pcs[which(pcs$sampleID %in% flip), 'PC1'] <- -1 * pcs[which(pcs$sampleID %in% flip), 'PC1']

# Color by sampleID
pdf('output/anova_pca.pdf', width = 5, height = 3)
ggplot() +
  geom_point(data = pcs, aes(x = time, y = PC1, color = sampleID)) +
  geom_line(data=pcs, aes(x = time,y = PC1, group = sampleID, color = sampleID)) +
  theme_classic() +
  theme(legend.position = 'none', 
        axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45)) +
  scale_x_discrete(expand = c(0.05, 0.05)) +
  labs(x = 'Time', y = 'PC1')
dev.off()

# Color by condition. We see no difference between ICU and nonICU
ggplot() +
  geom_point(data = pcs, aes(x = time, y = PC1, color = condition)) +
  geom_line(data=pcs, aes(x = time,y = PC1, group = sampleID, color = condition)) +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_x_discrete(expand = c(0.05, 0.05))


