rm(list = ls())
try(dev.off())

load('data/data.RData')

library(dplyr)
library(magrittr)
library(ggplot2)
library(ggsci)
library(RColorBrewer)
library(openxlsx)
library(pheatmap)

# Define the colours we want to use for ICU, non-ICU, postcovid and healthy
cols <- c('ICU' = '#BC3C29FF', 'non-ICU' = '#E18727FF', 'post-COVID-19' = '#0072B5FF', 'healthy' = '#1fc600')
theme_set(theme_classic() + theme(text = element_text(size = 12)))

# annot radboud unique patients
annot.radboud %<>% filter(time == 'W1T1')

annot <- rbind(annot.mhh %>% arrange(time.convalescence) %>% select(age, gender, condition), 
            annot.radboud %>% select(age, gender, condition))

# Heatmap over Radboud and MHH cohorts
common <- intersect(colnames(mhh), colnames(radboud))
annot.radboud %<>% filter(time == 'W1T1')

df <- rbind(
  mhh %>% select(all_of(common)), 
  radboud %>% select(all_of(common))
)

df <- na.omit(df)
annot <- annot[which(rownames(annot) %in% rownames(df)), ]

annot <- na.omit(annot)
df <- df[which(rownames(df) %in% rownames(annot)), ]

annot$condition[which(annot$condition == 'postcovid')] <- 'post-COVID-19'
annot$condition <- factor(annot$condition, levels = c('ICU', 'non-ICU', 'post-COVID-19', 'healthy'))
annot %<>% arrange(condition)
df %<>% arrange(match(rownames(df), rownames(annot)))

all(rownames(df) == rownames(annot))

# PCA using these samples
pca.res <- prcomp(df, center = T, scale. = T)
summary(pca.res)

pdf('output/pca_figure_1.pdf', width = 6, height = 4)
ggplot(data = data.frame(pca.res$x)) +
  geom_point(aes(x = PC1, y = PC2, color = annot$condition)) +
  # scale_color_manual(values = cols) + labs(x = 'PC1 [28.31%]', y = 'PC2 [10.62]', color = 'Condition') +
  theme(legend.position = 'top')
dev.off()

# Retrieve the healthy / postcovid samples and run PCA again
annot %<>% filter(condition %in% c('post-COVID-19', 'healthy'))
df <- df[which(rownames(df) %in% rownames(annot)), ]
pca.res <- prcomp(df, center = T, scale. = T)


pdf('output/pca_healthy_postcovid.pdf', width = 6, height = 4)
ggplot(data = data.frame(pca.res$x)) +
  geom_point(aes(x = PC1, y = PC2, color = annot$condition)) +
  scale_color_manual(values = cols) + labs(x = 'PC1 [18.6%]', y = 'PC2 [8.6%]', color = 'Condition') +
  theme(legend.position = 'none')
dev.off()

x <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)


# Heatmap of the differentially expressed genes
load('data/data.RData')

DE <- unique(c(
  read.xlsx('output/biomarkers.xlsx', sheet = 'ICU', colNames=F) %>% pull(X1), 
  read.xlsx('output/biomarkers.xlsx', sheet = 'nonICU', colNames=F) %>% pull(X1),
  read.xlsx('output/biomarkers.xlsx', sheet = 'postcovid', colNames=F) %>% pull(X1)
))

df %<>% select(all_of(DE))

pdf('output/heatmap_figure_1.pdf', width = 8, height = 4)
ph <- pheatmap::pheatmap(mat = t(df), 
         scale = 'row',
         show_colnames = F, 
         cluster_cols = F,
         annotation_col = annot %>% select(condition), 
         breaks = seq(-2, 2, length.out = 100),
         show_rownames = F, 
         annotation_colors = list(condition = cols), 
         annotation_names_col = F, 
         fontsize = 8,
         color = x, 
         cutree_rows = 4)
dev.off()


# Export excel sheet with clusters
prots <- sort(cutree(ph$tree_row, k=4))
names(prots) <- conv %>% filter(OlinkID %in% names(prots)) %>% arrange(match(OlinkID, names(prots))) %>% pull(Assay)
prots
table(prots)

df2 <- list('1' = names(prots[which(prots == 1)]), 
                  '2' = names(prots[which(prots == 3)]), 
                  '3' = names(prots[which(prots == 2)]), 
                  '4' = names(prots[which(prots == 4)]))
write.xlsx('output/clusters_heatmap.xlsx', x = df2)
