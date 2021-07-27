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
cols <- c('ICU' = '#BC3C29FF', 'non-ICU' = '#E18727FF', 'convalescent' = '#0072B5FF', 'healthy' = '#1fc600')
theme_set(theme_classic() + theme(text = element_text(size = 12)))

# Only first timepoint
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

annot$condition[which(annot$condition == 'postcovid')] <- 'convalescent'
annot$condition <- factor(annot$condition, levels = c('ICU', 'non-ICU', 'convalescent', 'healthy'))
annot %<>% arrange(condition)
df %<>% arrange(match(rownames(df), rownames(annot)))

all(rownames(df) == rownames(annot))

# PCA using these samples
pca.res <- prcomp(df, center = T, scale. = T)
summary(pca.res)

pdf('output/pca_figure_1.pdf', width = 6, height = 4)
ggplot(data = data.frame(pca.res$x)) +
  geom_point(aes(x = PC1, y = PC2, color = annot$condition)) +
  scale_color_manual(values = cols) + labs(x = 'PC1 [28%]', y = 'PC2 [11%]', color = 'Condition') +
  theme(legend.position = 'top')
dev.off()

# Retrieve the healthy / postcovid samples and run PCA again
annot %<>% filter(condition %in% c('convalescent', 'healthy'))
df <- df[which(rownames(df) %in% rownames(annot)), ]
pca.res <- prcomp(df, center = T, scale. = T)

summary(pca.res)

pdf('output/pca_healthy_postcovid.pdf', width = 6, height = 4)
ggplot(data = data.frame(pca.res$x)) +
  geom_point(aes(x = PC1, y = PC2, color = annot$condition)) +
  scale_color_manual(values = cols) + labs(x = 'PC1 [19%]', y = 'PC2 [9%]', color = 'Condition') +
  theme(legend.position = 'none')
dev.off()

# Heatmap of the differentially expressed genes
load('data/data.RData')

ICU <- read.xlsx('output/DE.xlsx', sheet = 'ICU_vs_healthy')
nonICU <- read.xlsx('output/DE.xlsx', sheet = 'nonICU_vs_healthy')
postcovid <- read.xlsx('output/DE.xlsx', sheet = 'postcovid_vs_healthy') %>% filter(Olink.panel != 'Olink NEUROLOGY')

ICU %<>% filter(significance == T) %>% pull(OlinkID)
nonICU %<>% filter(significance == T) %>% pull(OlinkID)
postcovid %<>% filter(significance == T) %>% pull(OlinkID)

DE <- unique(c(ICU, nonICU, postcovid))

DE <- DE[which(DE %in% colnames(df))]
df %<>% select(all_of(DE))

# Heatmap Fig. 2
# pdf('output/heatmap_figure_2.pdf', width = 8, height = 4)
png('output/heatmap_figure_2.png', width = 13, height = 6, units = 'in', res = 700)
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
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100), 
         cutree_rows = 4)
ph
dev.off()

# Excel sheets with the clusters
prots <- cutree(ph$tree_row, k = 4)

# Reorder based on row order of pheatmap
prots <- prots[ph$tree_row$order]
# names(prots) <- conv %>% filter(OlinkID %in% names(prots)) %>% arrange(match(OlinkID, names(prots))) %>% pull(Assay)
prots
data.frame(
  cluster = prots, 
  protein = names(prots)
) -> df

write.xlsx('output/clusters_heatmap.xlsx', x = df)
