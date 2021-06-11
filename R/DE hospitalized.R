# Differential expression analysis within the Breda and Radboud cohorts separately
# COVID-19 ICU vs COVID-19 non-ICU 
# Perform DE in both cohorts. Use Breda as discovery and Radboud as replication

rm(list = ls())
try(dev.off())

library(dplyr)
library(magrittr)
library(openxlsx)
library(ggpubr)
library(limma)
library(ggplot2)
library(RColorBrewer)

theme_set(theme_classic() + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5)))
# All plots 4x4

cols <- c('ICU' = '#BC3C29FF', 'non-ICU' = '#E18727FF')

load('data/data.RData')

# Radboud only uses T1
annot.radboud %<>% filter(time == 'W1T1')
radboud <- radboud[which(rownames(radboud) %in% rownames(annot.radboud)), ]

# PCA
annot <- rbind(
    annot.radboud %>% select(age, gender, condition, cohort), 
    annot.breda %>% select(age, gender, condition, cohort)
  ) %>% na.omit()

df <- rbind(
  radboud %>% select(intersect(colnames(radboud), colnames(breda))), 
  breda %>% select(intersect(colnames(breda), colnames(radboud)))
) %>% na.omit()

df <- df[which(rownames(df) %in% rownames(annot)), ]
annot <- annot[which(rownames(annot) %in% rownames(df)), ]
all(rownames(annot) == rownames(df))

pca <- prcomp(df, center = T, scale. = T)
summary(pca)$importance[1:3, 1:3]

pdf('output/pca_breda_radboud.pdf', width = 5, height = 4)
ggplot() +
  geom_point(data = data.frame(pca$x), aes(PC1, PC2, color = annot$condition)) + 
  scale_color_manual(values = cols) +
  labs(x ='PC1 [23.1%]', y = 'PC2 [9.2%]', color = 'Condition')
dev.off()


#### DE analysis
load('data/data.RData')

# Radboud only uses T1
annot.radboud %<>% filter(time == 'W1T1')
radboud <- radboud[which(rownames(radboud) %in% rownames(annot.radboud)), ]

##### Radboud ICU vs Radboud nonICU
sum(is.na(annot.radboud))
sum(is.na(radboud))

annot.radboud$condition <- factor(annot.radboud$condition, levels = c('non-ICU', 'ICU'))

# Limma model
design <- model.matrix(~ annot.radboud$condition + annot.radboud$age + annot.radboud$gender)
head(design)
colnames(design) <- c('Intercept', 'condition', 'age', 'gender_male')
fit <- lmFit(as.matrix(t(radboud)), design = design)
fit <- eBayes(fit)
ICU_vs_nonICU_radboud <- topTable(fit = fit, adjust.method = "BH", number = Inf, coef = 'condition')
ICU_vs_nonICU_radboud$OlinkID <- rownames(ICU_vs_nonICU_radboud)
ICU_vs_nonICU_radboud$significance <- ifelse(ICU_vs_nonICU_radboud$adj.P.Val < 0.05, TRUE, FALSE)
ICU_vs_nonICU_radboud$direction <- ifelse(ICU_vs_nonICU_radboud$logFC < 0.0, 'Downregulated', 'Upregulated')
ICU_vs_nonICU_radboud <- cbind(ICU_vs_nonICU_radboud, 
                               conv %>% 
                                 filter(OlinkID %in% ICU_vs_nonICU_radboud$OlinkID) %>% 
                                 arrange(match(OlinkID, ICU_vs_nonICU_radboud$OlinkID)) %>% 
                                 select(Assay, Uniprot.ID, Olink.panel))


###### Breda ICU vs Breda nonICU
sum(is.na(annot.breda))
sum(is.na(breda))

annot.breda %<>% select(condition, age, gender)
annot.breda <- na.omit(annot.breda)

annot.breda$condition <- factor(annot.breda$condition, levels = c('non-ICU', 'ICU'))
annot.breda$age <- as.numeric(annot.breda$age)

breda <- breda[which(rownames(breda) %in% rownames(annot.breda)), ]


# Limma model
design <- model.matrix(~ annot.breda$condition + annot.breda$age + annot.breda$gender)
head(design)
colnames(design) <- c('Intercept', 'condition', 'age', 'gender_male')
fit <- lmFit(as.matrix(t(breda)), design = design)
fit <- eBayes(fit)
ICU_vs_nonICU_breda <- topTable(fit = fit, adjust.method = "BH", number = Inf, coef = 'condition')
ICU_vs_nonICU_breda$OlinkID <- rownames(ICU_vs_nonICU_breda)
ICU_vs_nonICU_breda$significance <- ifelse(ICU_vs_nonICU_breda$adj.P.Val < 0.05, TRUE, FALSE)
ICU_vs_nonICU_breda$direction <- ifelse(ICU_vs_nonICU_breda$logFC < 0.0, 'Downregulated', 'Upregulated')
ICU_vs_nonICU_breda <- cbind(ICU_vs_nonICU_breda, 
                               conv %>% 
                                 filter(OlinkID %in% ICU_vs_nonICU_breda$OlinkID) %>% 
                                 arrange(match(OlinkID, ICU_vs_nonICU_breda$OlinkID)) %>% 
                                 select(Assay, Uniprot.ID, Olink.panel))

write.xlsx(x = ICU_vs_nonICU_breda, file = 'output/DE_ICU_vs_nonICU_breda.xlsx')
write.xlsx(x = ICU_vs_nonICU_radboud, file = 'output/DE_ICU_vs_nonICU_radboud.xlsx')

######### Visualization 

DE_Breda <- ggplot() +
    geom_point(data = ICU_vs_nonICU_breda %>% filter(significance == F) , aes(x = logFC, y = -log10(adj.P.Val)), color = 'lightgray', show.legend = F) +
    geom_point(data = ICU_vs_nonICU_breda %>% filter(significance == T) , aes(x = logFC, y = -log10(adj.P.Val), color = direction)) +
    labs(x =expression(paste(log[2], ' fold change')), y = expression(paste(-log[10], '(adj. p-value)')), color = '', title = 'ICU vs non-ICU\nCohort 1') +
    scale_color_manual(values = c('#0072B5FF', '#BC3C29FF')) +
    geom_hline(yintercept = -log10(0.05), lty = 'dotted') +
    geom_vline(xintercept = 0.0, lty = 'dotted') +
    xlim(c(-1.5, 2)) +
    ylim(c(0, 60)) +
    ggrepel::geom_label_repel(data = head(x = ICU_vs_nonICU_breda, 5), aes(label = Assay, x = logFC, y = -log10(adj.P.Val)), box.padding = 0.5, show.legend = F)
  
DE_Radboud <- ggplot() +
    geom_point(data = ICU_vs_nonICU_radboud %>% filter(significance == F) , aes(x = logFC, y = -log10(adj.P.Val)), color = 'lightgray', show.legend = F) +
    geom_point(data = ICU_vs_nonICU_radboud %>% filter(significance == T) , aes(x = logFC, y = -log10(adj.P.Val), color = direction)) +
    labs(x = expression(paste(log[2], ' fold change')), y = expression(paste(-log[10], '(adj. p-value)')), color = '', title = 'ICU vs non-ICU\nCohort 2') +
    scale_color_manual(values = c('#0072B5FF', '#BC3C29FF')) +
    geom_hline(yintercept = -log10(0.05), lty = 'dotted') +
    geom_vline(xintercept = 0.0, lty = 'dotted') +
    xlim(c(-1.5, 2)) +
    ylim(c(0, 60)) +
    ggrepel::geom_label_repel(data = head(x = ICU_vs_nonICU_radboud %>% arrange(adj.P.Val), 3), aes(label = Assay, x = logFC, y = -log10(adj.P.Val)), box.padding = 0.5, show.legend = F) 
  


pdf('output/DE_ICU_vs_nonICU_Breda_Nijmegen.pdf', width = 8, height = 4, onefile = F)
ggarrange(DE_Breda, DE_Radboud, nrow = 1, ncol=2, common.legend = T, legend='top') 
dev.off()


######### Replication analysis
common <- intersect(ICU_vs_nonICU_breda$OlinkID, ICU_vs_nonICU_radboud$OlinkID)
ICU_vs_nonICU_breda %<>% filter(OlinkID %in% common)
ICU_vs_nonICU_radboud %<>% filter(OlinkID %in% common)
ICU_vs_nonICU_radboud %<>% arrange(match(OlinkID, ICU_vs_nonICU_breda$OlinkID))

stopifnot(ICU_vs_nonICU_radboud$OlinkID == ICU_vs_nonICU_breda$OlinkID)


# How replicable?
df <- data.frame(
  'protein' = ICU_vs_nonICU_radboud$Assay, 
  'OlinkID' = ICU_vs_nonICU_radboud$OlinkID, 
  'UniprotID' = ICU_vs_nonICU_radboud$Uniprot.ID, 
  'Olink panel' = ICU_vs_nonICU_radboud$Olink.panel,
  'logFC.radboud' = ICU_vs_nonICU_radboud$logFC, 
  'pval.radboud' = ICU_vs_nonICU_radboud$P.Value,
  'adj.pval.radboud' = ICU_vs_nonICU_radboud$adj.P.Val, 
  'sig.radboud' = ifelse(ICU_vs_nonICU_radboud$adj.P.Val < 0.05, TRUE, FALSE),
  'logFC.breda' = ICU_vs_nonICU_breda$logFC, 
  'pvalue.breda' = ICU_vs_nonICU_breda$P.Value,
  'adj.pval.breda' = ICU_vs_nonICU_breda$adj.P.Val, 
  'sig.breda' = ifelse(ICU_vs_nonICU_breda$adj.P.Val < 0.05, TRUE, FALSE), 
  'Comparison' = 'ICU_vs_nonICU'
)

# We want to show only the proteins that are significant in Breda cohort -> Discovery
df %<>% filter(sig.breda == TRUE)

# 40 proteins significant in Breda
# Of these, we replicate 29.
table(df$adj.pval.radboud < 0.05)
table(df$logFC.radboud > 0, df$logFC.breda > 0)
cor.test(df$logFC.breda, df$logFC.radboud)

pdf('output/replication_DE_breda_radboud.pdf', width = 4.5, height = 3.5)
ggplot(data = df) +
  geom_vline(xintercept = 0, lty=2, color = 'darkgrey') +
  geom_hline(yintercept = 0, lty=2, color = 'darkgrey') +
  geom_point(aes(x = logFC.breda, y = logFC.radboud, color = sig.radboud)) +
  labs(x = expression(paste(log[2], ' fold change cohort 1')), 
       y = expression(paste(log[2], ' fold change cohort 2')), 
       color = 'Significance in\ncohort 2') + 
  scale_color_manual(values = c('darkgrey', 'black')) +
  theme(legend.position = c(.8, .2))
dev.off()




# Heatmap of the proteins we are able to replicate
library(pheatmap)
load('data/data.RData')

replicated <- df %>% filter(adj.pval.radboud < 0.05) %>% pull(OlinkID)
df %>% filter(adj.pval.radboud < 0.05)
# Radboud only uses T1
annot.radboud %<>% filter(time == 'W1T1')
radboud <- radboud[which(rownames(radboud) %in% rownames(annot.radboud)), ]

# Heatmap
annot <- rbind(
  annot.radboud %>% select(age, gender, condition, cohort), 
  annot.breda %>% select(age, gender, condition, cohort)
) %>% na.omit()

df <- rbind(
  radboud %>% select(intersect(colnames(radboud), colnames(breda))), 
  breda %>% select(intersect(colnames(breda), colnames(radboud)))
) %>% na.omit()

df <- df[which(rownames(df) %in% rownames(annot)), ]
annot <- annot[which(rownames(annot) %in% rownames(df)), ]
all(rownames(annot) == rownames(df))

df %<>% select(all_of(replicated))
colnames(df) <- conv %>% filter(OlinkID %in% colnames(df)) %>% arrange(match(OlinkID, colnames(df))) %>% pull(Assay)

annot %<>% arrange(condition)
df <- df[rownames(annot), ]
annot$cohort <- as.numeric(factor(annot$cohort, levels = c('Breda', 'Radboud')))
annot$cohort <- as.character(annot$cohort)


# pdf('output/heatmap_ICU_nonICU.pdf', width = 13, height = 4)
png('output/heatmap_ICU_nonICU.png', width = 13, height = 4, units = 'in', res = 700)
pheatmap::pheatmap(t(df), 
         scale = 'row', cluster_cols = T,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100), 
         annotation_col = annot %>% select(cohort, condition), clustering_distance_cols = 'correlation',
         show_colnames = F, 
         breaks = seq(-1.5, 1.5, length.out = 100), 
         annotation_colors = list(condition = cols, 
                                  cohort = c('1' = 'chartreuse3', '2' = 'purple')), fontsize = 8)
dev.off()
