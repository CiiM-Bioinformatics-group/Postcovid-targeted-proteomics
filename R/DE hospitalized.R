# Differential expression analysis within the Breda and Radboud cohorts separately
# COVID-19 ICU vs COVID-19 non-ICU 
# Perform DE in both cohorts. Use Breda as discovery and Radboud as replication
# Only use T1 for each patient in both cohorts

rm(list = ls())
try(dev.off())

library(dplyr)
library(magrittr)
library(openxlsx)
library(ggpubr)
library(limma)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

theme_set(theme_classic() + theme(text = element_text(size = 10), plot.title = element_text(hjust = 0.5)))
cols <- c('ICU' = '#BC3C29FF', 'non-ICU' = '#E18727FF')

load('data/data.RData')

# Filter for only T1 and remove the points that we do not use in the matrices
annot.radboud %<>% filter(time == 'W1T1')
radboud <- radboud[which(rownames(radboud) %in% rownames(annot.radboud)), ]

annot.breda %<>% filter(timepoint == 'T1')
breda <- breda[which(rownames(breda) %in% rownames(annot.breda)), ]

##### Radboud ICU vs Radboud nonICU
annot.radboud %<>% select(age, gender, condition) %>% na.omit()
radboud <- radboud[which(rownames(radboud) %in% rownames(annot.radboud)), ]
annot.radboud$condition[which(annot.radboud$condition == 'non-ICU')] <- 'nonICU'
stopifnot(all(rownames(radboud) == rownames(annot.radboud)))

condition <- annot.radboud$condition
age <- annot.radboud$age
gender <- annot.radboud$gender

design <- model.matrix(~0 + condition + age + gender)
fit <- lmFit(as.matrix(t(radboud)), design = design)
cont <- makeContrasts(cond = 'conditionICU-conditionnonICU', levels = design)
fit <- contrasts.fit(fit, cont)
fit <- eBayes(fit)

ICU_vs_nonICU_radboud <- topTable(fit = fit, adjust.method = "BH", number = Inf)
ICU_vs_nonICU_radboud$OlinkID <- rownames(ICU_vs_nonICU_radboud)
ICU_vs_nonICU_radboud$significance <- ifelse(ICU_vs_nonICU_radboud$adj.P.Val < 0.05, TRUE, FALSE)
ICU_vs_nonICU_radboud$direction <- ifelse(ICU_vs_nonICU_radboud$logFC < 0.0, 'Downregulated', 'Upregulated')
ICU_vs_nonICU_radboud <- cbind(ICU_vs_nonICU_radboud, 
                               conv %>% 
                                 filter(OlinkID %in% ICU_vs_nonICU_radboud$OlinkID) %>% 
                                 arrange(match(OlinkID, ICU_vs_nonICU_radboud$OlinkID)) %>% 
                                 select(Assay, Uniprot.ID, Olink.panel))

###### Breda ICU vs Breda nonICU
annot.breda %<>% select(condition, age, gender) %>% na.omit()
breda <- breda[which(rownames(breda) %in% rownames(annot.breda)), ]
stopifnot(all(rownames(breda) == rownames(annot.breda)))
annot.breda$condition[which(annot.breda$condition == 'non-ICU')] <- 'nonICU'

condition <- annot.breda$condition
age <- annot.breda$age
gender <- annot.breda$gender

design <- model.matrix(~0 + condition + age + gender)
fit <- lmFit(as.matrix(t(breda)), design = design)
cont <- makeContrasts(cond = 'conditionICU-conditionnonICU', levels = design)
fit <- contrasts.fit(fit, cont)
fit <- eBayes(fit)

ICU_vs_nonICU_breda <- topTable(fit = fit, adjust.method = "BH", number = Inf)
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
# Since we focus on the inflammatory proteins here, we show only the Inflammation panel in the Radboud cohort, 
# While there are more proteins that were measured. We export the full protein list with results.

DE_Breda <- ggplot() +
  geom_point(data = ICU_vs_nonICU_breda %>% filter(significance == F) , aes(x = logFC, y = -log10(adj.P.Val)), color = 'lightgray', show.legend = F) +
  geom_point(data = ICU_vs_nonICU_breda %>% filter(significance == T) , aes(x = logFC, y = -log10(adj.P.Val), color = direction)) +
  labs(x =expression(paste(log[2], ' fold change')), y = expression(paste(-log[10], '(adj. p-value)')), color = '', title = 'ICU vs non-ICU\nCohort 1') +
  scale_color_manual(values = c('#0072B5FF', '#BC3C29FF')) +
  geom_hline(yintercept = -log10(0.05), lty = 'dotted') +
  geom_vline(xintercept = 0.0, lty = 'dotted') +
  xlim(c(-1.5, 2)) +
  ylim(c(0, 15)) +
  ggrepel::geom_label_repel(data = head(x = ICU_vs_nonICU_breda, 8), aes(label = Assay, x = logFC, y = -log10(adj.P.Val)), box.padding = 0.5, show.legend = F)

DE_Radboud <- ggplot() +
  geom_point(data = ICU_vs_nonICU_radboud %>% filter(Olink.panel == 'Olink INFLAMMATION') %>% filter(significance == F) , aes(x = logFC, y = -log10(adj.P.Val)), color = 'lightgray', show.legend = F) +
  geom_point(data = ICU_vs_nonICU_radboud %>% filter(Olink.panel == 'Olink INFLAMMATION') %>% filter(significance == T) , aes(x = logFC, y = -log10(adj.P.Val), color = direction)) +
  labs(x = expression(paste(log[2], ' fold change')), y = expression(paste(-log[10], '(adj. p-value)')), color = '', title = 'ICU vs non-ICU\nCohort 2') +
  scale_color_manual(values = c('#0072B5FF', '#BC3C29FF')) +
  geom_hline(yintercept = -log10(0.05), lty = 'dotted') +
  geom_vline(xintercept = 0.0, lty = 'dotted') +
  xlim(c(-1.5, 2)) +
  ylim(c(0, 15)) +
  ggrepel::geom_label_repel(data = head(x = ICU_vs_nonICU_radboud %>% arrange(adj.P.Val), 8), aes(label = Assay, x = logFC, y = -log10(adj.P.Val)), box.padding = 0.5, show.legend = F) 

pdf('output/DE_ICU_vs_nonICU_Breda_Nijmegen.pdf', width = 8, height = 4, onefile = F)
ggarrange(DE_Breda, DE_Radboud, nrow = 1, ncol=2, common.legend = T, legend = 'none') 
dev.off()


######### Replication analysis
common <- intersect(ICU_vs_nonICU_breda$OlinkID, ICU_vs_nonICU_radboud$OlinkID)
ICU_vs_nonICU_breda %<>% filter(OlinkID %in% common)
ICU_vs_nonICU_radboud %<>% filter(OlinkID %in% common)
ICU_vs_nonICU_radboud %<>% arrange(match(OlinkID, ICU_vs_nonICU_breda$OlinkID))

stopifnot(ICU_vs_nonICU_radboud$OlinkID == ICU_vs_nonICU_breda$OlinkID)


# How replicable?
df <- data.frame(
  'OlinkID'          = ICU_vs_nonICU_radboud$OlinkID,
  'protein'          = ICU_vs_nonICU_radboud$Assay, 
  'logFC.radboud'    = ICU_vs_nonICU_radboud$logFC, 
  'adj.pval.radboud' = ICU_vs_nonICU_radboud$adj.P.Val, 
  'sig.radboud'      = ICU_vs_nonICU_radboud$significance,
  'logFC.breda'      = ICU_vs_nonICU_breda$logFC, 
  'adj.pval.breda'   = ICU_vs_nonICU_breda$adj.P.Val, 
  'sig.breda'        = ICU_vs_nonICU_breda$significance
)

df$sig.both <- ifelse(df$sig.breda & df$sig.radboud, T, F)
table(df$sig.both)

# We want to show only the proteins that are significant in Breda cohort -> Discovery
table(df$sig.breda) 
# FALSE  TRUE 
# 32    30 
df %<>% filter(sig.breda == TRUE)

# 30 proteins significant in Breda
# Of these, we replicate 26.
# 28/30 same direction of regulation
table(df$sig.radboud)
table(df$logFC.radboud > 0, df$logFC.breda > 0)

fit <- lm(df$logFC.breda ~ df$logFC.radboud)
summary(fit)
label.cor <- paste("italic(R) ^ 2 ==",signif(summary(fit)$adj.r.squared, 2))
label.P <- paste('italic(p) == 5 %*% 10^-12')


pdf('output/replication_DE_breda_radboud.pdf', width = 4.5, height = 3.5)
ggplot(data = df) +
  geom_vline(xintercept = 0, lty=2, color = 'darkgrey') +
  geom_hline(yintercept = 0, lty=2, color = 'darkgrey') +
  geom_point(aes(x = logFC.breda, y = logFC.radboud, color = sig.radboud)) +
  labs(x = expression(paste(log[2], ' fold change cohort 1')), 
       y = expression(paste(log[2], ' fold change cohort 2')), 
       color = 'Significance in\ncohort 2') + 
  scale_color_manual(values = c('darkgrey', 'black')) +
  theme(legend.position = c(.8, .2)) +
  annotate(geom = 'text', label = label.cor, x = -0.5, y = 1, parse=T) +
  annotate(geom = 'text', label = label.P, x = -0.5, y = 0.8, parse=T)
dev.off()


# Heatmap of the proteins we are able to replicate
load('data/data.RData')
replicated <- df %>% 
  filter(sig.radboud == T) %>% 
  filter(sig.breda == T) %>% 
  pull(OlinkID)

# Radboud and Breda only use T1
annot.radboud %<>% filter(time == 'W1T1')
radboud <- radboud[which(rownames(radboud) %in% rownames(annot.radboud)), ]

annot.breda %<>% filter(timepoint == 'T1')
breda <- breda[which(rownames(breda) %in% rownames(annot.breda)), ]

# Combine the cohorts
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
stopifnot(all(rownames(annot) == rownames(df)))

df %<>% select(all_of(replicated))
colnames(df) <- conv %>% filter(OlinkID %in% colnames(df)) %>% arrange(match(OlinkID, colnames(df))) %>% pull(Assay)

annot %<>% arrange(condition)
df <- df[rownames(annot), ]
annot$cohort <- as.numeric(factor(annot$cohort, levels = c('Breda', 'Radboud')))
annot$cohort <- as.character(annot$cohort)

# Export in png because pdf makes the colors ugly on the annotation bar
png('output/heatmap_ICU_nonICU.png', width = 13, height = 4, units = 'in', res = 700)
pheatmap(t(df), 
         scale = 'row', cluster_cols = T,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100), 
         annotation_col = annot %>% select(cohort, condition), clustering_distance_cols = 'correlation', clustering_distance_rows = 'correlation',
         show_colnames = F, 
         breaks = seq(-1.5, 1.5, length.out = 100), 
         annotation_colors = list(condition = cols, 
                                  cohort = c('1' = 'chartreuse3', '2' = 'purple')), fontsize = 8)
dev.off()
