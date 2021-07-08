# Differential expression between:
  #1. MHH healthy vs Radboud ICU
  #2. MHH healthy vs Radboud non-ICU
  #3. MHH healthy vs MHH postcovid
# using limma corrected for age and gender

# Visualization: volcano plots, Venn diagrams and heatmaps

library(openxlsx)
library(magrittr)
library(limma)
library(dplyr)
library(ggpubr)
library(ggVennDiagram)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

rm(list = ls())
load('data/data.RData')

# prepare data. Take only common proteins and the first timepoint for Nijmegen cohort
# Prepare the data. We take only the common proteins and patients who were in the ICU at their first timepoint for the Nijmegen cohort

annot <- rbind(
  annot.radboud %>% filter(time == 'W1T1') %>% select(age, gender, condition, cohort), 
  annot.mhh %>% select(age, gender, condition, cohort)
)

df <- rbind(
  mhh %>% select(intersect(colnames(mhh), colnames(radboud))), 
  radboud %>% select(intersect(colnames(mhh), colnames(radboud)))
)

# Only T1
df %<>% filter(rownames(df) %in% intersect(rownames(df), rownames(annot)))
annot %<>% filter(rownames(annot) %in% intersect(rownames(annot), rownames(df))) %>% arrange(match(rownames(annot), rownames(df)))


##### MHH healthy vs Radboud ICU
sub.annot <- annot %>% filter(condition %in% c('healthy', 'ICU'))
sub.df <- df %>% filter(row.names(df) %in% rownames(sub.annot))

sum(is.na(sub.annot))
sum(is.na(sub.df))

sub.annot$condition <- factor(sub.annot$condition, levels = c('healthy', 'ICU'))

str(sub.annot)

# Limma model
design <- model.matrix(~ sub.annot$condition + sub.annot$age + sub.annot$gender)
head(design)
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

##### MHH healthy vs Radboud non-ICU
sub.annot <- annot %>% filter(condition %in% c('healthy', 'non-ICU'))
sub.df <- df %>% filter(row.names(df) %in% rownames(sub.annot))

sum(is.na(sub.annot))
sum(is.na(sub.df))

sub.annot$condition <- factor(sub.annot$condition, levels = c('healthy', 'non-ICU'))

# Limma model
design <- model.matrix(~ sub.annot$condition + sub.annot$age + sub.annot$gender)
head(design)
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

##### Radboud nonICU vs Radboud ICU
sub.annot <- annot %>% filter(condition %in% c('non-ICU', 'ICU'))
sub.df <- df %>% filter(row.names(df) %in% rownames(sub.annot))

sum(is.na(sub.annot))
sum(is.na(sub.df))

sub.annot$condition <- factor(sub.annot$condition, levels = c('non-ICU', 'ICU'))
table(sub.annot$condition)
str(sub.annot)

# Limma model
design <- model.matrix(~ sub.annot$condition + sub.annot$age + sub.annot$gender)
head(design)
colnames(design) <- c('Intercept', 'condition', 'age', 'gender_male')
fit <- lmFit(as.matrix(t(sub.df)), design = design)
fit <- eBayes(fit)
nonICU_vs_ICU <- topTable(fit = fit, adjust.method = "BH", number = Inf, coef = 'condition')
nonICU_vs_ICU$OlinkID <- rownames(nonICU_vs_ICU)
nonICU_vs_ICU$significance <- ifelse(nonICU_vs_ICU$adj.P.Val < 0.05, TRUE, FALSE)
nonICU_vs_ICU$direction <- ifelse(nonICU_vs_ICU$logFC < 0.0, 'Downregulated', 'Upregulated')
nonICU_vs_ICU <- cbind(nonICU_vs_ICU, 
                       conv %>% 
                         filter(OlinkID %in% nonICU_vs_ICU$OlinkID) %>% 
                         arrange(match(OlinkID, nonICU_vs_ICU$OlinkID)) %>% 
                         select(Assay, Uniprot.ID, Olink.panel))


##### MHH healthy vs MHH postcovid
# Different data because we want to include the neurology panel as well
annot.mhh %<>% filter(condition %in% c('healthy', 'postcovid'))
annot.mhh$condition <- factor(annot.mhh$condition, levels = c('healthy', 'postcovid'))

mhh %<>% filter(rownames(mhh) %in% rownames(annot.mhh))

# Limma model
design <- model.matrix(~ annot.mhh$condition + annot.mhh$age + annot.mhh$gender)
head(design)
colnames(design) <- c('Intercept', 'condition', 'age', 'gender_male')


fit <- lmFit(as.matrix(t(mhh)), design = design)
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



##### Export results
res <- list(healthy_vs_ICU = healthy_vs_ICU, 
            healthy_vs_nonICU = healthy_vs_nonICU, 
            healthy_vs_postcovid = healthy_vs_postcovid, 
            nonICU_vs_ICU = nonICU_vs_ICU)

write.xlsx(x = res, file = 'output/DE.xlsx')

##### Visualize

# Volcano plots for all comparisons
theme_set(theme_classic() + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 12)))
n = 13

pdf('output/volcano_DE_radboud.pdf', width = 15, height = 4, onefile = F)

ggarrange(
  ggplot() +
    geom_point(data = healthy_vs_ICU %>% filter(significance == F) , aes(x = logFC, y = -log10(adj.P.Val)), color = 'lightgray', show.legend = F) +
    geom_point(data = healthy_vs_ICU %>% filter(significance == T) , aes(x = logFC, y = -log10(adj.P.Val), color = direction)) +
    labs(x = expression(paste(log[2], ' fold change')), y = expression(paste(-log[10], '(adj. p-value)')), color = '', title = 'ICU vs healthy') +
    scale_color_manual(values = c('#0072B5FF', '#BC3C29FF')) +
    geom_hline(yintercept = -log10(0.05), lty = 'dotted') +
    geom_vline(xintercept = 0.0, lty = 'dotted') + 
    xlim(c(-5, 6)) +
    ylim(c(0, 90)) +
    ggrepel::geom_label_repel(data = head(x = healthy_vs_ICU, 25), aes(label = Assay, x = logFC, y = -log10(adj.P.Val)), box.padding = 0.5, show.legend = F,max.overlaps = n),
  
  ggplot() +
    geom_point(data = healthy_vs_nonICU %>% filter(significance == F) , aes(x = logFC, y = -log10(adj.P.Val)), color = 'lightgray', show.legend = F) +
    geom_point(data = healthy_vs_nonICU %>% filter(significance == T) , aes(x = logFC, y = -log10(adj.P.Val), color = direction)) +
    labs(x = expression(paste(log[2], ' fold change')), y = expression(paste(-log[10], '(adj. p-value)')), color = '', title = 'non-ICU vs healthy') +
    scale_color_manual(values = c('#0072B5FF', '#BC3C29FF')) +
    geom_hline(yintercept = -log10(0.05), lty = 'dotted') +
    geom_vline(xintercept = 0.0, lty = 'dotted') +
    xlim(c(-5, 6)) +
    ylim(c(0, 90)) +
    ggrepel::geom_label_repel(data = head(x = healthy_vs_nonICU, 25), aes(label = Assay, x = logFC, y = -log10(adj.P.Val)), box.padding = 0.5, show.legend = F, max.overlaps = n), 
  
  ggplot() +
    geom_point(data = healthy_vs_postcovid %>% filter(significance == F) , aes(x = logFC, y = -log10(adj.P.Val)), color = 'lightgray', show.legend = F) +
    geom_point(data = healthy_vs_postcovid %>% filter(significance == T) , aes(x = logFC, y = -log10(adj.P.Val), color = direction)) +
    labs(x = expression(paste(log[2], ' fold change')), y = expression(paste(-log[10], '(adj. p-value)')), color = '', title = 'post-COVID-19 vs healthy') +
    scale_color_manual(values = c('#0072B5FF', '#BC3C29FF')) +
    geom_hline(yintercept = -log10(0.05), lty = 'dotted') +
    geom_vline(xintercept = 0.0, lty = 'dotted') +
    xlim(c(-5, 6)) +
    ylim(c(0, 90)) +
    ggrepel::geom_label_repel(data = head(x = healthy_vs_postcovid, 5), aes(label = Assay, x = logFC, y = -log10(adj.P.Val)), box.padding = 0.5, show.legend = F, max.overlaps = n),
  
  ncol = 3, nrow = 1, common.legend = T, legend =  'right'
)
dev.off()


#### What are the biomarkers for each of the disease conditions?
  # e.g. extract sigDE proteins for each of the conditions and export
# For postcovid, extract the Neurology panel and write this separately.

# What are the biomarkers for each of the conditions?
# For postcovid, extract the Neurology panel and write this separately
# For postcovid, filter to take only proteins that are included in the ICU and nonICU proteins.

healthy_vs_postcovid_filtered <- healthy_vs_postcovid[which(healthy_vs_postcovid$OlinkID %in% healthy_vs_ICU$OlinkID), ]

biomarkers <- list(
  'ICU' = healthy_vs_ICU %>% filter(significance == TRUE) %>% pull(OlinkID), 
  'nonICU' = healthy_vs_nonICU %>% filter(significance ==TRUE) %>% pull(OlinkID), 
  'postcovid' = healthy_vs_postcovid_filtered %>% filter(significance ==TRUE) %>% filter(Olink.panel != 'Olink NEUROLOGY') %>% pull(OlinkID), 
  'postcovid neurology' = healthy_vs_postcovid %>% filter(significance == TRUE) %>% filter(Olink.panel == 'Olink NEUROLOGY') %>% pull(OlinkID)
)

write.xlsx('output/biomarkers.xlsx', x = biomarkers)


length(biomarkers$ICU)
length(biomarkers$nonICU)
length(biomarkers$postcovid)
length(biomarkers$`postcovid neurology`)

# Venn diagram separated for up and downregulated biomarkers
ICU.up <- healthy_vs_ICU %>% filter(significance == TRUE) %>% filter(direction == 'Upregulated') %>% pull(OlinkID)
nonICU.up <- healthy_vs_nonICU %>% filter(significance == TRUE) %>% filter(direction == 'Upregulated') %>% pull(OlinkID)
postcovid.up <- healthy_vs_postcovid_filtered %>% filter(significance == TRUE) %>% filter(Olink.panel != 'Olink NEUROLOGY') %>% filter(direction == 'Upregulated') %>% pull(OlinkID)

ICU.down <- healthy_vs_ICU %>% filter(significance == TRUE) %>% filter(direction == 'Downregulated') %>% pull(OlinkID)
nonICU.down <- healthy_vs_nonICU %>% filter(significance == TRUE) %>% filter(direction == 'Downregulated') %>% pull(OlinkID)
postcovid.down <- healthy_vs_postcovid_filtered %>% filter(significance == TRUE) %>% filter(Olink.panel != 'Olink NEUROLOGY') %>% filter(direction == 'Downregulated') %>% pull(OlinkID)

res.up <- list('ICU' = ICU.up, 'non-ICU' = nonICU.up, 'postcovid' = postcovid.up)
res.down <- list('ICU' = ICU.down, 'non-ICU' = nonICU.down, 'postcovid' = postcovid.down)

# devtools::install_github('MZoodsma/ggVennDiagram')
pdf('output/venn_biomarkers_up.pdf', width = 4, height = 4)
ggVennDiagram(x = res.up, 
              label = 'count', 
              label_alpha = 0, 
              category = data.frame(x = c(2, -2.5, 7),
                                    y = c(8.5, -4, -4), 
                                    label = names(res.up))) +
  labs(title = 'Upregulated biomarkers') +
  scale_fill_gradient2(low = 'blue', high = 'red', limits = c(0, 70)) +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        legend.position = 'none') 
dev.off()

pdf('output/venn_biomarkers_down.pdf', width = 4, height = 4)
ggVennDiagram(x = res.down, 
              label = 'count', 
              label_alpha = 0, 
              category = data.frame(x = c(2, -2.5, 7),
                                    y = c(8.5, -4, -4), 
                                    label = names(res.up))) +
  labs(title = 'Downregulated biomarkers') +
  scale_fill_gradient2(low = 'blue', high = 'red', limits = c(0, 70)) +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        legend.position = 'none') 
dev.off()

# Heatmaps for the different upregulated and downregulated proteins
load('data/data.RData') 

make.heatmap <- function(proteins, color, lim, ...) {
  
  df <- rbind(mhh %>% select(all_of(proteins)), radboud %>% select(all_of(proteins)))
  annot <- rbind(annot.mhh %>% select(age, gender, condition, cohort), annot.radboud %>% select(age, gender, condition, cohort))
  
  annot %<>% filter(condition != 'bridge')
  df <- df[which(rownames(df) %in% rownames(annot)), ]
  
  df <- na.omit(df)
  annot <- annot[which(rownames(annot) %in% rownames(df)), ]
  annot$condition <- factor(annot$condition, levels = c('ICU', 'non-ICU', 'postcovid', 'healthy'))
  # 
  annot <- annot[order(annot$condition), ]

  df %<>% arrange(match(rownames(df), rownames(annot)))
  
  colnames(df) <- conv %>% filter(OlinkID %in% colnames(df)) %>% arrange(match(OlinkID, colnames(df))) %>% pull(Assay)
  cols <- list(condition = c(ICU = 'darkred', 'non-ICU' = 'orange', postcovid = 'yellow', healthy = 'green'))
  
  stopifnot(all(rownames(df) == rownames(annot)))
  
  ph <- pheatmap( t(df), 
                  scale = 'row', 
                  show_colnames = F, 
                  cluster_cols = F,
                  annotation_col = annot %>% select(condition),
                  annotation_colors = cols, 
                  color = color, 
                  breaks = seq(-lim, lim, length.out = 100), ..., )
  return(ph)
}

cols <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)

# Shared signature
# 27 shared upregulated and 6 shared downregulated
all_up <- Reduce(intersect, res.up)
all_down <- Reduce(intersect, res.down)
prots <- c(all_up, all_down)

pdf('output/heatmap_shared_biomarkers.pdf', width = 12)
make.heatmap(prots, color = cols, lim=1.5, cutree_rows = 2)
dev.off()

# Postcovid signature
# 25 specific downregulated, 8 specific upregulated
poscovid.signature <- c(setdiff(res.down$postcovid, c(res.down$ICU, res.down$`non-ICU`)), 
                        setdiff(res.up$postcovid, c(res.up$ICU, res.up$`non-ICU`)))
pdf('output/heatmap_postcovid_specific.pdf', width = 12)
make.heatmap(proteins = poscovid.signature, color = cols, lim = 2, cutree_rows = 2)
dev.off()

# 70 proteins that are recovered in postcovid.
recover.postcovid <- setdiff(intersect(res.up$ICU, res.up$`non-ICU`), res.up$postcovid)
pdf('output/heatmap_postcovid_recovered.pdf', width = 20, height = 12)
make.heatmap(recover.postcovid, color = cols, lim = 1.5)
dev.off()