# Differential expression between:
# 1. Nijmegen ICU vs Hannover healthy
# 2. Nijmegen nonICU vs Hannover healthy
# 3. Hannover convalescent vs Hannover healthy

# Visualization: volcano plots, Venn diagrams of the overlap between conditions
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

annot <- rbind(
  annot.radboud %>% filter(time == 'W1T1') %>% select(age, gender, condition, cohort), 
  annot.mhh %>% select(age, gender, condition, cohort)
)

df <- rbind(
  mhh %>% select(intersect(colnames(mhh), colnames(radboud))), 
  radboud %>% select(intersect(colnames(mhh), colnames(radboud)))
) # We are testing 220 proteins in total

# Only T1
df %<>% filter(rownames(df) %in% intersect(rownames(df), rownames(annot)))
annot %<>% filter(rownames(annot) %in% intersect(rownames(annot), rownames(df))) %>% arrange(match(rownames(annot), rownames(df)))



# Differential expression analyses
# 1. Nijmegen ICU vs Hannover healthy
sub.annot <- annot %>% filter(condition %in% c('healthy', 'ICU')) %>% na.omit()
sub.df <- df %>% filter(row.names(df) %in% rownames(sub.annot))

stopifnot(all(rownames(sub.annot) == rownames(sub.df)))

condition <- sub.annot$condition
age <- sub.annot$age
gender <- sub.annot$gender

design <- model.matrix(~0 + condition + age + gender)
fit <- lmFit(as.matrix(t(sub.df)), design = design)
cont <- makeContrasts(cond = 'conditionICU-conditionhealthy', levels = design)
fit <- contrasts.fit(fit, cont)
fit <- eBayes(fit)

ICU_vs_healthy <- topTable(fit = fit, adjust.method = "BH", number = Inf)
ICU_vs_healthy$OlinkID <- rownames(ICU_vs_healthy)
ICU_vs_healthy$significance <- ifelse(ICU_vs_healthy$adj.P.Val < 0.05, TRUE, FALSE)
ICU_vs_healthy$direction <- ifelse(ICU_vs_healthy$logFC < 0.0, 'Downregulated', 'Upregulated')
ICU_vs_healthy <- cbind(ICU_vs_healthy, 
                        conv %>% 
                          filter(OlinkID %in% ICU_vs_healthy$OlinkID) %>% 
                          arrange(match(OlinkID, ICU_vs_healthy$OlinkID)) %>% 
                          select(Assay, Uniprot.ID, Olink.panel))

# 2. Nijmegen nonICU vs Hannover healthy
sub.annot <- annot %>% filter(condition %in% c('healthy', 'non-ICU')) %>% na.omit()
sub.df <- df %>% filter(row.names(df) %in% rownames(sub.annot))

stopifnot(all(rownames(sub.annot) == rownames(sub.df)))

sub.annot$condition[which(sub.annot$condition == 'non-ICU')] <- 'nonICU'

condition <- sub.annot$condition
age <- sub.annot$age
gender <- sub.annot$gender

design <- model.matrix(~0 + condition + age + gender)
fit <- lmFit(as.matrix(t(sub.df)), design = design)
cont <- makeContrasts(cond = 'conditionnonICU-conditionhealthy', levels = design)
fit <- contrasts.fit(fit, cont)
fit <- eBayes(fit)

nonICU_vs_healthy <- topTable(fit = fit, adjust.method = "BH", number = Inf)
nonICU_vs_healthy$OlinkID <- rownames(nonICU_vs_healthy)
nonICU_vs_healthy$significance <- ifelse(nonICU_vs_healthy$adj.P.Val < 0.05, TRUE, FALSE)
nonICU_vs_healthy$direction <- ifelse(nonICU_vs_healthy$logFC < 0.0, 'Downregulated', 'Upregulated')
nonICU_vs_healthy <- cbind(nonICU_vs_healthy, 
                           conv %>% 
                             filter(OlinkID %in% nonICU_vs_healthy$OlinkID) %>% 
                             arrange(match(OlinkID, nonICU_vs_healthy$OlinkID)) %>% 
                             select(Assay, Uniprot.ID, Olink.panel))


# 3. Hannover convalescent vs Hannover healthy
# Different data because we want to include the neurology panel as well
annot.mhh %<>% filter(condition %in% c('healthy', 'postcovid')) 
annot.mhh$condition <- factor(annot.mhh$condition, levels = c('healthy', 'postcovid'))

mhh %<>% filter(rownames(mhh) %in% rownames(annot.mhh))
stopifnot(all(rownames(mhh) == rownames(annot.mhh)))

condition <- annot.mhh$condition
age <- annot.mhh$age
gender <- annot.mhh$gender

design <- model.matrix(~0 + condition + age + gender)
fit <- lmFit(as.matrix(t(mhh)), design = design)
cont <- makeContrasts(cond = 'conditionpostcovid-conditionhealthy', levels = design)
fit <- contrasts.fit(fit, cont)
fit <- eBayes(fit)

postcovid_vs_healthy <- topTable(fit = fit, adjust.method = "BH", number = Inf)
postcovid_vs_healthy$OlinkID <- rownames(postcovid_vs_healthy)
postcovid_vs_healthy$significance <- ifelse(postcovid_vs_healthy$adj.P.Val < 0.05, TRUE, FALSE)
postcovid_vs_healthy$direction <- ifelse(postcovid_vs_healthy$logFC < 0.0, 'Downregulated', 'Upregulated')
postcovid_vs_healthy <- cbind(postcovid_vs_healthy, 
                              conv %>% 
                                filter(OlinkID %in% postcovid_vs_healthy$OlinkID) %>% 
                                arrange(match(OlinkID, postcovid_vs_healthy$OlinkID)) %>% 
                                select(Assay, Uniprot.ID, Olink.panel))


##### Export results
res <- list(ICU_vs_healthy = ICU_vs_healthy, 
            nonICU_vs_healthy = nonICU_vs_healthy, 
            postcovid_vs_healthy = postcovid_vs_healthy)

write.xlsx(x = res, file = 'output/DE.xlsx')

##### Visualize

# Volcano plots for all comparisons
theme_set(theme_classic() + theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 12)))
n = 13
show = 25

pdf('output/volcano_DE_radboud.pdf', width = 15, height = 4, onefile = F)

pl <- ggarrange(
  ggplot() +
    geom_point(data = ICU_vs_healthy %>% filter(significance == F) , aes(x = logFC, y = -log10(adj.P.Val)), color = 'lightgray', show.legend = F) +
    geom_point(data = ICU_vs_healthy %>% filter(significance == T) , aes(x = logFC, y = -log10(adj.P.Val), color = direction)) +
    labs(x = expression(paste(log[2], ' fold change')), y = expression(paste(-log[10], '(adj. p-value)')), color = '', title = 'ICU vs healthy') +
    scale_color_manual(values = c('#0072B5FF', '#BC3C29FF')) +
    geom_hline(yintercept = -log10(0.05), lty = 'dotted') +
    geom_vline(xintercept = 0.0, lty = 'dotted') + 
    xlim(c(-5, 6)) +
    ylim(c(0, 90)) +
    ggrepel::geom_label_repel(data = head(x = ICU_vs_healthy, show), aes(label = Assay, x = logFC, y = -log10(adj.P.Val)), 
                              box.padding = 0.5, show.legend = F,max.overlaps = n),
  
  ggplot() +
    geom_point(data = nonICU_vs_healthy %>% filter(significance == F) , aes(x = logFC, y = -log10(adj.P.Val)), color = 'lightgray', show.legend = F) +
    geom_point(data = nonICU_vs_healthy %>% filter(significance == T) , aes(x = logFC, y = -log10(adj.P.Val), color = direction)) +
    labs(x = expression(paste(log[2], ' fold change')), y = expression(paste(-log[10], '(adj. p-value)')), color = '', title = 'non-ICU vs healthy') +
    scale_color_manual(values = c('#0072B5FF', '#BC3C29FF')) +
    geom_hline(yintercept = -log10(0.05), lty = 'dotted') +
    geom_vline(xintercept = 0.0, lty = 'dotted') +
    xlim(c(-5, 6)) +
    ylim(c(0, 90)) +
    ggrepel::geom_label_repel(data = head(x = nonICU_vs_healthy, show), aes(label = Assay, x = logFC, y = -log10(adj.P.Val)), 
                              box.padding = 0.5, show.legend = F, max.overlaps = n), 
  
  ggplot() +
    geom_point(data = postcovid_vs_healthy %>% filter(significance == F) , aes(x = logFC, y = -log10(adj.P.Val)), color = 'lightgray', show.legend = F) +
    geom_point(data = postcovid_vs_healthy %>% filter(significance == T) , aes(x = logFC, y = -log10(adj.P.Val), color = direction)) +
    labs(x = expression(paste(log[2], ' fold change')), y = expression(paste(-log[10], '(adj. p-value)')), color = '', title = 'Convalescent vs healthy') +
    scale_color_manual(values = c('#0072B5FF', '#BC3C29FF')) +
    geom_hline(yintercept = -log10(0.05), lty = 'dotted') +
    geom_vline(xintercept = 0.0, lty = 'dotted') +
    xlim(c(-5, 6)) +
    ylim(c(0, 90)) +
    ggrepel::geom_label_repel(data = head(x = postcovid_vs_healthy, show), aes(label = Assay, x = logFC, y = -log10(adj.P.Val)), 
                              box.padding = 0.5, show.legend = F, max.overlaps = 35),
  
  ncol = 3, nrow = 1, common.legend = T, legend =  'none'
)
pl
dev.off()

# High def for poster
jpeg(filename = 'volcanos_conv_highres.jpeg', width = 15, height = 4, units = 'in', res = 1200)
pl
dev.off()


#### What are the biomarkers for each of the disease conditions?
# What are the biomarkers for each of the conditions?
# For postcovid, extract the Neurology panel and write this separately
# For postcovid, filter to take only proteins that are included in the ICU and nonICU proteins.
biomarkers <- list(
  'ICU' = ICU_vs_healthy %>% filter(significance == TRUE) %>% pull(OlinkID), 
  'nonICU' = nonICU_vs_healthy %>% filter(significance ==TRUE) %>% pull(OlinkID), 
  'postcovid' = postcovid_vs_healthy %>% filter(significance ==TRUE) %>% filter(Olink.panel != 'Olink NEUROLOGY') %>% pull(OlinkID), 
  'postcovid neurology' = postcovid_vs_healthy %>% filter(significance == TRUE) %>% filter(Olink.panel == 'Olink NEUROLOGY') %>% pull(OlinkID)
)

write.xlsx('output/biomarkers.xlsx', x = biomarkers)


length(biomarkers$ICU)
length(biomarkers$nonICU)
length(biomarkers$postcovid)
length(biomarkers$`postcovid neurology`)

# Venn diagram separated for up and downregulated biomarkers
ICU.up <- ICU_vs_healthy %>% filter(significance == TRUE) %>% filter(direction == 'Upregulated') %>% pull(OlinkID)
nonICU.up <- nonICU_vs_healthy %>% filter(significance == TRUE) %>% filter(direction == 'Upregulated') %>% pull(OlinkID)
postcovid.up <- postcovid_vs_healthy %>% filter(significance == TRUE) %>% filter(Olink.panel != 'Olink NEUROLOGY') %>% filter(direction == 'Upregulated') %>% pull(OlinkID)

ICU.down <- ICU_vs_healthy %>% filter(significance == TRUE) %>% filter(direction == 'Downregulated') %>% pull(OlinkID)
nonICU.down <- nonICU_vs_healthy %>% filter(significance == TRUE) %>% filter(direction == 'Downregulated') %>% pull(OlinkID)
postcovid.down <- postcovid_vs_healthy %>% filter(significance == TRUE) %>% filter(Olink.panel != 'Olink NEUROLOGY') %>% filter(direction == 'Downregulated') %>% pull(OlinkID)

res.up <- list('ICU' = ICU.up, 'non-ICU' = nonICU.up, 'Convalescent' = postcovid.up)
res.down <- list('ICU' = ICU.down, 'non-ICU' = nonICU.down, 'Convalescent' = postcovid.down)

pdf('output/venn_biomarkers_up.pdf', width = 4, height = 4)
ggVennDiagram(x = res.up, label_alpha = 0, label = 'count') +
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
              label_alpha = 0) +
  labs(title = 'Downregulated biomarkers') +
  scale_fill_gradient2(low = 'blue', high = 'red', limits = c(0, 70)) +
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        legend.position = 'none') 
dev.off()
