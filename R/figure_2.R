rm(list = ls())
try(dev.off())

library(ggplot2)
library(ggpubr)

setwd('/Users/martijnzoodsma/Documents/PhD/corona/Postcovid-targeted-proteomics/')
load('data/data.RData')

make.plot <- function(prot, nrow = NULL) {
  load('data/data.RData')
  std <- function(x) sd(x, na.rm = T)/sqrt(length(x))
  
  ##################
  # Radboud data
  radboud %<>% dplyr::select(all_of(prot))
  colnames(radboud) <- conv %>% filter(OlinkID %in% colnames(radboud)) %>% arrange(match(OlinkID, colnames(radboud))) %>% pull(Assay) %>% make.unique()
  
  radboud <- cbind(radboud, annot.radboud)
  radboud <- reshape2::melt(radboud, id.vars = colnames(annot.radboud))
  # radboud %<>% filter(time %in% c('W1T1', 'W1T2', 'W1T3', 'W2T1', 'W2T2', 'W2T3'))
  radboud %<>% filter(time %in% c('W1T1', 'W1T2', 'W1T3', 'W2T1', 'W2T2'))
  # radboud %<>% filter(time %in% c('W1T1', 'W1T2', 'W1T3'))
  
  radboud %<>% na.omit(radboud)
  
  radboud %<>% group_by(variable, time, condition) %>% summarise(mean = mean(value, na.rm = T), se = std(value))
  # radboud$time <- factor(radboud$time, levels = c('W1T1', 'W1T2', 'W1T3', 'W2T1', 'W2T2', 'W2T3', 'postcovid', 'healthy'))
  radboud$time <- factor(radboud$time, levels = c('W1T1', 'W1T2', 'W1T3', 'W2T1', 'W2T2', 'postcovid', 'healthy'))
  # radboud$time <- factor(radboud$time, levels = c('W1T1', 'W1T2', 'W1T3', 'postcovid', 'healthy'))
  
  
  
  ##################
  # MHH data
  mhh %<>% dplyr::select(all_of(prot))
  colnames(mhh) <- conv %>% filter(OlinkID %in% colnames(mhh)) %>% arrange(match(OlinkID, colnames(mhh))) %>% pull(Assay) %>% make.unique()
  mhh <- cbind(mhh, annot.mhh %>% dplyr::select(condition))
  
  mhh <- reshape2::melt(mhh, id.vars = 'condition')
  mhh %<>% filter(condition != 'bridge')
  mhh %<>% na.omit(mhh)
  # mhh$condition <- factor(mhh$condition, levels = c('W1T1', 'W1T2', 'W1T3', 'W2T1', 'W2T2', 'W2T3', 'postcovid', 'healthy'))
  mhh$condition <- factor(mhh$condition, levels = c('W1T1', 'W1T2', 'W1T3', 'W2T1', 'W2T2', 'postcovid', 'healthy'))
  # mhh$condition <- factor(mhh$condition, levels = c('W1T1', 'W1T2', 'W1T3', 'postcovid', 'healthy'))
  
  # Set the right colors.
  radboud$color <- ''
  radboud.1 <- radboud %>% filter(time %in% c('W1T1', 'W1T2', 'W1T3'))
  radboud.1[which(radboud.1$condition == 'ICU'), 'color']<- '#BC3C29FF'
  radboud.1[which(radboud.1$condition == 'non-ICU'), 'color'] <- '#E18727FF'
  
  radboud.2 <- radboud %>% filter(time %in% c('W2T1', 'W2T2'))
  radboud.2[which(radboud.2$condition == 'ICU'), 'color']<- 'gray'
  radboud.2[which(radboud.2$condition == 'non-ICU'), 'color'] <- 'gray'
  
  # Manually construct a dataframe for the connecting lines from W1T3 to W2T1 in gray
  x <- rbind(
    radboud.1 %>% filter(time =='W1T3'), 
    radboud.2 %>% filter(time == 'W2T1')
  ) %>% mutate(color = 'gray')
  x
  pd <- position_dodge(width = 0.2)
  
  pl <- ggplot(position = pd) +
    scale_color_manual(values = radboud.1$color, 
                       labels = radboud.1$condition) +
    geom_point(data = radboud.1, 
               aes(x = time, y = mean, color = color), position = pd) +
    geom_line(data = radboud.1, 
              aes(x = time, y = mean, color = color, group = condition), position = pd) +
    geom_errorbar(data = radboud.1, 
                  aes(x = time, y = mean, color = color, ymax = mean + se, ymin = mean - se), position = pd) +
    scale_x_discrete(drop = F) +
    geom_boxplot(data = mhh, aes(x = condition, y = value)) +
    geom_point(data = radboud.2, 
               aes(x = time, y = mean), color = radboud.2$color) +
    geom_line(data = radboud.2, 
              aes(x = time, y = mean, group = condition), color = radboud.2$color, lty = 2) +
    geom_line(data = x, 
              aes(x = time, y = mean, group = condition), color = x$color, lty = 2) +
    geom_errorbar(data = radboud.2, 
                  aes(x = time, y = mean, ymax = mean + se, ymin = mean - se), position = pd, color = radboud.2$color) +
    theme_classic() +
    labs(y = 'Protein level', 
         color = 'Condition') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
          axis.title.x = element_blank()) +
    facet_wrap(~ variable, scales = 'free_y', nrow = nrow)
  
  return(pl)
}

# Pie chart and representative example for each of the categories for ANOVA
pvals <- read.csv('output/pvalues_anova.csv', row.names=1, header=T)

n <- nrow(pvals)
n.sig.time <- pvals %>% filter(pval.time < 0.05) %>% nrow()
n.sig.condition <- pvals %>% filter(pval.condition < 0.05) %>% nrow()
n.sig.interaction <- pvals %>% filter(pval.interaction < 0.05) %>% nrow()

time <- data.frame(
  class = c('Significant', 'Non-significant'), 
  n = c(n, n), 
  prop = c(n.sig.time / n, (n - n.sig.time) / n)
) %>% mutate(label.pos = cumsum(prop) - 0.5 * prop)

time$label <- paste0((round(time$prop, 2) * 100), '%')
time$label <- paste0(time$class, ':\n', time$label)

pie.time <- ggplot(time, aes(x = '', y = prop, fill = class)) +
  geom_bar(stat = 'identity', color = 'white', width = 1) +
  coord_polar('y', start = 0) +
  theme_void() +
  labs(fill = '') +
  geom_text(aes(y = label.pos, label = label), color = 'black', size=3) + 
  scale_fill_manual(values = col <- c('lightgrey', "darkgreen")) +
  theme(legend.position =  'none', text = element_text(size  = 20))

condition <- data.frame(
  class = c('Significant', 'Non-significant'), 
  n = c(n, n), 
  prop = c(n.sig.condition / n, (n - n.sig.condition) / n)
) %>% mutate(label.pos = cumsum(prop) - 0.5 * prop)


condition$label <- paste0((round(condition$prop, 2) * 100), '%')
condition$label <- paste0(condition$class, ':\n', condition$label)

pie.condition <- ggplot(condition, aes(x = '', y = prop, fill = class)) +
  geom_bar(stat = 'identity', color = 'white', width = 1) +
  coord_polar('y', start = 0) +
  theme_void() +
  labs(fill = '') +
  geom_text(aes(y = label.pos, label = label), color = 'black', size=3) + 
  scale_fill_manual(values = col <- c('lightgrey', "darkgreen")) +
  theme(legend.position =  'none')


interaction <- data.frame(
  class = c('Significant', 'Non-significant'), 
  n = c(n, n), 
  prop = c(n.sig.interaction / n, (n - n.sig.interaction) / n)
) %>% mutate(label.pos = cumsum(prop) - 0.5 * prop)


interaction$label <- paste0((round(interaction$prop, 2) * 100), '%')
interaction$label <- paste0(interaction$class, ':\n', interaction$label)

pie.interaction <- ggplot(interaction, aes(x = '', y = prop, fill = class)) +
  geom_bar(stat = 'identity', color = 'white', width = 1) +
  coord_polar('y', start = 0) +
  theme_void() +
  labs(fill = '') +
  geom_text(aes(y = label.pos, label = label), color = 'black', size=3) + 
  scale_fill_manual(values = col <- c('lightgrey', "darkgreen")) +
  theme(legend.position =  'none')

# Plot the pies next to three examples
ex.time <- 'OID00535'
ex.condition <- 'OID00408'
ex.interaction <- 'OID01305'

ex.time %in% (pvals %>% filter(pval.time < 0.05) %>% pull(biomarker))
ex.condition %in% (pvals %>% filter(pval.condition < 0.05) %>% pull(biomarker))
ex.interaction %in% (pvals %>% filter(pval.interaction < 0.05) %>% pull(biomarker))

widths <- c(1, 1.5)

pdf('output/example_time.pdf', width = 8, height = 4)
ggarrange(pie.time, make.plot(ex.time), nrow = 1, widths = widths)
dev.off()

pdf('output/example_condition.pdf', width = 8, height = 4)
ggarrange(pie.condition, make.plot(ex.condition), nrow = 1, widths = widths)
dev.off()

pdf('output/example_interaction.pdf', width = 8, height = 4)
ggarrange(pie.interaction, make.plot(ex.interaction), nrow = 1, widths = widths)
dev.off()


# Manually picked proteins Fig. 2
# With per individual plot: 3
# Height: 3
prots <- c('TWEAK', 'AXIN1', 'ITGB1BP2', 'STAMBP')
prots <- conv %>% filter(Assay %in% prots) %>% pull(OlinkID)

pdf('output/heterogene_proteins_fig2.pdf', width = 12, height = 3)
make.plot(prots, nrow = 1)
dev.off()


prots <- c('BMP-6', 'EFEMP1', 'MMP-1')
prots <- conv %>% filter(Assay %in% prots) %>% pull(OlinkID)

pdf('output/fibrotic_proteins_fig2.pdf', width = 9, height = 3)
make.plot(prots, nrow = 1)
dev.off()


prots <- c('TWEAK', 'TRAIL', 'CD40-L', 'TNF')
prots <- conv %>% filter(Assay %in% prots) %>% pull(OlinkID)

pdf('output/apoptotic_proteins_fig2.pdf', width = 12, height = 3)
make.plot(prots, nrow = 1)
dev.off()
