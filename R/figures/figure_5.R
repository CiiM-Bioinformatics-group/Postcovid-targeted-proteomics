rm(list = ls())
try(dev.off())

library(pheatmap)
library(reshape2)
library(dplyr)
library(magrittr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

###################
### Load data and setup
###################

load('data/data.RData')


# Contruct entire df and annot
annot <- rbind(
  annot.mhh %>% select(age, gender, condition, cohort), 
  annot.radboud %>% select(age, gender, condition, cohort)
)

df <- rbind(
  mhh %>% select(intersect(colnames(mhh), colnames(radboud))), 
  radboud %>% select(intersect(colnames(mhh), colnames(radboud)))
)

df <- na.omit(df)
annot <- na.omit(annot)

df %<>% filter(rownames(df) %in% intersect(rownames(df), rownames(annot)))
annot %<>% filter(rownames(annot) %in% intersect(rownames(annot), rownames(df))) 
annot %<>% arrange(match(rownames(annot), rownames(df)))

all(rownames(annot) == rownames(df))

#####################
# Split in four conditions and determine clustering order for the four conditions
#####################

annot.ICU <- annot %>% filter(condition == 'ICU')
annot.nonICU <- annot %>% filter(condition == 'non-ICU')
annot.postcovid <- annot %>% filter(condition == 'postcovid')
annot.healthy <- annot %>% filter(condition == 'healthy')

df.ICU <- df[which(rownames(df) %in% rownames(annot.ICU)), ]
df.nonICU <- df[which(rownames(df) %in% rownames(annot.nonICU)), ]
df.postcovid <- df[which(rownames(df) %in% rownames(annot.postcovid)), ]
df.healthy <- df[which(rownames(df) %in% rownames(annot.healthy)), ]

cor.df.ICU <- cor(df.ICU, method = 'pearson')
cor.df.nonICU <- cor(df.nonICU, method = 'pearson')
cor.df.postcovid <- cor(df.postcovid, method = 'pearson')
cor.df.healthy <- cor(df.healthy, method = 'pearson')

ord.ICU <- pheatmap(cor.df.ICU)$tree_row$order
ord.nonICU <- pheatmap(cor.df.nonICU)$tree_row$order
ord.postcovid <- pheatmap(cor.df.postcovid)$tree_row$order
ord.healthy <- pheatmap(cor.df.healthy)$tree_row$order

dev.off()

do.heatmap <- function(cor.df.ICU, cor.df.nonICU, cor.df.postcovid, cor.df.healthy, ord) {
  
  cor.df.ICU <- cor.df.ICU[ord, ord]
  cor.df.nonICU <- cor.df.nonICU[ord, ord]
  cor.df.postcovid <- cor.df.postcovid[ord, ord]
  cor.df.healthy <- cor.df.healthy[ord, ord]
  
  # Set lower triangles including diagonal to NA and remove these rows
  cor.df.ICU[upper.tri(cor.df.ICU, diag = T)] <- NA
  cor.df.ICU.melt <- na.omit(melt(cor.df.ICU))
  
  cor.df.nonICU[upper.tri(cor.df.nonICU, diag = T)] <- NA
  cor.df.nonICU.melt <- na.omit(melt(cor.df.nonICU))
  
  cor.df.postcovid[upper.tri(cor.df.postcovid, diag = T)] <- NA
  cor.df.postcovid.melt <- na.omit(melt(cor.df.postcovid))
  
  cor.df.healthy[upper.tri(cor.df.healthy, diag = T)] <- NA
  cor.df.healthy.melt <- na.omit(melt(cor.df.healthy))
  
  
  stopifnot(all(cor.df.ICU.melt$Var1 == cor.df.nonICU.melt$Var1))
  stopifnot(all(cor.df.ICU.melt$Var2 == cor.df.nonICU.melt$Var2))
  
  stopifnot(all(cor.df.ICU.melt$Var1 == cor.df.postcovid.melt$Var1))
  stopifnot(all(cor.df.ICU.melt$Var2 == cor.df.postcovid.melt$Var2))
  
  stopifnot(all(cor.df.ICU.melt$Var1 == cor.df.healthy.melt$Var1))
  stopifnot(all(cor.df.ICU.melt$Var2 == cor.df.healthy.melt$Var2))
  
  lims = c(-0.5, 0.5)
  colours = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
  theme_set(theme_classic() +
              theme(axis.text.x = element_blank(), 
                    axis.ticks = element_blank(),
                    axis.text.y = element_blank(), 
                    plot.title = element_text(hjust = 0.5)))
  
  p <- ggarrange(
    ggplot() +
      geom_tile(data = cor.df.ICU.melt, aes(Var2, Var1, fill = value)) +
      labs(x = 'Proteins', y = 'Proteins', title = 'COVID-19 ICU') +
      scale_fill_gradientn(colours = colours, limits = lims, oob = scales::squish),
    
    ggplot() +
      geom_tile(data = cor.df.nonICU.melt, aes(Var2, Var1, fill = value)) +
      labs(x = 'Proteins', y = 'Proteins', title = 'COVID-19 non-ICU') +
      scale_fill_gradientn(colours = colours, limits = lims, oob = scales::squish),
    
    ggplot() +
      geom_tile(data = cor.df.postcovid.melt, aes(Var2, Var1, fill = value)) +
      labs(x = 'Proteins', y = 'Proteins', title = 'Postcovid') +
      scale_fill_gradientn(colours = colours, limits = lims, oob = scales::squish),
    
    ggplot() +
      geom_tile(data = cor.df.healthy.melt, aes(Var2, Var1, fill = value)) +
      labs(x = 'Proteins', y = 'Proteins', title = 'Healthy') +
      scale_fill_gradientn(colours = colours, limits = lims, oob = scales::squish),
    nrow = 1, ncol = 4, common.legend = T, legend = 'right'
  )
  return(p)
}

pdf('output/corr_heatmap_ICU.pdf', width = 15, height = 3, onefile = F)
do.heatmap(cor.df.ICU = cor.df.ICU, cor.df.nonICU = cor.df.nonICU, cor.df.postcovid = cor.df.postcovid, cor.df.healthy = cor.df.healthy, ord = ord.ICU)
dev.off()

pdf('output/corr_heatmap_nonICU.pdf', width = 15, height = 3, onefile = F)
do.heatmap(cor.df.ICU = cor.df.ICU, cor.df.nonICU = cor.df.nonICU, cor.df.postcovid = cor.df.postcovid, cor.df.healthy = cor.df.healthy, ord = ord.nonICU)
dev.off()

pdf('output/corr_heatmap_postcovid.pdf', width = 15, height = 3, onefile = F)
do.heatmap(cor.df.ICU = cor.df.ICU, cor.df.nonICU = cor.df.nonICU, cor.df.postcovid = cor.df.postcovid, cor.df.healthy = cor.df.healthy, ord = ord.postcovid)
dev.off()

pdf('output/corr_heatmap_healthy.pdf', width = 15, height = 3, onefile = F)
do.heatmap(cor.df.ICU = cor.df.ICU, cor.df.nonICU = cor.df.nonICU, cor.df.postcovid = cor.df.postcovid, cor.df.healthy = cor.df.healthy, ord = ord.healthy)
dev.off()





