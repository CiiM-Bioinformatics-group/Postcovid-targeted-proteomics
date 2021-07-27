rm(list = ls())
dev.off()

library(pheatmap)
library(reshape2)
library(dplyr)
library(magrittr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(Hmisc)
library(igraph)
library(RCy3)

###################
### Load data and setup
###################

setwd('/Users/martijnzoodsma/Documents/PhD/corona/Postcovid-targeted-proteomics/')
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

annot$condition[which(annot$condition == 'postcovid')] <- 'Convalescent'

annot.ICU <- annot %>% filter(condition == 'ICU')
annot.nonICU <- annot %>% filter(condition == 'non-ICU')
annot.postcovid <- annot %>% filter(condition == 'Convalescent')
annot.healthy <- annot %>% filter(condition == 'healthy')

df.ICU <- df[which(rownames(df) %in% rownames(annot.ICU)), ]
df.nonICU <- df[which(rownames(df) %in% rownames(annot.nonICU)), ]
df.postcovid <- df[which(rownames(df) %in% rownames(annot.postcovid)), ]
df.healthy <- df[which(rownames(df) %in% rownames(annot.healthy)), ]

####### ipgrah function
do.plot <- function(cormatrix) {
  
  dups <- c('CCL3', 'CXCL1', 'FGF-21', 'IL18', 'SCF')
  
  rem <- conv %>% filter(Assay %in% dups) %>% filter(Olink.panel == 'Olink INFLAMMATION') %>% pull(OlinkID)
  
  cormatrix %<>% select(-all_of(rem)) 
  
  test <- rcorr(as.matrix(cormatrix), type = 'spearman')$r

  colnames(test) <- conv %>% filter(OlinkID %in% colnames(test)) %>% arrange(match(OlinkID, colnames(test))) %>% pull(Assay)
  rownames(test) <- conv %>% filter(OlinkID %in% rownames(test)) %>% arrange(match(OlinkID, rownames(test))) %>% pull(Assay)

  g1 <- graph.adjacency(test, weighted = T, mode = 'lower')
  
  g1 <- simplify(g1)
  
  V(g1)$size <- 8
  V(g1)$frame.color <- 'black'
  V(g1)$color <- 'orange'
  V(g1)$label <- colnames(test)
  E(g1)$arrow.mode <- 0
  
  E(g1)$weight <- abs(E(g1)$weight)
  
  cut.off <- quantile(E(g1)$weight, probs = c(0.975))
  g1 <- delete_edges(g1, E(g1)[weight < cut.off])
  g1 <- simplify(g1)
  
  isolated = which(degree(g1)==0)
  g1 = delete.vertices(g1, isolated)
  
  hs <- hub_score(g1, weights=NA)$vector
  hubs <- names(tail(sort(degree(g1)), 3))
  
  hs[which(!names(hs) %in% hubs)] <- 1
  hs[which(names(hs) %in% hubs)] <- 15
  return (g1)
  
  # plot (g1, vertex.size=hs, vertex.label= ifelse(degree(g1) >= degree(g1)[hubs[1]], names(degree(g1)), NA))
  
  
}

ICU.network <- do.plot(df.ICU)
nonICU.network <- do.plot(df.nonICU)
postcovid.network <- do.plot(df.postcovid)
healthy.network <- do.plot(df.healthy)

createNetworkFromIgraph(ICU.network, 'ICU')
createNetworkFromIgraph(nonICU.network, 'nonICU')
createNetworkFromIgraph(postcovid.network, 'postcovid')
createNetworkFromIgraph(healthy.network, 'healthy')


pdf('output/ICU_network.pdf', width = 10, height = 10)
do.plot(df.ICU)
dev.off()

pdf('output/nonICU_network.pdf', width = 10, height = 10)
do.plot(df.nonICU)
dev.off()

pdf('output/postcovid_network.pdf', width = 10, height = 10)
do.plot(df.postcovid)
dev.off()

pdf('output/healthy_network.pdf', width = 10, height = 10)
do.plot(df.healthy)
dev.off()

###### end function 

cor.df.ICU <- cor(df.ICU, method = 'pearson')
cor.df.nonICU <- cor(df.nonICU, method = 'pearson')
cor.df.postcovid <- cor(df.postcovid, method = 'pearson')
cor.df.healthy <- cor(df.healthy, method = 'pearson')

ord.ICU <- pheatmap(cor.df.ICU)$tree_row$order
ord.nonICU <- pheatmap(cor.df.nonICU)$tree_row$order
ord.postcovid <- pheatmap(cor.df.postcovid)$tree_row$order
ord.healthy <- pheatmap(cor.df.healthy)$tree_row$order

dev.off()

do.heatmap <- function(cor.df.ICU, cor.df.nonICU, cor.df.postcovid, cor.df.healthy, ord, legend) {
  
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
      labs(x = 'Proteins', y = 'Proteins', title = 'COVID-19 ICU', fill = 'Correlation') +
      scale_fill_gradientn(colours = colours, limits = lims, oob = scales::squish),
    
    ggplot() +
      geom_tile(data = cor.df.nonICU.melt, aes(Var2, Var1, fill = value)) +
      labs(x = 'Proteins', y = 'Proteins', title = 'COVID-19 non-ICU', fill = 'Correlation') +
      scale_fill_gradientn(colours = colours, limits = lims, oob = scales::squish),
    
    ggplot() +
      geom_tile(data = cor.df.postcovid.melt, aes(Var2, Var1, fill = value)) +
      labs(x = 'Proteins', y = 'Proteins', title = 'Convalescent', fill = 'Correlation') +
      scale_fill_gradientn(colours = colours, limits = lims, oob = scales::squish),
    
    ggplot() +
      geom_tile(data = cor.df.healthy.melt, aes(Var2, Var1, fill = value)) +
      labs(x = 'Proteins', y = 'Proteins', title = 'Healthy', fill = 'Correlation') +
      scale_fill_gradientn(colours = colours, limits = lims, oob = scales::squish),
    nrow = 1, ncol = 4, common.legend = T, legend = legend
  )
  return(p)
}

pdf('output/corr_heatmap_ICU.pdf', width = 15, height = 3, onefile = F)
p1 <- do.heatmap(cor.df.ICU = cor.df.ICU, 
           cor.df.nonICU = cor.df.nonICU, 
           cor.df.postcovid = cor.df.postcovid, 
           cor.df.healthy = cor.df.healthy, 
           ord = ord.ICU, legend = 'none')
dev.off()

pdf('output/corr_heatmap_nonICU.pdf', width = 15, height = 3, onefile = F)
p2 <- do.heatmap(cor.df.ICU = cor.df.ICU, 
           cor.df.nonICU = cor.df.nonICU, 
           cor.df.postcovid = cor.df.postcovid, 
           cor.df.healthy = cor.df.healthy, 
           ord = ord.nonICU, legend = 'none')
dev.off()

pdf('output/corr_heatmap_postcovid.pdf', width = 15, height = 3, onefile = F)
p3 <- do.heatmap(cor.df.ICU = cor.df.ICU, 
           cor.df.nonICU = cor.df.nonICU, 
           cor.df.postcovid = cor.df.postcovid, 
           cor.df.healthy = cor.df.healthy, 
           ord = ord.postcovid, legend = 'none')
dev.off()

pdf('output/corr_heatmap_healthy.pdf', width = 15, height = 3, onefile = F)
p4 <- do.heatmap(cor.df.ICU = cor.df.ICU, 
           cor.df.nonICU = cor.df.nonICU, 
           cor.df.postcovid = cor.df.postcovid, 
           cor.df.healthy = cor.df.healthy, 
           ord = ord.healthy, legend = 'none')
dev.off()

# Empty plot we use just for creation of the legend grob
legend.grob <- get_legend(ggplot() + 
  geom_point(data = melt(cor.df.ICU),
             aes(Var1, Var2, fill = value)) +
  scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100), 
                       limits = c(-0.5, 0.5), breaks = c(-0.5, 0.0, 0.5)) +
    theme(legend.position="bottom") +
labs(fill = 'Correlation'))




pdf('output/heatmaps_combined.pdf', width = 8, height = 8)
ggarrange(p1, p2, p3, p4, ncol = 1, nrow = 4, legend.grob = legend.grob, legend = 'bottom')
dev.off()

# Write the correlation values to .csv files for cytoscape import
# cor.df.ICU, cor.df.nonICU, cor.df.postcovid and cor.df.healthy contain the correlations for each of the proteins.
# Remove the values that are below or above the 10th and 90th percentiles

# ICU
colnames(cor.df.ICU) <- conv %>% filter(OlinkID %in% colnames(cor.df.ICU)) %>% arrange(match(colnames(cor.df.ICU), OlinkID)) %>% pull(Assay)
rownames(cor.df.ICU) <- conv %>% filter(OlinkID %in% rownames(cor.df.ICU)) %>% arrange(match(rownames(cor.df.ICU), OlinkID)) %>% pull(Assay)

cor.df.ICU[upper.tri(cor.df.ICU, diag = T)] <- NA
cor.df.ICU.melt <- na.omit(melt(cor.df.ICU))
quants <- quantile(cor.df.ICU.melt$value, prob = c(0.05, 0.95))
cor.df.ICU.melt %<>% filter( (value < quants[1]) | (value > quants[2]) )
write.csv('output/gene_coexpr_network_ICU.csv', x = cor.df.ICU.melt)

# nonICU
colnames(cor.df.nonICU) <- conv %>% filter(OlinkID %in% colnames(cor.df.nonICU)) %>% arrange(match(colnames(cor.df.nonICU), OlinkID)) %>% pull(Assay)
rownames(cor.df.nonICU) <- conv %>% filter(OlinkID %in% rownames(cor.df.nonICU)) %>% arrange(match(rownames(cor.df.nonICU), OlinkID)) %>% pull(Assay)

cor.df.nonICU[upper.tri(cor.df.nonICU, diag = T)] <- NA
cor.df.nonICU.melt <- na.omit(melt(cor.df.nonICU))
quants <- quantile(cor.df.nonICU.melt$value, prob = c(0.10, 0.90))
cor.df.nonICU.melt %<>% filter( (value < quants[1]) | (value > quants[2]) )
write.csv('output/gene_coexpr_network_nonICU.csv', x = cor.df.nonICU.melt)

# Postcovid
colnames(cor.df.postcovid) <- conv %>% filter(OlinkID %in% colnames(cor.df.postcovid)) %>% arrange(match(colnames(cor.df.postcovid), OlinkID)) %>% pull(Assay)
rownames(cor.df.postcovid) <- conv %>% filter(OlinkID %in% rownames(cor.df.postcovid)) %>% arrange(match(rownames(cor.df.postcovid), OlinkID)) %>% pull(Assay)

cor.df.postcovid[upper.tri(cor.df.postcovid, diag = T)] <- NA
cor.df.postcovid.melt <- na.omit(melt(cor.df.postcovid))
quants <- quantile(cor.df.postcovid.melt$value, prob = c(0.10, 0.90))
cor.df.postcovid.melt %<>% filter( (value < quants[1]) | (value > quants[2]) )
write.csv('output/gene_coexpr_network_postcovid.csv', x = cor.df.postcovid.melt)

# Healthy
colnames(cor.df.healthy) <- conv %>% filter(OlinkID %in% colnames(cor.df.healthy)) %>% arrange(match(colnames(cor.df.healthy), OlinkID)) %>% pull(Assay)
rownames(cor.df.healthy) <- conv %>% filter(OlinkID %in% rownames(cor.df.healthy)) %>% arrange(match(rownames(cor.df.healthy), OlinkID)) %>% pull(Assay)

cor.df.healthy[upper.tri(cor.df.healthy, diag = T)] <- NA
cor.df.healthy.melt <- na.omit(melt(cor.df.healthy))
quants <- quantile(cor.df.healthy.melt$value, prob = c(0.10, 0.90))
cor.df.healthy.melt %<>% filter( (value < quants[1]) | (value > quants[2]) )
write.csv('output/gene_coexpr_network_healthy.csv', x = cor.df.healthy.melt)





