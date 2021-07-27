rm(list = ls())
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(openxlsx)

load('data/data.RData')
rm(annot.breda, breda, radboud, annot.radboud)

# How are the AB levels correlated to each other?
cor.test(annot.mhh$antibody.ncp, annot.mhh$antibody.S1, exact = T, method = 'pearson')

correlations.time <- matrix(nrow = ncol(mhh), ncol = 2, data = NA, dimnames = list( colnames(mhh), c('cor', 'pvalue')))

# Time of conv
for (protein in colnames(mhh)){
  
  df <- cbind(mhh %>% select(protein), 
              annot.mhh %>% select(time.convalescence))

  df <- df[complete.cases(df), ]
  x <- cor.test(df %>% pull(protein), 
                df %>% pull(time.convalescence), method = 'pearson')
  
  correlations.time[protein, 'cor'] <- x$estimate
  correlations.time[protein, 'pvalue'] <- x$p.value
}

# S1 antibody
correlations.S1 <- matrix(nrow = ncol(mhh), ncol = 2, data = NA, dimnames = list(colnames(mhh), c('cor', 'pvalue')))

for (protein in colnames(mhh)){
    
  df <- cbind(mhh %>% select(protein), 
              annot.mhh %>% select(antibody.S1))
  df <- df[complete.cases(df), ]
  x <- cor.test(df %>% pull(protein), 
                df %>% pull(antibody.S1))
  
  correlations.S1[protein, 'cor'] <- x$estimate
  correlations.S1[protein, 'pvalue'] <- x$p.value
}


# NCP antibody
correlations.NCP <- matrix(nrow = ncol(mhh), ncol = 2, data = NA, dimnames = list(colnames(mhh), c('cor', 'pvalue')))

for (protein in colnames(mhh)){
  
  df <- cbind(mhh %>% select(protein), 
              annot.mhh %>% select(antibody.ncp))
  
  df <- df[complete.cases(df), ]
  x <- cor.test(df %>% pull(protein), 
                df %>% pull(antibody.ncp))
  
  correlations.NCP[protein, 'cor'] <- x$estimate
  correlations.NCP[protein, 'pvalue'] <- x$p.value
}

# Posthoc correction
correlations.time <- data.frame(correlations.time)
correlations.S1 <- data.frame(correlations.S1)
correlations.NCP <- data.frame(correlations.NCP)

correlations.time$adj.p.val <- p.adjust(correlations.time$pvalue)
correlations.S1$adj.p.val <- p.adjust(correlations.S1$pvalue)
correlations.NCP$adj.p.val <- p.adjust(correlations.NCP$pvalue)

correlations.time['OID00443', ]


# Get the sig proteins
sig.time <- correlations.time %>% filter(adj.p.val < 0.05) %>% row.names()
sig.S1 <- correlations.S1 %>% filter(adj.p.val < 0.05) %>% row.names()
sig.NCP <- correlations.NCP %>% filter(adj.p.val < 0.05) %>% row.names()

# What are the values of theses proteins in ICU and nonICU samples?
DE.ICU <- openxlsx::read.xlsx('output/DE.xlsx', colNames=T, sheet = 'healthy_vs_ICU')
DE.nonICU <- openxlsx::read.xlsx('output/DE.xlsx', colNames=T, sheet = 'healthy_vs_nonICU')
DE.postcovid <- openxlsx::read.xlsx('output/DE.xlsx', colNames=T, sheet = 'healthy_vs_postcovid')

DE.ICU %>% filter(OlinkID %in% sig.time)
DE.ICU %>% filter(OlinkID %in% sig.S1)
DE.ICU %>% filter(OlinkID %in% sig.NCP)

DE.nonICU %>% filter(OlinkID %in% sig.time)
DE.nonICU %>% filter(OlinkID %in% sig.S1)
DE.nonICU %>% filter(OlinkID %in% sig.NCP)

DE.postcovid %>% filter(OlinkID %in% sig.time)
DE.postcovid %>% filter(OlinkID %in% sig.S1)
DE.postcovid %>% filter(OlinkID %in% sig.NCP)

# For each of the sig proteins, check correlation to age and gender

age_pvalues <- c()
gender_pvalues <- c()
vec <- c(sig.time, sig.S1, sig.NCP)
vec <- c(sig.S1, sig.NCP)
vec <- unique(vec)
for (prot in vec) {
  
  df <- cbind(mhh %>% select(prot), 
              annot.mhh %>% select(age, gender))
  
  # Age
  res <- cor.test(df %>% pull(prot), 
                  df %>% pull(age))
  if (res$p.value < 0.05) { print(paste0(prot, '/ age'))}
  age_pvalues <- c(res$p.value, age_pvalues)
  
  # gender
  df$gender <- as.numeric(as.factor(df$gender))
  res <- cor.test(df %>% pull(prot), 
                  df %>% pull(gender))
  if (res$p.value < 0.05) { print(paste0(prot, '/ gender'))}
  gender_pvalues <- c(res$p.value, gender_pvalues)
}
pvals <- data.frame(age = age_pvalues, gender = gender_pvalues, row.names = vec)
rownames(pvals) <- conv %>% filter(OlinkID %in% rownames(pvals)) %>% arrange(match(OlinkID, rownames(pvals))) %>% pull(Assay)
pvals$age_adj <- p.adjust(pvals$age, method = 'BH')
pvals$gender_adj <- p.adjust(pvals$gender, method = 'BH')

sum(pvals$age_adj < 0.05)
sum(pvals$gender_adj < 0.05)

pvals < 0.05





# Heatmap within the postcovid cohort only
# Time of conv.
load('data/data.RData')

# Select only patients that we know the time of conv for.
annot.mhh <- annot.mhh[which(!is.na(annot.mhh$time.convalescence)), ]
annot.mhh <- annot.mhh %>% arrange(desc(time.convalescence))
mhh <- mhh %>% na.omit()

annot.mhh <- annot.mhh[which(rownames(annot.mhh) %in% rownames(mhh)), ]
annot.mhh$sex <- annot.mhh$gender

mhh <- mhh[which(rownames(mhh) %in% rownames(annot.mhh)), ]
mhh <-   mhh %>%  arrange(match(rownames(mhh), rownames(annot.mhh)))

sub <- mhh %>% select(all_of(sig.time))

colnames(sub) <- conv %>% 
  filter(OlinkID %in% colnames(sub)) %>% 
  arrange(match(OlinkID, colnames(sub))) %>% 
  pull(Assay)

color_age <- colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100)
color_time <- colorRampPalette(brewer.pal(n = 7, name = "Greys"))(100)



pdf('output/heatmap_proteins_sig_time_conv.pdf', width = 20, height = 5)
pheatmap(mat = t(sub), 
         show_colnames = F, 
         scale = 'row', 
         breaks = seq(-1.5, 1.5, length.out = 100), 
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100), 
         cluster_cols = F, 
         annotation_col = annot.mhh %>% select(age, sex, time.convalescence) %>% 
           rename('Time of convalescence' = time.convalescence, 
                  'Sex' = sex, 
                  'Age' = age), 
         annotation_colors = list(
           'Sex' = c(male = 'steelblue', female = 'red'), 
           'Age' = color_age, 
           'Time of convalescence' = color_time
         ))
dev.off()

# Antibody levels
load('data/data.RData')

annot.mhh <- annot.mhh[which(!is.na(annot.mhh$antibody.S1), !is.na(annot.mhh$antibody.ncp)), ]

annot.mhh <- annot.mhh %>% arrange(antibody.ncp)
mhh <- mhh %>% na.omit()
annot.mhh$sex <- annot.mhh$gender

annot.mhh <- annot.mhh[which(rownames(annot.mhh) %in% rownames(mhh)), ]

mhh <- mhh[which(rownames(mhh) %in% rownames(annot.mhh)), ]
mhh <-   mhh %>%  arrange(match(rownames(mhh), rownames(annot.mhh)))

sub <- mhh %>% select(all_of(unique(c(sig.S1, sig.NCP))))

colnames(sub) <- conv %>% 
  filter(OlinkID %in% colnames(sub)) %>% 
  arrange(match(OlinkID, colnames(sub))) %>% 
  pull(Assay)

color_age <- colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100)
color_NCP <- colorRampPalette(brewer.pal(n = 7, name = "BuPu"))(100)
color_S1 <- colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100)

pdf('output/heatmap_proteins_sig_antibodies.pdf', width = 20, height = 5)
pheatmap(mat = t(sub), 
         show_colnames = F,
         scale = 'row', 
         breaks = seq(-1.5, 1.5, length.out = 100), 
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100), 
         cluster_cols = F, 
         annotation_col = annot.mhh %>% 
           select(age, sex, antibody.S1, antibody.ncp) %>% 
           rename('Antibody (S1)' = antibody.S1, 
                  'Antibody (NCP)' = antibody.ncp, 
                  'Sex' = sex, 
                  'Age' = age), 
         annotation_colors = list(
           'Sex' = c(male = 'steelblue', female = 'red'), 
           'Age' = color_age, 
           'Antibody (S1)' = color_S1, 
           'Antibody (NCP)' = color_NCP 
         )
)

dev.off()
