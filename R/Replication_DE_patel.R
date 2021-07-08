rm(list = ls())
try(dev.off())

library(openxlsx)
library(dplyr)
library(magrittr)
library(ggplot2)

conv <- read.xlsx('data/convert_olink_names.xlsx')


# Patel et al data
patel <- read.csv('data/Patel et al//tabula-supplementary information Patel et al 2021 _ subset.csv', header=T, sep = ',')
severe <- patel %>% filter(Comparison == "Control_v_Severe") 
critical <- patel %>% filter(Comparison == "Control_v_Critical") 

# Own DE results
ICU <- openxlsx::read.xlsx('output/DE.xlsx', sheet = 'healthy_vs_ICU')
nonICU <- openxlsx::read.xlsx('output/DE.xlsx', sheet = 'healthy_vs_nonICU')

# Comparisons:
  # ICU vs Healthy compared to Critical vs Healthy
  # nonICU vs healthy compared to Severe vs healthy

# ICU vs critical
# How many of the significant proteins in ICU can we replicate significantly and with the same direction in Patel data?
ICU %<>% filter(OlinkID %in% critical$OLINK.ID) %>% arrange(match(OlinkID, critical$OLINK.ID))
critical %<>% filter(OLINK.ID %in% ICU$OlinkID) %>% arrange(match(OLINK.ID, ICU$OlinkID))
all(ICU$OlinkID == critical$OLINK.ID)

# How many same direction? Without significance threshold
table(ICU$logFC < 0, critical$logFC < 0)
#        FALSE TRUE
# FALSE    69   29
# TRUE     13   34
# 103 same direction, while 42 different
# 59 / 41%

ICU %<>% filter(significance == TRUE)
critical %<>% filter(OLINK.ID %in% ICU$OlinkID) %>% arrange(match(OLINK.ID, ICU$OlinkID))

table(ICU$logFC < 0, critical$logFC < 0)
#         FALSE TRUE
# FALSE    62   20
# TRUE      6   25
# 87 same dir, 26 not
# 61% same dir, 29% not

table(critical$adj.P.Val < 0.05) # 58 FALSE, 55 TRUE
# 55 out of 113 = 48%

#################
# non-ICU vs severe
# How many of the significant proteins in ICU can we replicate significantly and with the same direction in Patel data?
nonICU %<>% filter(OlinkID %in% severe$OLINK.ID) %>% arrange(match(OlinkID, severe$OLINK.ID))
severe %<>% filter(OLINK.ID %in% nonICU$OlinkID) %>% arrange(match(OLINK.ID, severe$OLINK.ID))
all(nonICU$OlinkID == severe$OLINK.ID)

# How many same direction without significance threshold?
table(nonICU$logFC < 0, severe$logFC < 0)
#         FALSE TRUE
# FALSE    42   56
# TRUE     10   37


nonICU %<>% filter(significance == TRUE)
severe %<>% filter(OLINK.ID %in% nonICU$OlinkID) %>% arrange(match(OLINK.ID, nonICU$OlinkID))

table(nonICU$logFC < 0, severe$logFC < 0)
#         FALSE TRUE
# FALSE    37   38
# TRUE      7   21
# 58 same dir, 45 not
# 56% same dir, 44% not

table(severe$adj.P.Val < 0.05)
# 12 out of 113 = 11%

