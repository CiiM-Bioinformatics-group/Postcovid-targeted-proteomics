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
ICU <- openxlsx::read.xlsx('output/DE.xlsx', sheet = 'ICU_vs_healthy')
nonICU <- openxlsx::read.xlsx('output/DE.xlsx', sheet = 'nonICU_vs_healthy')

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
#       FALSE TRUE
# FALSE    70   29
# TRUE     12   34
# 104 same direction, while 41 different
# 71% / 29%

# With sig. threshold. How many same dir?
ICU %<>% filter(significance == TRUE)
critical %<>% filter(OLINK.ID %in% ICU$OlinkID) %>% arrange(match(OLINK.ID, ICU$OlinkID))

table(ICU$logFC < 0, critical$logFC < 0)
#         FALSE TRUE
# FALSE    62   20
# TRUE      6   26
# 88 same dir, 26 not
# 77% same dir, 23% not

table(critical$adj.P.Val < 0.05) # 58 FALSE, 56 TRUE
# 56 out of 114 = 49%

#################
# non-ICU vs severe
# How many of the significant proteins in ICU can we replicate significantly and with the same direction in Patel data?
nonICU %<>% filter(OlinkID %in% severe$OLINK.ID) %>% arrange(match(OlinkID, severe$OLINK.ID))
severe %<>% filter(OLINK.ID %in% nonICU$OlinkID) %>% arrange(match(OLINK.ID, severe$OLINK.ID))
all(nonICU$OlinkID == severe$OLINK.ID)

# How many same direction without significance threshold?
table(nonICU$logFC < 0, severe$logFC < 0)
#         FALSE TRUE
# FALSE    43   54
# TRUE     9   39
# 82 same dir, 63 not
# 56% same dir, 45% not

nonICU %<>% filter(significance == TRUE)
severe %<>% filter(OLINK.ID %in% nonICU$OlinkID) %>% arrange(match(OLINK.ID, nonICU$OlinkID))

table(nonICU$logFC < 0, severe$logFC < 0)
#         FALSE TRUE
# FALSE    37   39
# TRUE      8   22
# 58 same dir, 45 not
# 54% same dir, 46% not

table(severe$adj.P.Val < 0.05)
# 12 out of 113 = 11%
