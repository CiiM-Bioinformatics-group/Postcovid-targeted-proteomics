rm(list = ls())

library(openxlsx)
library(dplyr)


# Inflammation panel
pl1 <- read.xlsx('/Users/martijnzoodsma/Documents/PhD/corona/data/Nijmegen cohort/integration/inflammation panel/integration_INF_pl1.xlsx', rowNames = T, colNames = T, sheet = 'Radboud final', sep.names = ' ')
pl2 <- read.xlsx('/Users/martijnzoodsma/Documents/PhD/corona/data/Nijmegen cohort/integration/inflammation panel/integration_INF_pl2.xlsx', rowNames = T, colNames = T, sheet = 'Radboud final', sep.names = ' ')
pl3 <- read.xlsx('/Users/martijnzoodsma/Documents/PhD/corona/data/Nijmegen cohort/integration/inflammation panel/integration_INF_pl3-7.xlsx', rowNames = T, colNames = T, sheet = 'Radboud final', sep.names = ' ')

# Common proteins between all three panels only
total.common <- Reduce(intersect, list(colnames(pl1), colnames(pl2), colnames(pl3)))

# Subset and reorder
pl1 <- pl1[, which(colnames(pl1) %in% total.common)]
pl2 <- pl2[, which(colnames(pl2) %in% total.common)]
pl3 <- pl3[, which(colnames(pl3) %in% total.common)]

pl1 <- pl1[total.common]
pl2 <- pl2[total.common]
pl3 <- pl3[total.common]

panel_inflammation <- rbind(pl1, pl2, pl3)


# Cardiometabolic panel
pl1 <- read.xlsx('/Users/martijnzoodsma/Documents/PhD/corona/data/Nijmegen cohort/integration/cardiometabolic panel/integration_radboud_CM_plate1.xlsx', rowNames=T, colNames=T, sheet = 'Radboud final', sep.names = ' ')
pl2 <- read.xlsx('/Users/martijnzoodsma/Documents/PhD/corona/data/Nijmegen cohort/integration/cardiometabolic panel/integration_radboud_CM_plate2.xlsx', rowNames=T, colNames=T, sheet = 'Radboud final', sep.names = ' ')
pl3 <- read.xlsx('/Users/martijnzoodsma/Documents/PhD/corona/data/Nijmegen cohort/integration/cardiometabolic panel/integration_radboud_CM_plate3-7.xlsx', rowNames=T, colNames=T, sheet = 'Radboud final', sep.names = ' ')

# Common proteins between all three panels only
total.common <- Reduce(intersect, list(colnames(pl1), colnames(pl2), colnames(pl3)))

# Subset and reorder
pl1 <- pl1[, which(colnames(pl1) %in% total.common)]
pl2 <- pl2[, which(colnames(pl2) %in% total.common)]
pl3 <- pl3[, which(colnames(pl3) %in% total.common)]

pl1 <- pl1[total.common]
pl2 <- pl2[total.common]
pl3 <- pl3[total.common]

panel_cardiometabolic <- rbind(pl1, pl2, pl3)


# Cardiovascular panel
pl1 <- read.xlsx('/Users/martijnzoodsma/Documents/PhD/corona/data/Nijmegen cohort/integration/cardiovascular panel/integration_radboud_CV_plate1.xlsx', rowNames = T, colNames = T, sheet = 'Radboud final', sep.names = ' ')
pl2 <- read.xlsx('/Users/martijnzoodsma/Documents/PhD/corona/data/Nijmegen cohort/integration/cardiovascular panel/integration_radboud_CV_plate2.xlsx', rowNames = T, colNames = T, sheet = 'Radboud final', sep.names = ' ')
pl3 <- read.xlsx('/Users/martijnzoodsma/Documents/PhD/corona/data/Nijmegen cohort/integration/cardiovascular panel/integration_radboud_CV_plates3-7.xlsx', rowNames = T, colNames = T, sheet = 'Radboud final', sep.names = ' ')


# Common proteins between all three panels only
total.common <- Reduce(intersect, list(colnames(pl1), colnames(pl2), colnames(pl3)))

# Subset and reorder
pl1 <- pl1[, which(colnames(pl1) %in% total.common)]
pl2 <- pl2[, which(colnames(pl2) %in% total.common)]
pl3 <- pl3[, which(colnames(pl3) %in% total.common)]

pl1 <- pl1[total.common]
pl2 <- pl2[total.common]
pl3 <- pl3[total.common]

panel_cardiovascular <- rbind(pl1, pl2, pl3)

# output
write.xlsx(x = panel_inflammation, file='/Users/martijnzoodsma/Documents/PhD/corona/Postcovid-targeted-proteomics/data/Radboud/merged_data/inflammation.xlsx', rowNames=T)
write.xlsx(x = panel_cardiovascular, file = '/Users/martijnzoodsma/Documents/PhD/corona/Postcovid-targeted-proteomics/data/Radboud/merged_data/cardiovascular.xlsx', rowNames=T)
write.xlsx(x = panel_cardiometabolic, file = '/Users/martijnzoodsma/Documents/PhD/corona/Postcovid-targeted-proteomics/data/Radboud/merged_data/cardiometabolic.xlsx', rowNames=T)
