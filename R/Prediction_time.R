rm(list = ls())
dev.off()

library(pROC)
library(caret)
library(glmnet)
library(dplyr)
library(magrittr)

# Prediction for the time of convalescence within MHH postcovid samples
load('../../../data/data.RData')

mhh <- na.omit(mhh)
annot.mhh %<>% filter(condition == 'postcovid') %>% 
  select(time.convalescence, age, gender) %>% 
  mutate(age = as.numeric(age)) %>% 
  mutate(gender = as.numeric(as.factor(gender)))
annot.mhh %<>% na.omit() 

mhh <- mhh[which(rownames(mhh) %in% rownames(annot.mhh)), ]
annot.mhh <- annot.mhh[which(rownames(annot.mhh) %in% rownames(mhh)), ]

all(rownames(mhh) == rownames(annot.mhh))

n <- 1
res <- matrix(NA, nrow = n, ncol = 2)
colnames(res) <- c('cor', 'pval')
coefficients <- data.frame(row.names = c('Intercept', colnames(mhh), 'age', 'gender'))

for(i in 1:n) {
  splitSample <- createDataPartition(annot.mhh$time.convalescence, p = 0.7, list = FALSE)
  
  # Train data
  train_x <- mhh[splitSample, ]
  train_x <- scale(train_x, center=F, scale=T)
  stopifnot(all(rownames(train_x) == rownames(annot.mhh[splitSample, ])))
  train_x <- cbind(train_x, 
                   as.matrix(annot.mhh[splitSample, c('age', 'gender')]))

  train_y <- annot.mhh[splitSample, ] %>% pull(time.convalescence)
  
  # Validation data
  validation_x <- mhh[-splitSample, ]
  validation_x <- scale(validation_x, center=F, scale = T)
  stopifnot(all(rownames(validation_x) == rownames(annot.mhh[-splitSample, ])))
  validation_x <- cbind(validation_x, 
                   as.matrix(annot.mhh[-splitSample, c('age', 'gender')]))
  
  validation_y <- annot.mhh[-splitSample, ] %>% pull(time.convalescence)
  
  netFit <- train(x = train_x, 
                  y = train_y,
                  method = "glmnet", metric = 'RMSE',
                  tuneGrid=expand.grid(.alpha = seq(0.1,0.9, by=0.1), 
                                       .lambda = seq(0,1,by=0.01)),
                  trControl = trainControl(method="repeatedcv", 
                                           number=10, 
                                           repeats=5))

  prediction_y <- predict(object = netFit, newdata = validation_x, type = 'raw')
  
  cor <- cor.test(validation_y, prediction_y)
  res[i, 'cor'] <- cor$estimate
  res[i, 'pval'] <- cor$p.value
  
  # Get the final model and extract coefficients
  final.model <- netFit$finalModel
  coefficients_time <- coef(object = final.model, s = netFit$bestTune$lambda)
  coefficients_time <- matrix(coefficients_time)
  coefficients_time <- data.frame(coefficients_time)
  colnames(coefficients_time) <- paste0('run', i)
  coefficients <- cbind(coefficients, coefficients_time)
}

saveRDS(object = netFit, file = 'output/prediction_model_convalescence_time.RDS')
write.csv(res, file = 'output/pred_convalescence_res.csv')
write.csv(coefficients, file = 'output/pred_convalescence_coefficients.csv')

# Correlation plot for the last run
pdf('output/correlation_time_conv.pdf', width = 3, height = 3)
corplot <- ggplot() +
  theme_classic() +
  geom_point(aes(x = validation_y, y = prediction_y)) +
  labs(x = 'Reported convalescence time', y = 'Predicted convalescence time') +
  xlim(c(15, 76)) +
  ylim(c(15, 50))
corplot
dev.off()
