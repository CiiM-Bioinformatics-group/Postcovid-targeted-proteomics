rm(list = ls())
dev.off()

library(pROC)
library(caret)
library(glmnet)
library(dplyr)
library(magrittr)

# prediction for the separation between ICU and non-ICU
# Training: Breda cohort
# Test: Radboud cohort
setwd('/Users/martijnzoodsma/Documents/PhD/corona/Postcovid-targeted-proteomics/')
load('data/data.RData')

common <- intersect(colnames(breda), colnames(radboud))
breda %<>% select(all_of(common)) %>% na.omit()
radboud %<>% select(all_of(common)) %>% na.omit()

annot.radboud <- annot.radboud[which(rownames(annot.radboud) %in% rownames(radboud)), ]
annot.breda <- annot.breda[which(rownames(annot.breda) %in% rownames(breda)), ]

annot.breda %<>% na.omit()
annot.radboud %<>% filter(time == 'W1T1') %>% na.omit()

annot.breda[which(annot.breda$condition == 'non-ICU'), 'condition'] <- 'nonICU'
annot.radboud[which(annot.radboud$condition == 'non-ICU'), 'condition'] <- 'nonICU'

breda <- breda[which(rownames(breda) %in% rownames(annot.breda)), ]
radboud <- radboud[which(rownames(radboud) %in% rownames(annot.radboud)), ]

##########
# Construct training and validation data
# We use Breda as training cohort then validate in Radboud cohort
# Protein NPX values are scaled manually afterwards age and gender are added.
train_x <- breda %>% as.matrix()
train_x <- scale(train_x, center = F, scale = T)
annot.breda %>%
  select(age, gender) %>%
  mutate(age = as.numeric(age)) %>%
  mutate(gender = as.numeric(as.factor(gender))) %>%
  as.matrix() -> breda_age_gender
train_x <- cbind(train_x, breda_age_gender)
train_y <- factor(annot.breda$condition, levels = c('ICU', 'nonICU'))

validation_x <- radboud %>% as.matrix()
validation_x <- scale(validation_x, center = F, scale = T)
annot.radboud %>%
  select(age, gender) %>%
  mutate(age = as.numeric(age)) %>%
  mutate(gender = as.numeric(as.factor(gender))) %>%
  as.matrix() -> radboud_age_gender
validation_x <- cbind(validation_x, radboud_age_gender)
validation_y <- factor(annot.radboud$condition, levels = c('ICU', 'nonICU'))

## Make the model
# Train 100 times. Each time, save the coefficients so we can see which proteins were set to zero and which are included in the model.

n <- 1
res <- matrix(NA, nrow = n, ncol = 1)
coefficients <- data.frame(row.names = c('Intercept', colnames(train_x)))

for (i in 1:n) {
  netFit <- train(x = train_x, 
                  y = train_y,
                  method = "glmnet", 
                  metric = "ROC", 
                  tuneGrid=expand.grid(.alpha = seq(0.1,0.9, by=0.1), 
                                       .lambda = seq(0,1,by=0.01)),
                  trControl = trainControl(method="repeatedcv", 
                                           number=10, 
                                           repeats=5, 
                                           classProbs = T, 
                                           summaryFunction=twoClassSummary))

  # Accuracy and training perf
  print(paste0('Training performance: ', getTrainPerf(netFit)))
  
  predict_validation <- predict(object = netFit, newdata = validation_x, type = 'raw')
  conf <- confusionMatrix(data = predict_validation, reference = validation_y)
  acc <- conf$overall['Accuracy']
  res[i, 1] <- acc

  # Get the final model and extract coefficients
  final.model <- netFit$finalModel
  coefficients_condition <- coef(object = final.model, s = netFit$bestTune$lambda)
  coefficients_condition <- matrix(coefficients_condition)
  coefficients_condition <- data.frame(coefficients_condition)
  colnames(coefficients_condition) <- paste0('run', i)
  coefficients <- cbind(coefficients, coefficients_condition)
}

saveRDS(object = netFit, file = 'output/prediction_model_condition.RDS')
write.csv(res, file = 'output/pred_condition_res.csv')
write.csv(coefficients, file = 'output/pred_condition_coefficients.csv')

# In both training and validation phenotypes:
  # ICU = 1, non-ICU = 2. We set 2 as the reference level since this is also what we do in the Limma DE

### Training performance
# AUC on the training data
training.predict <- predict(object = netFit, newdata = train_x)
conf <- confusionMatrix(data = training.predict, reference = train_y)
conf
training_conf <- conf$table %>% as.data.frame()

roc.data <- roc(as.numeric(as.factor(training.predict)), as.numeric(as.factor(train_y)))
auc(roc.data)
ci.auc(roc.data)

### Validation data
# ROC curve and get AUC information along with prediction accuracy
pred <- predict(object = netFit, newdata = validation_x, type = 'prob')[[2]]
roc.data <- roc(as.numeric(validation_y), pred)

pred.raw <- predict(netFit, validation_x, type = 'raw')
conf <- confusionMatrix(data = pred.raw, reference = validation_y)
conf
validation_conf <- conf$table %>% as.data.frame()
auc(roc.data) #0.87
auc(roc.data, )
pdf('output/roc_prediction_ICU_nonICU.pdf', width = 3, height = 3)
roc <- ggroc(roc.data) +
  theme_classic() +
  geom_abline(intercept = 1, slope = 1, lty = 2) +
  labs(x = '1 - Specificity', y = 'Sensitity') +
  annotate("text", x=0.5, y=0.75, label="AUC: 0.87\nCI: 0.79 - 0.94", color = "black", size = 3)

roc
dev.off()
