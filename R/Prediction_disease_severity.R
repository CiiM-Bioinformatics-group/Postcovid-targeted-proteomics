rm(list = ls())
try(dev.off())

library(pROC)
library(caret)
library(glmnet)
library(dplyr)
library(magrittr)

# prediction for the separation between ICU and non-ICU
# Training: Breda cohort
# Test: Radboud cohort

load('data.RData')

annot.radboud %<>% select(age, gender, condition, time, cohort)
annot.breda %<>% select(age, gender, condition, timepoint, cohort)

common <- intersect(colnames(breda), colnames(radboud))
breda %<>% select(all_of(common)) %>% na.omit()
radboud %<>% select(all_of(common)) %>% na.omit()

annot.radboud <- annot.radboud[which(rownames(annot.radboud) %in% rownames(radboud)), ]
annot.breda <- annot.breda[which(rownames(annot.breda) %in% rownames(breda)), ]

annot.breda %<>% filter(timepoint == 'T1') %>% na.omit()
annot.radboud %<>% filter(time == 'W1T1') %>% na.omit()

annot.breda[which(annot.breda$condition == 'non-ICU'), 'condition'] <- 'nonICU'
annot.radboud[which(annot.radboud$condition == 'non-ICU'), 'condition'] <- 'nonICU'

breda <- breda[which(rownames(breda) %in% rownames(annot.breda)), ]
radboud <- radboud[which(rownames(radboud) %in% rownames(annot.radboud)), ]

#########
# Construct training and validation data
# Use Breda as training cohort, validation in Radboud cohort
# Scale using caret preprocess and use the preprocess model to train the Radboud cohort too
annot.breda$age <- as.numeric(annot.breda$age)
annot.breda$gender <- as.factor(annot.breda$gender)

train_x <- cbind(breda, 
                 annot.breda %>% select(age,gender))

prep_model <- preProcess(x = train_x) # Gender is ignored in the model since it is a 2-level factor.
train_x <- predict(prep_model, newdata = train_x)
train_y <- factor(annot.breda$condition, levels = c('ICU', 'nonICU'))


annot.radboud$age <- as.numeric(annot.radboud$age)
annot.radboud$gender <- as.factor(annot.radboud$gender)
validation_x <- cbind(radboud, 
                      annot.radboud %>% select(age, gender))

stopifnot(all(colnames(train_x) == colnames(validation_x)))
validation_x <- predict(prep_model, newdata = validation_x)
validation_y <- factor(annot.radboud$condition, levels = c('ICU', 'nonICU'))

train_x$gender <- as.numeric(train_x$gender)
validation_x$gender <- as.numeric(validation_x$gender)


## Make the model
n <- 1
res <- matrix(NA, nrow = n, ncol = 1)
coefficients <- data.frame(row.names = c('Intercept', colnames(train_x)))

for (i in 1:n) {
  print(i)
  train(x = train_x, 
        y = train_y, 
        method = 'glmnet', 
        metric = 'ROC', 
        trControl = trainControl(classProbs = T, 
                                 summaryFunction = twoClassSummary, 
                                 method = 'repeatedcv', 
                                 number = 5, 
                                 repeats = 10), 
        tuneGrid=expand.grid(.alpha = seq(0.1,0.9, by=0.1),
                             .lambda = seq(0,1,by=0.01))
  ) -> netFit
  
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

saveRDS(object = netFit, file = 'prediction_model_condition.RDS')
write.csv(res, file = 'pred_condition_res.csv')
write.csv(coefficients, file = 'pred_condition_coefficients.csv')

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

### Validation data
# ROC curve and get AUC information along with prediction accuracy
pred <- predict(object = netFit, newdata = validation_x, type = 'prob')[[2]]
roc.data <- roc(as.numeric(validation_y), pred)

pred.raw <- predict(netFit, validation_x, type = 'raw')
conf <- confusionMatrix(data = pred.raw, reference = validation_y)
conf
validation_conf <- conf$table %>% as.data.frame()
print('AUC: ')
print(auc(roc.data)) #0.87

pdf('roc_prediction_ICU_nonICU.pdf', width = 3, height = 3)
roc <- ggroc(roc.data) +
  theme_classic() +
  geom_abline(intercept = 1, slope = 1, lty = 2) +
  labs(x = '1 - Specificity', y = 'Sensitity') +
  annotate("text", x=0.65, y=0.750, label="AUC: 0.87\n[0.79 - 0.95]", color = "black", size = 3)

roc
dev.off()


# Coefficients plot
df <- read.csv('/vol/projects/mzoodsma/postcovid/new_rebuttal/dis_sev/pred_condition_coefficients.csv', header=T, row.names=1)
# df <- coefficients
df <- df[which(!rownames(df) == 'Intercept'), ]

df2 <- data.frame(
  'mean' = apply(df, 1, mean),
  'sd' = apply(df, 1, sd),
  'abs_sd' = apply(abs(df), 1, sd),
  'abs_mean' = abs(apply(df, 1, mean))
)

df <- cbind(df, df2)

df <- df %>%
  arrange(desc(abs_mean)) %>%
  select(mean, sd, abs_mean) %>%
  head(10)

conv[which(conv$OlinkID %in% rownames(df)), ]
df
rownames(df) <- c('HGF', 'TGF-alpha', 'Age', 'CXCL9', 'TWEAK', 'CXCL11', 'CST5', 'uPA', 'TNFRSF9', 'MMP-10')


df$protein <- rownames(df)
df %<>% arrange(abs_mean)
df$protein <- factor(df$protein, levels = unique(df$protein))

pdf('/vol/projects/mzoodsma/postcovid/new_rebuttal/output/coefficients_pred_disease.pdf', width =3, height = 3)
coefs <- ggplot(data = df) +
  geom_point(aes(x = abs_mean, y = protein), size= 2) +
  geom_errorbar(aes(x = abs_mean, y = protein, xmax = abs_mean-sd, xmin=abs_mean + sd, width = 0.2)) +
  theme_classic() +
  labs(x = 'Absolute coefficient') +
  xlim(c(-.2, NA))
coefs
dev.off()

pdf('/vol/projects/mzoodsma/postcovid/new_rebuttal/output/disease_pred_combplot.pdf', width = 6, height = 3)
ggpubr::ggarrange(roc, coefs, align = 'h')
dev.off()
