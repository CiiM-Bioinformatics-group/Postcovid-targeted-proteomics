rm(list = ls())
try(dev.off())

library(pROC)
library(caret)
library(glmnet)
library(dplyr)
library(ggpubr)
library(magrittr)
library(doMC)

registerDoMC(cores = 25)

# Prediction for the time of convalescence within MHH cohort: postcovid samples
setwd('/vol/projects/mzoodsma/postcovid/new_rebuttal/conv/')
load('data.RData')

# we only need the convalescent samples for here. Filter
# We also dont have the time of conv for 11 samples. Take these out
annot.mhh %<>% filter(condition == 'postcovid') %>%
  select(time.convalescence, age, gender) %>%
  mutate(age = as.numeric(age)) %>%
  mutate(gender = as.factor(gender))
annot.mhh <- na.omit(annot.mhh)
sum(is.na(annot.mhh))

mhh <- mhh[which(rownames(mhh) %in% rownames(annot.mhh)), ]
stopifnot(all(rownames(mhh) == rownames(annot.mhh)))

sum(is.na(mhh)) #525 missing values. Knn imputation using caret preProcess
data_total <- cbind(mhh, 
                    annot.mhh %>% select(age, gender))
# pred_model <- preProcess(x = data_total, method = 'knnImpute')

# data_total <- predict(pred_model, newdata = data_total)
phenotypes <- annot.mhh$time.convalescence

n <- 100
res <- matrix(NA, nrow = n, ncol = 2)
colnames(res) <- c('cor', 'pval')
coefficients <- data.frame(row.names = c('Intercept', colnames(mhh), 'age', 'gender'))

for(i in 1:n) {
  splitSample <- createDataPartition(phenotypes, p = 0.7, list = FALSE)
  
  # Train data
  train_x <- data_total[splitSample, ]
  train_y <- phenotypes[splitSample]
  
  # Validation data
  validation_x <- data_total[splitSample, ]
  validation_y <- phenotypes[splitSample]
  
  # Preprocessing
  prep_model <- preProcess(train_x, method = 'knnImpute')
  
  train_x <- predict(prep_model, train_x)
  validation_x <- predict(prep_model, validation_x)
  
  train_x$gender <- as.numeric(train_x$gender)
  validation_x$gender <- as.numeric(validation_x$gender)
  
  netFit <- train(x = train_x,
                  y = train_y,
                  method = "glmnet", metric = 'RMSE',
                  tuneGrid=expand.grid(.alpha = seq(0.1,0.9, by=0.1),
                                       .lambda = seq(0,1,by=0.01)),
                  trControl = trainControl(method="repeatedcv",
                                           number=5,
                                           repeats=30))
  
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
  print('done')
}

#saveRDS(object = netFit, file = 'prediction_model_convalescence_time.RDS')
#write.csv(res, file = 'pred_convalescence_res.csv')
#write.csv(coefficients, file = 'pred_convalescence_coefficients.csv')

res <- data.frame(res)

# Correlation plot for the last run
pdf('/vol/projects/mzoodsma/postcovid/new_rebuttal/output/correlation_conv.pdf', width = 3, height = 3)
corplot <- ggplot() +
  theme_classic() +
  geom_point(aes(x = validation_y, y = prediction_y)) +
  labs(x = 'Reported convalescence time', y = 'Predicted convalescence time') +
  annotate(x = 50, y = 30, geom = 'text', label= paste0("Corr: ", round(cor(validation_y, prediction_y), 3)))

corplot +
  theme(aspect.ratio = 1)
dev.off()

cors <- ggplot(data = res) +
  geom_boxplot(aes(y = cor, x = as.factor('test')), width = 0.5) +
  geom_jitter(aes(y = cor, x = as.factor('test')), alpha = 0.3, width = 0.1) +
  theme_classic() +
  labs(y = 'Correlation reported / predicted') +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  ylim(c(0.5, 1))

# Coefficients plot
df <- coefficients
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

rownames(df) <- conv %>%
  filter(OlinkID %in% rownames(df)) %>%
  arrange(match(OlinkID, rownames(df))) %>%
  pull(Assay)

df$protein <- rownames(df)
df %<>% arrange(abs_mean)
df$protein <- factor(df$protein, levels = unique(df$protein))

coefs <- ggplot(data = df) +
  geom_point(aes(x = abs_mean, y = protein), size= 2) +
  geom_errorbar(aes(x = abs_mean, y = protein, xmax = abs_mean-sd, xmin=abs_mean + sd, width = 0.2)) +
  theme_classic() +
  labs(x = 'Absolute coefficient')

pdf('/vol/projects/mzoodsma/postcovid/new_rebuttal/output/all_conv.pdf', width = 8, height = 3)
ggarrange(corplot, cors, coefs, align = 'h', nrow = 1, widths = c(2, 1, 2))
dev.off()