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
# Prediction for the time of convalescence within MHH postcovid samples
load('data.RData')

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

saveRDS(object = netFit, file = 'prediction_model_convalescence_time.RDS')
write.csv(res, file = 'pred_convalescence_res.csv')
write.csv(coefficients, file = 'pred_convalescence_coefficients.csv')

res <- data.frame(res)

# Correlation plot for the last run
corplot <- ggplot() +
  theme_classic() +
  geom_point(aes(x = validation_y, y = prediction_y)) +
  labs(x = 'Reported convalescence time', y = 'Predicted convalescence time')

cors <- ggplot(data = res) +
  geom_boxplot(aes(y = cor, x = as.factor('test')), width = 0.5) +
  geom_jitter(aes(y = cor, x = as.factor('test')), alpha = 0.3, width = 0.1) +
  theme_classic() +
  labs(y = 'Correlation reported / predicted') +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())

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

pdf('all.pdf', width = 10, height = 3)
ggarrange(corplot, cors, coefs, align = 'h', nrow = 1, widths = c(2, 1, 2))
dev.off()
