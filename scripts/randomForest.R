library(randomForest)
library(tidyverse)
library(data.table)
library(pROC)
library(forestplot)


set.seed(20240625) 
indices <- sample(1:nrow(data), nrow(data)*0.7) 
trainingset <- data[indices, ]
testingset <- data[-indices, ]


n = ncol(trainingset) - 1

rate = 1
for (i in 1:n) {
    model = randomForest(label~., data = trainingset, mtry = i, importance = TRUE, ntree = 1000)
    rate[i] = mean(model$err.rate) 
}
rate
min(rate)

rf_model <- randomForest(label ~ ., data = trainingset, 
      ntree = 1000,
      mtry = 9,
        importance = TRUE,
        proximity = TRUE)

varImpPlot(rf_model)
rf_model$importance

plot(rf_model, main = "Random Forest Model Error", type = "o")


pre_ran <- predict(rf_model,newdata=trainingset)
obs_p_ran = data.frame(Pred=pre_ran,True=trainingset$label)
table(obs_p_ran) %>% t()


pre_ran <- predict(rf_model,newdata=testingset, type = "prob")
obs_p_ran = data.frame(Pred=pre_ran,True=testingset$label)
tail(obs_p_ran)


ran_roc <- roc(obs_p_ran$True, obs_p_ran$Pred.abTconst_gdTflux)

save(rf_model, file = "/work/haoq/TCGA/randomForest/random_forest_model_0.984.RData")