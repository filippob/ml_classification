## Rscript used for the development of the machine learning classifiers used to analyse the data for the article
## "Machine learning classification of archaea and bacteria identifies novel predictive genomic features"
## submitted to BMC Genomics 2024

# Load libraries
library("e1071")
library("caret")
library("glmnet")
library("ggplot2")
library("tidyverse")
library("data.table")
library("doParallel")
library("randomForest")
library("caretEnsemble")


library(formatR)

library(psych)
library(ggfortify)


library(party)

library(kernlab)
library(readxl)
library(flextable)

library(MLmetrics)

library(pROC)
library(glmnet)

# Import dataset

basefolder = "/home/filippo/Documents/tania/archaea_bacteria"
nproc = 4
nrepeats_rfe = 3
nrepeats_cv = 10
#set.seed(113)

fname = file.path(basefolder, "data/total_final_nozero.csv")
dataset <- fread(fname, dec = ",")

### THINNING #######################
## comment  the code below
## to use the entire dataset
dataset <- dataset |> 
  sample_frac(size = 0.1, replace = FALSE)
######################################

## first four columns are metadata
ids = dataset[,c(1:4)] ## saving metadata
dataset <- dataset %>% 
  mutate(Outcome = as.factor(Outcome)) %>%
  select(-c(Organism, Taxon)) %>%
  as.data.frame()

## sanity check
if(sum(is.na(dataset)) == 0) print("No missing data in the dataset: OK!")

# Setup parallel back end to use multiple processors
cl <- makeCluster(nproc)
registerDoParallel(cl)

###########################
# Machine learning analysis
###########################

## 1) Creation of train/test and validation datasets
training_index <- createDataPartition(dataset$Outcome, p=0.80, list=FALSE)
# Select 20% of the data for test
test <- dataset[-training_index,]
# Use the remaining 80% of data to train and test the models
training <- dataset[training_index,]

## sanity check on numbers
test %>% count(Outcome) |> print()
training %>% count(Outcome) |> print()

## 2) Feature selection (RFE): automatically selecting a subset of the most predictive
writeLines(" - running RFE")

# 10-fold cross validation (cv) repeated 100 times
# Use different seed for cv in rfe and train
rfectrl_rf <- rfeControl(functions=rfFuncs, method="repeatedcv", number=10, repeats=nrepeats_rfe, 
                         verbose = FALSE, returnResamp = "final", allowParallel = TRUE) 

# After indicating predictors(3-90) and outcome (1), excluding id (2),
# check from 1 to n possible combinations of predictors to define the best subset
nvar = ncol(training)
results_rfe <- rfe(training[,3:nvar], training$Outcome, sizes=c(1:(nvar-2)), rfeControl=rfectrl_rf)
print(results_rfe$fit)

# Plot the results
#plot(results_rfe, type=c("g", "o"))
png(file.path(basefolder, "rfe.png"))
plot(results_rfe, type=c("g", "o"))
dev.off()

# List the chosen features
print("Features selected from RFE")
print(predictors(results_rfe))

## 3) Prepare datasets (train/test= dataset3 and validation), excluding non-selected features
selected <- predictors(results_rfe) # list of selected features
training_reduced <- training[,c("Outcome", "ID", selected)] # keep only selected features from dataset3 
training_reduced <- as.data.frame(training_reduced) # save as dataframe
test_reduced <- test[,c("Outcome", "ID", selected)] # keep only selected features from validation 
test_reduced <- as.data.frame(test_reduced) #save as dataframe

## 4) Model training
writeLines(" - Model training")
res.list <- list()

for(i in 1:nrepeats_cv) {
  
  print(paste("repetition", i))	
  folds <- createMultiFolds(training_reduced$Outcome, k = 10, times = 1) #stratification
  
  trainctrl <- trainControl(method="repeatedcv", number=10, repeats=1, index = folds,
                          classProbs= TRUE, savePredictions = TRUE, allowParallel = TRUE, sampling = "up") 
  metric <- "Accuracy" 
  metodi <- c('glmnet', 'rf', 'svmRadial', 'nnet')
  
  fit.models <- caretList(Outcome ~ ., data=training_reduced[,-2], trControl=trainctrl, preProcess = c("center", "scale"),
                          methodList=metodi, metric=metric) 
  
  res.list[[i]] <- fit.models 
  
}


## 5) gathering results
# Create a data.frame with accuracy and Kappa of the 10 fold of all replicates for each method
ALL_RESULTS <- NULL

for (i in seq(length(res.list))){
  print(i)
  summary(resamples(res.list[[i]]))[[1]] %>% mutate(REPLICA=i)-> X
  bind_rows(ALL_RESULTS,X) ->ALL_RESULTS
}

# Descriptive statistics for accuracy and kappa of each method
ALL_RESULTS %>%
  select(-REPLICA) %>%
  summarise(
    across(everything(),list(media=~mean(.x,na.rm=T), 
                             median=~median(.x,na.rm=T),
                             min=~min(.x,na.rm=T),
                             max=~max(.x,na.rm=T),
                             Q1=~quantile(.x,.25,na.rm=T),
                             Q3=~quantile(.x,.75,na.rm=T)),
           
           .names = "{.col}_{.fn}")
  ) -> STAT_ALL

# Final tables of descriptive statistics for accuracy e kappa
STAT_ALL %>%
  select(contains('Accuracy')) %>%
  pivot_longer(everything(),values_to = 'stat',names_to = 'param') %>%
  separate(param,into = c('Method','Param','statistic')) %>%
  pivot_wider(names_from = statistic,
              values_from = stat) ->  ACCURACY_TABLE

STAT_ALL %>%
  select(contains('Kappa')) %>%
  pivot_longer(everything(),values_to = 'stat',names_to = 'param') %>%
  separate(param,into = c('Method','Param','statistic')) %>%
  pivot_wider(names_from = statistic,
              values_from = stat) -> KAPPA_TABLE


# Calculate average value for each replicate using ALL_RESULTS, then use contains and across to create an object with mean values
# This new object is "wide" with 4 columns (REPLICA and the 3 methods): change it with pivot longer
# Now order the mean values (arrange) by method and -ACC, group them by method and keep only the first with (filter(row_number()==1)
# You will have an object with the number of the best replicate (higher accuracy) for each method

ALL_RESULTS %>%
  select(c(REPLICA,contains('Accuracy'))) %>%
  group_by(REPLICA) %>%
  summarise(across(contains('Accuracy'),mean)) -> medie

medie %>%
  pivot_longer(contains('Accuracy'),
               names_to = 'metodo',
               values_to = 'ACC') %>%
  arrange(metodo,-ACC) %>%
  group_by(metodo) %>%
  filter(row_number()==1) %>%
  ungroup() %>%
  separate(metodo,into = c('method','param'),sep='~')-> keep

# Select starting list
finalModel <- NULL
keep %>%
  filter(method==keep$method[1]) %>%
  select(REPLICA) %>%
  pull() -> x1

res.list[[x1]] -> finalModel

# Select the replicates of the other methods and add them to finalModel 

for (i in keep$method[-1]){
  keep %>%
    filter(method==i) %>%
    select(REPLICA) %>%
    pull() -> x1
  print(i)
  
  res.list[[x1]][i]-> A1
  finalModel[i]<- A1
}


# Comparison of final model of different methods
# Summarize accuracy and kappa for each method

results <- resamples(finalModel)
summary(results)

# Plot and compare accuracy of models
pdf(file.path(basefolder,"acc_methods.pdf"), width=7, height=5)
dotplot(results)
dev.off()


# Summarize all models (info about tuning parameters)

print(finalModel)

# See saved predictions on testing

fit.pred <- NULL
for (i in seq(length(finalModel))){
  print(i)
  finalModel[[i]]$pred %>% mutate(metodo=i) -> p
  bind_rows(fit.pred,p) -> fit.pred
}

# Feature importance

# Variable importance on testing
featImp_glm <- varImp(finalModel$glmnet,useModel = FALSE, scale=TRUE)
featImp_glm
#plot(featImp_glm)
pdf(file.path(basefolder,"featImp_glm.pdf"), width=7, height=5)
plot(featImp_glm)
dev.off()

featImp_rf <- varImp(finalModel$rf,useModel = FALSE, scale=TRUE)
featImp_rf
#plot(featImp_rf)
pdf(file.path(basefolder,"featImp_rf.pdf"), width=7, height=5)
plot(featImp_rf)
dev.off()

featImp_svm <- varImp(finalModel$svmRadial,useModel = FALSE, scale=TRUE)
featImp_svm
#plot(featImp_svm)
pdf(file.path(basefolder,"featImp_svm.pdf"), width=7, height=5)
plot(featImp_svm)
dev.off()

featImp_nn <- varImp(finalModel$nnet,useModel = FALSE, scale=TRUE)
featImp_nn
#plot(featImp_nn)
pdf(file.path(basefolder,"featImp_nn.pdf"), width=7, height=5)
plot(featImp_nn)
dev.off()

featImp_glm2 <- varImp(finalModel$glmnet,useModel = TRUE, scale=TRUE)
featImp_glm2
#plot(featImp_glm2)
pdf(file.path(basefolder,"featImp_glm2.pdf"), width=7, height=5)
plot(featImp_glm2)
dev.off()

featImp_rf2 <- varImp(finalModel$rf,useModel = TRUE, scale=TRUE)
featImp_rf2
#plot(featImp_rf2)
pdf(file.path(basefolder,"featImp_rf2.pdf"), width=7, height=5)
plot(featImp_rf2)
dev.off()

featImp_svm2 <- varImp(finalModel$svmRadial,useModel = TRUE, scale=TRUE)
featImp_svm2
#plot(featImp_svm2)
pdf(file.path(basefolder,"featImp_svm2.pdf"), width=7, height=5)
plot(featImp_svm2)
dev.off()

featImp_nn2 <- varImp(finalModel$nnet,useModel = TRUE, scale=TRUE)
featImp_nn2
#plot(featImp_nn2)
pdf(file.path(basefolder,"featImp_nn2.pdf"), width=7, height=5)
plot(featImp_nn2)
dev.off()

# Estimate skill of the models on the test dataset

writeLines(" - Model testing")
## predicted outcome 

pred.val <- NULL
for (i in seq(length(finalModel))){
  print(i)
  as.data.frame(predict(finalModel[[i]],test_reduced))%>% mutate(metodo=i)%>% 
    mutate(pred = predict(finalModel[[i]],test_reduced)) %>% 
    select(pred,metodo) -> pv
  bind_rows(pred.val,pv) -> pred.val
}

## predicted probability (from 0 to 1) 

pred.val2 <- NULL
for (i in seq(length(finalModel))){
  print(i)
  predict(finalModel[[i]],test_reduced, type='prob')%>% mutate(metodo=i) -> pv2
  bind_rows(pred.val2,pv2) -> pred.val2
}


# False positive, false negative and total error across methods
# Matthew's Correlation Coefficient (MCC)


error <- NULL
MCC <- NULL
alg <- unique(pred.val$metodo)

for (i in 1:length(alg)){
  print(i)
  pred.val %>%
    filter(metodo==i) -> pred.val3
  pred.val4 <- cbind(test_reduced, pred.val3$pred)
  
  p.val <- as.matrix(table(pred.val4$Outcome, pred.val4$pred))
  p.val
  d.val <- data.frame(P=sum(p.val[2,]),FN=p.val[2,1],N=sum(p.val[1,]),FP=p.val[1,2], TN=p.val[1,1], TP=p.val[2,2])
  d.val 
  FPR <- (d.val$FP/d.val$N)*100
  FNR <- (d.val$FN/d.val$P)*100
  TER <- (d.val$FN+d.val$FP)/(d.val$P+d.val$N)*100
  df.val <- data.frame(methods= i,error=c("FPR", "FNR", "TER"),
                       percentage=c(FPR, FNR, TER))
  
  d.val$FN <- as.numeric(d.val$FN)
  d.val$FP <- as.numeric(d.val$FP)
  d.val$TN <- as.numeric(d.val$TN)
  d.val$TP <- as.numeric(d.val$TP)
  MCC.val <- (d.val$TP * d.val$TN - d.val$FP * d.val$FN) /
    sqrt ((d.val$TP + d.val$FP) * (d.val$TP + d.val$FN) * (d.val$TN + d.val$FP) * (d.val$TN + d.val$FN))
  MCC.val
  as.data.frame(MCC.val) %>% mutate(metodo=i) -> MCC.val2
  
  bind_rows(error,df.val) -> error
  bind_rows(MCC,MCC.val2) -> MCC
}

# Plot the 3 calculated errors of each method

error %>% 
  mutate(methods=case_when(
    methods == 1 ~ 'RLR',
    methods == 2 ~ 'RF',
    methods == 3 ~ 'SVM',
    TRUE ~ 'NN'))-> error

pdf(file.path(basefolder, "errors.pdf"), width=7, height=2.5)
ggplot(error, aes(fill=error, y=percentage, x=methods)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))
dev.off()

png(file.path(basefolder,"errors.png"))
ggplot(error, aes(fill=error, y=percentage, x=methods)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))
dev.off()

# ROC curves 

ROC <- list()
alg <- unique(pred.val$metodo)
alg2 <- c('GLM', 'RF', 'SVM', 'NN')
ind <- 1 #index to access the methods name in the metodi vector.

pdf(file.path(basefolder,"roc.pdf"))
par(mfrow=c(2,2))
par(pty = "s")
for (i in 1:length(alg)){
  print(i)
  pred.val2 %>% 
    filter(metodo==i) -> pred.val2_metodo
  met <- alg2[ind]
  ind <- ind+1
  roc <- roc(test_reduced$Outcome, pred.val2_metodo$Bacteria, plot=TRUE, legacy.axes=TRUE, 
             percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", 
             col="#377eb8", lwd=2, print.auc = TRUE, ci = TRUE, print.auc.y = 10, 
             print.auc.x = 80, main= paste0("ROC curve "," ", met))
  
  ROC[[i]] <- roc 
}
dev.off()
par(pty = "m")


png(file.path(basefolder,"roc.png"))
par(mfrow=c(2,2))
par(pty = "s")
for (i in 1:length(alg)){
  print(i)
  pred.val2 %>% 
    filter(metodo==i) -> pred.val2_metodo
  met <- alg2[ind]
  ind <- ind+1
  roc <- roc(test_reduced$Outcome, pred.val2_metodo$Bacteria, plot=TRUE, legacy.axes=TRUE, 
             percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", 
             col="#377eb8", lwd=2, print.auc = TRUE, ci = TRUE, print.auc.y = 10, 
             print.auc.x = 80, main= paste0("ROC curve "," ", met))
  
  ROC[[i]] <- roc 
}
dev.off()
par(pty = "m")




#########################################################################

#ggplot(pred.val2, aes(y=Bacteria, x=as.factor(metodo))) + 
#  geom_boxplot(aes(fill=as.factor(metodo))) 
  

df <- data.frame(test$ID,test_reduced$Outcome, pred.val$pred)


############################################################################


# Confusion matrix e statistics (sensitivity, specificity...) validation set

CM <- list()
alg <- unique(pred.val$metodo)

for (i in 1:length(alg)){
  print(i)
  pred.val %>% 
    filter(metodo==i) -> pred.val_metodo
  cm <- confusionMatrix(pred.val_metodo$pred, test_reduced$Outcome, positive = "Bacteria", mode = "everything")
  CM[[i]] <- cm 
}

CM

# Save workspace, MCC and errors

write.csv(MCC, file = file.path(basefolder,"MCC.csv"), row.names=FALSE)
write.csv2(error, file = file.path(basefolder,"error.csv"), row.names=FALSE)
save.image(file.path(basefolder,"ML_archaea.RData"))

write.csv2(df, file = file.path(basefolder,"df_2.csv"), row.names=FALSE)


