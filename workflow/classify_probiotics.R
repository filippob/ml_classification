
# Load libraries
library("pROC")
library("caret")
library("glmnet")
library("ggplot2")
library("tidyverse")
library("doParallel")
library("data.table")
library("caretEnsemble")



######################################################
## parameters
basefolder = "/home/biscarinif/probiotics"
ncpu = 8
seed = 7
seed_rfe = 325
training_pct = 0.85
nfolds = 10
nrepeats = 100

# Remove warnings
options(warn=-1)

# Setup parallel back end to use multiple processors
cl = makeCluster(ncpu)
registerDoParallel(cl)
######################################################

# Import dataset
fname = file.path(basefolder, "batteri_filtered_TB.csv")
dataset <- fread(fname)

dataset <- dataset %>% 
  mutate(Outcome = as.factor(Label)) %>% 
  select(!c(Label,Class, Taxon, Locus_ID, Definition)) %>% 
  select(Outcome, everything())

group_by(dataset, Outcome) |>
  summarise(N = n()) |>
  print()


###############################
## Machine learning analysis ##
###############################
# Creation of train/test and validation datasets
# We will split the dataset into two, 85% of which will be used to train and test our models 
# and 15% that we will be hold back as a validation dataset

dataset <- data.frame(dataset) # to solve problems with rfe function

# set.seed(seed)
validation_index <- createDataPartition(dataset$Outcome, p=training_pct, list=FALSE)

# Select 15% of the data for validation
validation <- dataset[-validation_index,]

# Use the remaining 85% of data to train and test the models
training <- dataset[validation_index,]

validation %>% 
  count(Outcome) 

training %>% 
  count(Outcome) 

##################################################################################
# Feature selection (rfe): automatically selecting a subset of the most predictive
# features on train/test dataset
##################################################################################

# 10-fold cross validation (cv) repeated 100 times
# Use different seed for cv in rfe and train

writeLines(" - RFE: recursive feature elimination")
# set.seed(seed_rfe)
rfectrl_rf <- rfeControl(functions=rfFuncs, method="repeatedcv", number=nfolds, repeats=nrepeats, 
                         verbose = FALSE, returnResamp = "final", allowParallel = TRUE) 

# After indicating predictors(3-170) and outcome (1), excluding id (2),
# check from 1 to n possible combinations of predictors to define the best subset 
results_rfe <- rfe(training[,3:170], training[,1], sizes=c(1:168), rfeControl=rfectrl_rf)
print(results_rfe$fit)

# Plot the results

# plot(results_rfe, type=c("g", "o"))

fname = file.path(basefolder, "rfe.pdf")
pdf(fname, width=7, height=5)
plot(results_rfe, type=c("g", "o"))
dev.off()

fname = file.path(basefolder, "rfe.png")
png(filename=fname)
plot(results_rfe, type=c("g", "o"))
dev.off()


# List the chosen features
features = predictors(results_rfe)
print(features)
fname = file.path(basefolder, "features.RData")
save(features, file = fname)

# Prepare datasets (train/test= dataset3 and validation), excluding non-selected features
training_rfe <- training[,c("Outcome", "Organism", features)] # keep only selected features from dataset3 
training_rfe <- as.data.frame(training_rfe) # save as dataframe
validation_rfe <- validation[,c("Outcome", "Organism", features)] # keep only selected features from validation 
validation_rfe <- as.data.frame(validation_rfe) #save as dataframe

# Set up test harness, creation and evaluation of different models to predict outcome 

# Set up test harness, creation and evaluation of different models to predict outcome 
# Run algorithms using stratified 10-fold cross validation, repeated 100 times (repeated k-fold cv)
# Stratification will be used in case of imbalanced output classes
# Data are centered and scaled within cv
# Exclude ID from dataset_rfe:[,-2]

my_list <- list()

for(i in 1:nrepeats) {
  
  folds <- createMultiFolds(training_rfe$Outcome, k = nfolds, times = 1) #stratification
  
  trainctrl <- trainControl(method="repeatedcv", number=nfolds, repeats=1, index = folds,
                          classProbs= TRUE, savePredictions = TRUE, allowParallel = TRUE, sampling = "up") 
  metric = "Accuracy" 
  metodi = c('glmnet', 'rf', 'svmRadial', 'nnet')
  
  # set.seed(i)
  fit.models <- caretList(Outcome ~ ., data=training_rfe[,-2], trControl=trainctrl, preProcess = c("center", "scale"),
                          methodList=metodi, metric=metric) 
  
  my_list[[i]] <- fit.models 
  
}

# Create a data.frame with accuracy and Kappa of the 10 fold of all replicates for each method

ALL_RESULTS = NULL

for (i in seq(length(my_list))){
  print(i)
  summary(resamples(my_list[[i]]))[[1]] %>% mutate(REPLICA=i)-> X
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
              values_from = stat)-> ACCURACY_TABLE

STAT_ALL %>%
  select(contains('Kappa')) %>%
  pivot_longer(everything(),values_to = 'stat',names_to = 'param') %>%
  separate(param,into = c('Method','Param','statistic')) %>%
  pivot_wider(names_from = statistic,
              values_from = stat)-> KAPPA_TABLE


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

my_list[[x1]] -> finalModel

# Select the replicates of the other methods and add them to finalModel 

for (i in keep$method[-1]){
  keep %>%
    filter(method==i) %>%
    select(REPLICA) %>%
    pull() -> x1
  print(i)
  
  my_list[[x1]][i]-> A1
  finalModel[i]<- A1
}


# Comparison of final model of different methods
# Summarize accuracy and kappa for each method

results <- resamples(finalModel)
summary(results)

# Plot and compare accuracy of models

fname = file.path(basefolder, "acc_methods.pdf")
pdf(fname, width=7, height=5)
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

# set.seed(7)
featImp_glm <- varImp(finalModel$glmnet,useModel = FALSE, scale=TRUE)
print(featImp_glm)

fname = file.path(basefolder, "featImp_glm.pdf")
pdf(fname, width=7, height=5)
plot(featImp_glm)
dev.off()

# set.seed(7)
featImp_rf <- varImp(finalModel$rf,useModel = FALSE, scale=TRUE)
print(featImp_rf)

fname = file.path(basefolder, "featImp_rf.pdf")
pdf(fname, width=7, height=5)
plot(featImp_rf)
dev.off()

# set.seed(7)
featImp_svm <- varImp(finalModel$svmRadial,useModel = FALSE, scale=TRUE)
print(featImp_svm)

fname = file.path(basefolder, "featImp_svm.pdf")
pdf(fname, width=7, height=5)
plot(featImp_svm)
dev.off()

# set.seed(7)
featImp_nn <- varImp(finalModel$nnet,useModel = FALSE, scale=TRUE)
print(featImp_nn)

fname = file.path(basefolder, "featImp_nn.pdf")
pdf(fname, width=7, height=5)
plot(featImp_nn)
dev.off()

# set.seed(7)
featImp_glm2 <- varImp(finalModel$glmnet,useModel = TRUE, scale=TRUE)
print(featImp_glm2)

fname = file.path(basefolder, "featImp_glm2.pdf")
pdf(fname, width=7, height=5)
plot(featImp_glm2)
dev.off()

# set.seed(7)
featImp_rf2 <- varImp(finalModel$rf,useModel = TRUE, scale=TRUE)
print(featImp_rf2)

fname = file.path(basefolder, "featImp_rf2.pdf")
pdf(fname, width=7, height=5)
plot(featImp_rf2)
dev.off()

# set.seed(7)
featImp_svm2 <- varImp(finalModel$svmRadial,useModel = TRUE, scale=TRUE)
print(featImp_svm2)

fname = file.path(basefolder, "featImp_svm2.pdf")
pdf(fname, width=7, height=5)
plot(featImp_svm2)
dev.off()

# set.seed(7)
featImp_nn2 <- varImp(finalModel$nnet,useModel = TRUE, scale=TRUE)
print(featImp_nn2)

fname = file.path(basefolder, "featImp_nn2.pdf")
pdf(fname, width=7, height=5)
plot(featImp_nn2)
dev.off()

# Estimate model performance on the validation dataset
## predicted outcome 
pred.val <- NULL
for (i in seq(length(finalModel))){
  print(i)
  as.data.frame(predict(finalModel[[i]],validation_rfe))%>% mutate(metodo=i)%>% 
    mutate(pred = predict(finalModel[[i]],validation_rfe)) %>% 
    select(pred,metodo) -> pv
  bind_rows(pred.val,pv) -> pred.val
}

## predicted probability (from 0 to 1) 
prob.val <- NULL
for (i in seq(length(finalModel))){
  print(i)
  predict(finalModel[[i]],validation_rfe, type='prob')%>% mutate(metodo=i) -> pv2
  bind_rows(prob.val,pv2) -> prob.val
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
  pred.val4 <- cbind(validation_rfe, pred.val3$pred)
  
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

fname = file.path(basefolder, "errors.pdf")
pdf(fname, width=7, height=2.5)
ggplot(error, aes(fill=error, y=percentage, x=methods)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))
dev.off()

fname = file.path(basefolder, "errors.png")
png(fname)
ggplot(error, aes(fill=error, y=percentage, x=methods)) + 
  geom_bar(position="dodge", stat="identity") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))
dev.off()

fname = file.path(basefolder, "errors.csv")
fwrite(x = error, file = fname)

fname = file.path(basefolder, "mcc.csv")
fwrite(x = MCC, file = fname)

# ROC curves 
ROC <- list()
alg <- unique(pred.val$metodo)
alg2 <- c('RLR', 'RF', 'SVM', 'NN')
ind <- 1 #index to access the methods name in the metodi vector.

fname = file.path(basefolder, "roc.pdf")
pdf(fname)
par(mfrow=c(2,2))
par(pty = "s")
for (i in 1:length(alg)){
  print(i)
  prob.val %>% 
    filter(metodo==i) -> prob.val_metodo
  met <- alg2[ind]
  ind <- ind+1
  roc <- roc(validation_rfe$Outcome, prob.val_metodo$Probiotico, plot=TRUE, legacy.axes=TRUE, 
             percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", 
             col="#377eb8", lwd=2, print.auc = TRUE, ci = TRUE, print.auc.y = 10, 
             print.auc.x = 80, main= paste0("ROC curve "," ", met))
  
  ROC[[i]] <- roc 
}
dev.off()
par(pty = "m")


fname = file.path(basefolder, "roc.png")
png(fname)
par(mfrow=c(2,2))
par(pty = "s")
for (i in 1:length(alg)){
  print(i)
  prob.val %>% 
    filter(metodo==i) -> prob.val_metodo
  met <- alg2[ind]
  ind <- ind+1
  roc <- roc(validation_rfe$Outcome, prob.val_metodo$Probiotico, plot=TRUE, legacy.axes=TRUE, 
             percent=TRUE, xlab="False Positive Percentage", ylab="True Postive Percentage", 
             col="#377eb8", lwd=2, print.auc = TRUE, ci = TRUE, print.auc.y = 10, 
             print.auc.x = 80, main= paste0("ROC curve "," ", met))  
  
  ROC[[i]] <- roc 
}
dev.off()
par(pty = "m")



# Confusion matrix e statistics (sensitivity, specificity...) validation set
CM <- list()
alg <- unique(pred.val$metodo)

for (i in 1:length(alg)){
  print(i)
  pred.val %>% 
    filter(metodo==i) -> pred.val_metodo
  cm <- confusionMatrix(pred.val_metodo$pred, validation_rfe$Outcome, positive = "Probiotico", mode = "everything")
  CM[[i]] <- cm 
}

fname = file.path(basefolder, "CM.RData")
save(CM, file = fname)

# Save workspace, MCC and errors
# save.image("prob_bits.RData")

print("DONE!")

