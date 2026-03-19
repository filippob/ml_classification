
library("themis")
library("kernlab")
library("yardstick")
library("tidymodels")
library("data.table")

## Parameters
basefolder <- "/home/filippo/Documents/tania/probiotics"
tuned_model <- "results/oneclass_tuned_model.RData"
train_set = "splits/train_set.csv"
test_set = "splits/filtered_merged_bits26_testset.csv"
outdir = "results"
nproc <- 4
id_vars = c("Organism", "Taxon", "Definition")
normal_class <- "Probiotic"
target_var = "Label"
flag_manual = TRUE ## for manual explicit workflow
flag_evaluation = FALSE

## load tuned model
writeLines(" - loading tuned model")
fname = file.path(basefolder, tuned_model)
load(fname)

fine_tune_res <- to_save[[1]]
svm_recipe = to_save[[2]]

## 6) select best hyperparameters
best_params <- fine_tune_res %>%
  arrange(desc(score)) %>%
  slice(1)

print("best model")
print(best_params)

## Import training data
writeLines(" - reading the training dataset")
fname <- file.path(basefolder, train_set)
ext = gsub("^.*\\.","",fname)
if (ext %in% c("xlsx")) {
  
  dataset <- readxl::read_xlsx(fname)
} else dataset <- fread(fname)

print(paste("Training data size (n. records):",nrow(dataset)))

## Data preparation
writeLines(" - preparing the training data")
dataset <- dataset %>% 
  mutate({{target_var}} := as.factor(.data[[target_var]]))

dt_train_normal <- subset(dataset, Label == normal_class)
train_normal_class_label <- dataset$Label

additional_test <- subset(dataset, Label != normal_class)
test_label <- additional_test$Label

dt_train_normal$Label <- NULL
additional_test$Label <- NULL

dt_train_normal <- dt_train_normal |> select(!all_of(id_vars))
additional_test_reduced <- additional_test |> select(!all_of(id_vars))

prep_dt <- prep(svm_recipe)
print(prep_dt)

training_set <- juice(prep_dt)

## 7) train final model
writeLines(" - training the oneclass SVM model")
final_model <- ksvm(
  x = as.matrix(training_set),
  type = "one-svc",
  kernel = rbfdot(sigma = best_params$sigma),
  nu = best_params$nu,
  scaled = TRUE
)

## test data
writeLines(" - read test data")
fname = file.path(basefolder, test_set)
test <- fread(fname)

test_name = basename(test_set)

## recipe (preprocessing)
test_set <- bake(prep_dt, new_data = test)

#####################################################
#####################################################
## 8) evaluate trained model
#####################################################
## Decision value	Meaning
## > 0	inside the learned region (inlier / normal)
## ≈ 0	near the boundary
## < 0	outside the region (outlier / anomaly)
#####################################################

preds <- predict(final_model, test_set, type = "decision")
summary(preds)
test$decision = preds

vec <- which(preds < 0)

if(flag_evaluation) {
  
  writeLines(" - model evaluation: decision")
  ndec = ncol(test)
  nonprobs <- test[vec, c(1,2,3,4,ndec), with = FALSE]
  nonprobs$predicted_label = "Nonprobiotic"
  
  probs <- test[-vec, c(1,2,3,4,ndec), with = FALSE]
  probs$predicted_label = "Probiotic"
  
  preds <- probs |> bind_rows(nonprobs)
  preds$Label = as.factor(preds$Label)
  preds$predicted_label = as.factor(preds$predicted_label)
  
  # Large positive values → strongly normal
  # Values near 0 → borderline
  # Large negative values → strong anomalies
  
  print("Confusion Matrix")
  cm <- preds |>
    conf_mat(!!target_var, predicted_label)
  print(cm)
  
  fname <- file.path(basefolder, outdir, "oneclass-confusion_matrix.png")
  g <- autoplot(cm, type = "heatmap") + scale_fill_gradientn(colours = c("lightyellow", "yellow", "orange", "red"))
  print(g)
  ggsave(filename = fname, plot = g, device = "png", width = 5, height = 4.5)
  
  print("Performance metrics")
  # auc = roc_auc(preds, truth = Label, .pred_Nonprobiotic)
  acc = accuracy(preds, truth = Label, estimate = predicted_label)
  mcc = mcc(preds, truth = Label, estimate = predicted_label)
  # brier_score = brier_class(preds, truth = Label, .pred_Probiotic)
  
  df_metrics <- bind_rows(acc,mcc)
  print(df_metrics)
}

#####################################################
#####################################################

## error analysis
writeLines(" - error analysis")
if(flag_evaluation) {
  
  vec <- which(preds$Label != preds$predicted_label)
  errors <- preds[vec,]
  print(paste("N. of mismatches between labels from the test set and from collected predictions", length(vec)))
} else {
  
  errors <- test[-vec,] |> select(c(all_of(id_vars), "decision"))
}


writeLines(" - saving results")
temp = paste("oneclass-all_predictions-", file_path_sans_ext(test_name), ".csv", sep="")
fname = file.path(basefolder, outdir, temp)
fwrite(x = preds, file = fname, sep = "\t")

temp = paste("oneclass-errors-", file_path_sans_ext(test_name), ".csv", sep="")
fname = file.path(basefolder, outdir, temp)
fwrite(x = errors, file = fname, sep = "\t")

print("DONE!!")