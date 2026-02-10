
library("dplyr")
library("yardstick")
library("tidyverse")
library("tidymodels")
library("data.table")


###################
# script parameters
###################
base_folder <- "~/Documents/tania/vertrebrates_invertebrates/"
test_data = "cleaned_kmer_frequencies.csv"
trained_model = "trained_model_lasso.RDS"
outdir = "results"
method = "lasso" ## "rf", "lasso"
class_weights = TRUE

## preprocessing
collinearity_threshold = 0.95
nproc <- 4

### write out chosen parameters
print("###################################")
print(paste("Trained model file is:", trained_model))
print(paste("Method is :", method))
temp = ifelse(class_weights, "Yes", "No")
print(paste("Class weights:", temp))
print("###################################")

# create output folder
fname = file.path(base_folder, outdir)
if (!dir.exists(fname)) {
  dir.create(fname, recursive = TRUE, showWarnings = FALSE)
}

## input data
writeLines(" - reading input data")
fname = file.path(base_folder, test_data)
test_data <- fread(fname)

## preprocess test data
writeLines(" - preprocessing test data ")

preprocess_recipe <- recipe(Outcome ~ ., data = test_data) %>%
  update_role(Species, new_role = "Species") %>%  # Species will not be used as predictor
  step_zv(all_numeric(), -all_outcomes()) %>%
  step_corr(all_predictors(), threshold = collinearity_threshold)

prepped_recipe <- prep(preprocess_recipe)
testdata_preprocessed <- bake(prepped_recipe, new_data = NULL)

## trained model
writeLines(" - reading trained model")
fname = file.path(base_folder, outdir, trained_model)
trained_model <- readRDS(fname)
# print("Loaded model (from training):")
# print(trained_model)

## predictions
writeLines(" - calculating predictions")

preds <- bind_cols(
  predict(trained_model, testdata_preprocessed),
  predict(trained_model, testdata_preprocessed, type = "prob")
)

temp <- testdata_preprocessed |>
  select(Species, Outcome) |>
  mutate(Outcome = as.factor(Outcome))

preds <- bind_cols(temp, preds)
preds <- preds %>%
  mutate(.pred_class = factor(.pred_class, levels = levels(Outcome)))

model_acc = accuracy(preds, truth = Outcome, estimate = .pred_class)
model_roc = roc_auc(preds, truth = Outcome, .pred_invertebrate)
model_mcc = mcc(preds, truth = Outcome, .pred_class)

model_performance <- bind_rows(model_acc, model_roc, model_mcc)
print("Model performance on test data")
print(model_performance)

writeLines(" - saving model performance")
temp = paste("test_metrics_", method, ".csv", sep = "")
fname = file.path(base_folder, outdir, temp)
fwrite(x = model_performance, file = fname, sep = ",", col.names = TRUE)

writeLines(" - saving predictions")
temp = paste("test_predictions_", method, ".csv", sep = "")
fname = file.path(base_folder, outdir, temp)
fwrite(x = preds, file = fname, sep = ",", col.names = TRUE)

print("DONE!!")