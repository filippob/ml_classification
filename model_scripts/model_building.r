# ============================================================================
## script for model building and selection through cross-validation 
## (no test set)
# ============================================================================

library("dplyr")
library("tidyverse")
library("tidymodels")
library("data.table")
library("randomForest") 


###################
# script parameters
###################
base_folder <- "~/Documents/tania/vertrebrates_invertebrates/"
infile = "cleaned_kmer_frequencies.csv"
outdir = "results"

## preprocessing
collinearity_threshold = 0.95
nproc <- 4

## cross-validation
k_folds <- 5
nrepeats_cv <- 2

## random forest
mtry_range = c(3,8) ## n. of variables to be samples randomly in each tree
min_n_range = c(3,8) ## min number of obs in terminal nodes
ntrees = 200
class_weights = TRUE

# create output folder
fname = file.path(base_folder, outdir)
if (!dir.exists(fname)) {
  dir.create(fname, recursive = TRUE, showWarnings = FALSE)
}

## input data
writeLines(" - reading input data")

fname = file.path(base_folder, infile)
inpdata <- fread(fname)


## DESCRIPTION
## average value of numerical variables
print("average values of features")
inpdata |>
  summarise_if(is.numeric, mean, na.rm = TRUE)

print("standard deviations of features")
inpdata |>
  summarise_if(is.numeric, sd, na.rm = TRUE)


# ============================================================================
# 4. PREPROCESSING RECIPE ()
# ============================================================================
## residual preprocessing to avoid data quirks specific to individual data partitions
## e.g. zero-variance or correlated data in data subsets
writeLines(" - residual preprocessing ... ")

preprocess_recipe <- recipe(Outcome ~ ., data = inpdata) %>%
  update_role(Species, new_role = "Species") %>%  # Species will not be used as predictor
  step_zv(all_numeric(), -all_outcomes()) %>%
  step_corr(all_predictors(), threshold = collinearity_threshold)


# ============================================================================
# 5. MODEL SPECIFICATION : rf/lasso/etc.
# ============================================================================
writeLines(" - model building and selection ...")

# Calculate weights
if (class_weights) {
  
  class_freq <- table(inpdata$Outcome)
  class_weights <- sum(class_freq) / (length(class_freq) * class_freq)
  weight_vector = c(as.numeric(class_weights["invertebrate"]), as.numeric(class_weights["vertebrate"]))
  
} else weight_vector =c(1,1)

print("class weights")
print(weight_vector)


rf_spec <- rand_forest(
  mtry = tune(),
  trees = ntrees,
  min_n = tune()
) %>%
  set_mode("classification") %>%
  set_engine("randomForest", 
             classwt = weight_vector)

print("model specification")
print(rf_spec)

tune_rf_wf <- workflow() %>%
  add_recipe(preprocess_recipe) %>%   
  add_model(rf_spec)

print("tuning workflow")
print(tune_rf_wf)

writeLines(" - CV for fine-tuning of the hyperparameters ..")
cv_folds <- vfold_cv(inpdata, v = k_folds, repeats = nrepeats_cv, strata = Outcome)  # (GROUP_VFOLD-CV BY ORDINE?)


rf_grid <- grid_regular(
  mtry(range = mtry_range),
  min_n(range = min_n_range),
  levels = c(max(mtry_range),max(min_n_range))
)

rf_fine_tune_res <- tune_grid(
  tune_rf_wf,
  metrics = metric_set(roc_auc, accuracy, mcc),
  resamples = cv_folds,
  grid = rf_grid
)

autoplot(rf_fine_tune_res)

# Plot tuning results of mcc

p1 <- rf_fine_tune_res %>%
  collect_metrics() %>%
  filter(.metric == "mcc") %>%
  ggplot(aes(x = mtry, y = mean, color = factor(min_n), group = min_n)) +
  geom_line(alpha = 0.5, size = 1.5) +
  geom_point(size = 2) +
  labs(y = "MCC", color = "min_n", title = "MCC vs mtry") +
  theme_minimal()

print(p1)

# Collect all metrics
all_tuning_metrics <- collect_metrics(rf_fine_tune_res)

## Save all hyperparameters
writeLines(" - saving results from model tuning")
fname = file.path(base_folder, outdir, "all_tuning_metrics.csv")
fwrite(x = all_tuning_metrics, file = fname, sep = ",", col.names = TRUE)

# Select best hyperparameters
writeLines(" - selecte best hyperparameters")
show_best(rf_fine_tune_res, metric = "mcc")
rf_best_params <- select_best(rf_fine_tune_res, metric = "mcc")
print(rf_best_params)

# Extrrf_best_params# Extract best hyperparameters
best_mtry <- rf_best_params$mtry
best_min_n <- rf_best_params$min_n

# Filter by best penalty and collect the metrics
best_penalty_metrics <- all_tuning_metrics %>%
  filter(mtry == best_mtry, min_n == best_min_n) %>%
  select(.metric, mean, std_err) %>%
  pivot_wider(names_from = .metric, values_from = c(mean, std_err), names_sep = "_")

# Create a full dataframe
best_model <- data.frame(
  mtry = best_mtry,
  min_n = best_min_n,
  mcc_mean = best_penalty_metrics$mean_mcc,
  mcc_std_err = best_penalty_metrics$std_err_mcc,
  roc_auc_mean = best_penalty_metrics$mean_roc_auc,
  roc_auc_std_err = best_penalty_metrics$std_err_roc_auc,
  accuracy_mean = best_penalty_metrics$mean_accuracy,
  accuracy_std_err = best_penalty_metrics$std_err_accuracy
)

## Save best hyperparameters
writeLines(" - saving best model parameters")
fname = file.path(base_folder, outdir, "best_hyperparameters_tuning_rf.csv")
fwrite(x = best_model, file = fname, sep = ",", col.names = TRUE)

print("DONE!")

