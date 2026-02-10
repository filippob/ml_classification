## script to evaluate model through cross-validation

# library("yaml")
library("dplyr")
library("tidyverse")
library("tidymodels")
library("data.table")
# library("tidypredict")

###################
# script parameters
###################
base_folder <- "~/Documents/tania/vertrebrates_invertebrates/"
infile = "cleaned_kmer_frequencies.csv"
paramfile = "best_hyperparameters_tuning_lasso.csv"
outdir = "results"

## model
## cross-validation
k_folds <- 5
nrepeats_cv <- 10
method = "lasso" ## "rf", "lasso"
class_weights = TRUE

## preprocessing
collinearity_threshold = 0.95
nproc <- 4

# create output folder
fname = file.path(base_folder, outdir)
if (!dir.exists(fname)) {
  dir.create(fname, recursive = TRUE, showWarnings = FALSE)
}

## input data
writeLines(" - reading input data")

fname = file.path(base_folder, infile)
inpdata <- fread(fname)

## model parameters from fine-tuning
fname = file.path(base_folder, outdir, paramfile)
params = fread(fname)
print("Loaded model parameters (from fine-tuning):")
print(params)


## 
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

if (method == "lasso") {
  
  mod_spec <- logistic_reg(
    penalty = params$lambda, 
    mixture = 1 ## 1: lasso; 0: ridge
  ) |>
    set_mode("classification") |>
    set_engine("glmnet",
               classwt = weight_vector)
} else if (method == "rf") {
  
  mod_spec <- rand_forest(
    mtry = params$mtry,
    trees = params$ntrees,
    min_n = params$min_n
  ) %>%
    set_mode("classification") %>%
    set_engine("randomForest", 
               classwt = weight_vector)
}

print("model specification")
print(mod_spec)

model_wf <- workflow() %>%
  add_recipe(preprocess_recipe) %>%   
  add_model(mod_spec)

print("workflow")
print(model_wf)

#####
writeLines(" - CV for fine-tuning of the hyperparameters ..")
cv_folds <- vfold_cv(inpdata, v = k_folds, repeats = nrepeats_cv, strata = Outcome)  # (GROUP_VFOLD-CV BY ORDINE?)

writeLines(" - fitting model to resampled data")
model_fit_cv <- 
  model_wf %>% 
  fit_resamples(cv_folds, metrics = metric_set(roc_auc, accuracy, mcc))

all_performances <- collect_metrics(model_fit_cv, summarize = FALSE)
print("summary of model performance")
collect_metrics(model_fit_cv, summarize = TRUE)

p <- ggplot(all_performances, aes(x = .metric, y = .estimate)) + 
  geom_jitter(aes(color = .metric), width = 0.3)

